"""Headroom scan for stage 1 of direction one: where is a short-horizon certified MPC clearly suboptimal?

Learned terminal values can only help where the certified parallel MPC is infeasible or noticeably
worse than J*. The examples of the Parallel MPC paper leave <= 0.9% (pr_mpc_repro.py), so this script
measures the headroom of candidate stage-1 settings before any learning is added.

Candidates here: the paper's Examples 4.1 and 4.2 with the paper's certified units, run with shorter
lookaheads (exact QPs). On a grid of initial states in the box |x|_inf <= 5:
  coverage = share of states with finite J* where the certified parallel MPC is feasible,
  gap      = closed-loop cost / J* - 1 on the covered states.
J* is the converged lower bound (T^N 0)(x0) (QP for 4.1, MIQP for 4.2).

The pendulum swing-up candidate is in pendulum_headroom.py.
"""
import argparse

import numpy as np

import pr_mpc_repro as pr


def grid_states(n):
    g = np.linspace(-pr.XMAX, pr.XMAX, n)
    return [np.array([a, b]) for a in g for b in g]


def summarize(name, jstar, covered, costs):
    fin = np.isfinite(jstar)
    cov = np.sum(fin & covered) / np.sum(fin)
    pos = fin & covered & (jstar > 1e-9)
    gap = costs[pos] / jstar[pos] - 1
    if not len(gap):
        print(f"  {name:<34} coverage {cov:6.1%}")
        return
    print(f"  {name:<34} coverage {cov:6.1%}   gap mean {np.mean(gap):6.2%}  median {np.median(gap):6.2%}  "
          f"max {np.max(gap):6.2%}  (> 5% on {np.mean(gap > 0.05):5.1%} of covered states)")


def linear_part(n_grid, T=60):
    print("=" * 100)
    print(f"paper examples with shorter lookaheads, {n_grid}x{n_grid} grid of x0, closed loop T = {T}")

    # Example 4.1 (double integrator), the four certified units of the paper (polytopic terminal sets)
    units = pr.lq_units()
    xs = grid_states(n_grid)
    jstar = np.array([pr.lq_lower_bound(x) for x in xs])
    print(f"\n  Example 4.1 (LQ): J* finite on {np.mean(np.isfinite(jstar)):.1%} of the box")
    configs = [((ell,) * 4, f"l = {ell}, LQR unit only", [0]) for ell in (1, 2, 3)]
    configs += [(ells, f"l = {ells}, all 4 units", [0, 1, 2, 3])
                for ells in [(1, 1, 1, 1), (2, 2, 2, 2), (2, 2, 3, 3), (3, 3, 3, 3)]]
    for ells, label, idx in configs:
        probs = [pr.Lookahead([(pr.A_LQ, pr.B_LQ)] * ells[i], 1.0, units[i]["P"], term=units[i]["S"])
                 for i in idx]

        def policy(x, probs=probs):
            return pr.best(probs, x)[:2]

        covered = np.array([np.isfinite(pr.best(probs, x)[0]) for x in xs])
        costs = np.array([pr.lq_closed_loop(policy, x, T)[0] if c else np.inf for x, c in zip(xs, covered)])
        summarize(label, jstar, covered, costs)

    # Example 4.2 (PWA), the four certified units of Table I
    pw = pr.pwa_units()
    jstar = np.array([pr.pwa_lower_bound(x, 15) for x in xs])
    print(f"\n  Example 4.2 (PWA): J* finite on {np.mean(np.isfinite(jstar)):.1%} of the box")
    for ell in (1, 2):
        probs = {n: pr.pwa_lookaheads(pw[n]["P"], pw[n]["S"], ell) for n in pr.UNITS_PWA}
        for label, names in [("zero unit only", ["zero"]), ("LQR unit only", ["LQR"]), ("all 4 units", pr.UNITS_PWA)]:
            def policy(x, names=names, probs=probs):
                vals = [pr.pwa_unit_value(probs[n], x) for n in names]
                return vals[int(np.argmin([v for v, _ in vals]))]

            covered = np.array([np.isfinite(policy(x)[0]) for x in xs])
            costs = np.array([pr.pwa_closed_loop(policy, x, T) if c else np.inf for x, c in zip(xs, covered)])
            summarize(f"l = {ell}, {label}", jstar, covered, costs)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", type=int, default=11, help="grid points per axis")
    args = ap.parse_args()
    linear_part(args.grid)
