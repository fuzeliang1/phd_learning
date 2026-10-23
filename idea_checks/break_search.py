"""Do certificate breaks actually happen in certified_parallel_mpc.py, and do they matter?

A learned unit E (lookahead 2) is added to certified units on the scalar example, with E's
terminal value taken from one of several error families built from J*:

  hole      J*(x) * (1 - a exp(-((x - c)/w)^2))      spurious minimum (like unit C)
  sin       J*(x) * (1 + e sin(k x + p))              smooth relative error (like unit D)
  deadzone  max(0, J*(x) - b)                        value is 0 near the goal, offset b elsewhere
  band      J*(x) - d on r1 <= |x| <= r2 (smoothed)   locally optimistic offset
  scaled    s J*(x)                                   uniformly optimistic

For each alpha and certified set it reports, per x0 in a grid:
  genuine breaks   v_t > v_{t-1} - alpha g_{t-1} + BREAK_TOL with |x_t| > 1e-3,
  lock vs no lock  closed-loop cost with and without locking after a break,
  certificates     whether the cost exceeds v_0/alpha (Theorem 1, valid only without breaks)
                   and whether it exceeds the a posteriori bound of Theorem 2,
and, once per family, the cost of the naive minimum over the same units without any test.
Runtime about 30 minutes on one core.
"""
import itertools
import time

import numpy as np

import certified_parallel_mpc as cpm

X0S = np.linspace(-0.8, 0.8, 9)
T = 30
js = cpm.jstar


def families():
    for a, c, w in itertools.product([0.5, 0.8, 0.95], [-0.6, -0.3, 0.3, 0.6], [0.08, 0.2, 0.4]):
        yield "hole", f"a={a} c={c} w={w}", lambda x, a=a, c=c, w=w: js(x) * (1 - a * np.exp(-((np.asarray(x) - c) / w) ** 2))
    for e, k, p in itertools.product([0.3, 0.6], [3, 10, 25], [0.0, np.pi / 2, np.pi]):
        yield "sin", f"e={e} k={k} p={p:.2f}", lambda x, e=e, k=k, p=p: js(x) * (1 + e * np.sin(k * np.asarray(x) + p))
    for b in [0.01, 0.03, 0.1, 0.3, 0.6, 1.0]:
        yield "deadzone", f"b={b}", lambda x, b=b: np.maximum(0.0, js(x) - b)
    for d, r1, r2 in [(0.1, 0.1, 0.4), (0.2, 0.15, 0.5), (0.3, 0.2, 0.7), (0.15, 0.05, 0.3)]:
        def band(x, d=d, r1=r1, r2=r2):
            ax = np.abs(np.asarray(x))
            window = 1 / ((1 + np.exp(-(ax - r1) / 0.01)) * (1 + np.exp((ax - r2) / 0.01)))
            return np.maximum(0.0, js(x) - d * window)
        yield "band", f"d={d} [{r1},{r2}]", band
    for s in [0.3, 0.6]:
        yield "scaled", f"s={s}", lambda x, s=s: s * js(x)


def ratio(cost, x0):
    j = float(js(x0))
    return cost / j if j > 1e-9 else (1.0 if cost < 1e-9 else np.inf)


def certified_value(x, cert):
    """V_F(x): smallest predicted cost over the certified units (used in Theorem 2)."""
    return min(float(cpm.solve_unit(x, cpm.UNITS[n][0], cpm.UNITS[n][1])[0]) for n in cert)


def run_cert(cert, x0, lock):
    log = []
    cost = cpm.simulate(float(x0), cert + ["E"], True, T=T, lock=lock, log=log)
    b = next((t for t in range(1, len(log))
              if log[t][2] > log[t - 1][2] - cpm.ALPHA * log[t - 1][3] + cpm.BREAK_TOL and abs(log[t][0]) > 1e-3), None)
    n_breaks = sum(1 for t in range(1, len(log))
                   if log[t][2] > log[t - 1][2] - cpm.ALPHA * log[t - 1][3] + cpm.BREAK_TOL and abs(log[t][0]) > 1e-3)
    out = dict(cost=cost, breaks=n_breaks, v0=log[0][2] if log else np.inf, viol=0.0, thm2=None)
    if b is not None:
        out["viol"] = log[b][2] - (log[b - 1][2] - cpm.ALPHA * log[b - 1][3])
        v0 = log[0][2]
        out["thm2"] = (v0 - log[b - 1][2]) / cpm.ALPHA + log[b - 1][3] + certified_value(log[b][0], cert)
    return out


def main():
    t0 = time.time()
    fams = list(families())
    print(f"{len(fams)} error families x {len(X0S)} initial states, T={T}\n")

    # naive minimum (no test) vs certified, certified units {A, B}, alpha = 0.5
    cpm.ALPHA = 0.5
    worst_naive, worst_cert = {}, {}
    for kind, name, J in fams:
        cpm.UNITS["E"] = (2, J, False)
        rn = [ratio(cpm.simulate(float(x0), ["A", "B", "E"], False, T=T), x0) for x0 in X0S]
        rc = [ratio(cpm.simulate(float(x0), ["A", "B", "E"], True, T=T), x0) for x0 in X0S]
        worst_naive.setdefault(kind, []).append(max(r for r in rn if np.isfinite(r)))
        worst_cert.setdefault(kind, []).append(max(r for r in rc if np.isfinite(r)))
    print("Naive minimum vs certified, units {A,B,E}, alpha=0.5: worst J/J* over x0, then max over the family")
    print(f"{'family':10s} {'members':>7s} {'naive: max worst':>17s} {'members with worst>1.5':>23s} {'cert: max worst':>16s}")
    for kind in worst_naive:
        wn, wc = worst_naive[kind], worst_cert[kind]
        print(f"{kind:10s} {len(wn):7d} {max(wn):17.3f} {sum(w > 1.5 for w in wn):23d} {max(wc):16.3f}")
    print()

    print("Genuine certificate breaks")
    print(f"{'alpha':>5s} {'cert':>4s} {'families w/ breaks':>18s} {'runs w/ breaks':>14s} {'max violation':>13s} "
          f"{'max |cost lock - no lock|':>25s} {'cost > v0/alpha':>15s} {'cost > Thm 2':>12s} {'cert worst J/J*':>15s}")
    for alpha, cert in itertools.product([0.5, 0.9], [["A"], ["B"], ["A", "B"]]):
        cpm.ALPHA = alpha
        fam_hits = runs = exceed1 = exceed2 = 0
        max_viol = max_diff = 0.0
        worst = 1.0
        hit_kinds = set()
        for kind, name, J in fams:
            cpm.UNITS["E"] = (2, J, False)
            hit = False
            for x0 in X0S:
                r = run_cert(cert, x0, lock=True)
                if np.isfinite(r["cost"]):
                    worst = max(worst, ratio(r["cost"], x0))
                if r["breaks"] == 0:
                    continue
                hit = True
                runs += 1
                max_viol = max(max_viol, r["viol"])
                free = run_cert(cert, x0, lock=False)
                if np.isfinite(r["cost"]) and np.isfinite(free["cost"]):
                    max_diff = max(max_diff, abs(r["cost"] - free["cost"]))
                exceed1 += r["cost"] > r["v0"] / alpha + 1e-9
                exceed2 += r["cost"] > r["thm2"] + 1e-6
            if hit:
                fam_hits += 1
                hit_kinds.add(kind)
        kinds = ",".join(sorted(hit_kinds)) if hit_kinds else "-"
        print(f"{alpha:5.1f} {''.join(cert):>4s} {fam_hits:7d} ({kinds:>9s}) {runs:14d} {max_viol:13.4f} "
              f"{max_diff:25.2e} {exceed1:15d} {exceed2:12d} {worst:15.3f}", flush=True)
    print(f"\ndone in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
