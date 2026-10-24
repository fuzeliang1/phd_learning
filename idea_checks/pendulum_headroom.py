"""Headroom of a short-horizon certified NMPC on pendulum swing-up (stage 1 of direction one).

System: Gym Pendulum-v1 dynamics (m = l = 1, g = 10, dt = 0.05, semi-implicit Euler)
    w+ = w + (15 sin th + 3 u) dt,   th+ = th + w+ dt,   |u| <= u_max,   th = 0 upright.
Two changes make the NMPC smooth: the speed limit |w| <= 8 is a constraint instead of Gym's clipping,
and the angle cost is 2 (1 - cos th) (equal to th^2 near upright) instead of Gym's wrapped th^2.
    g(x, u) = 2 (1 - cos th) + 0.1 w^2 + 0.001 u^2.

What is computed:
  J_ref   grid value iteration (periodic th, bilinear interpolation, speed limit as a penalty), checked
          by simulating its greedy policy; the reference is the greedy closed-loop cost, an achievable
          cost and so an upper bound on J*. Gaps measured against it underestimate the true gaps.
  units   LQR around upright for three design weights; terminal cost (1 + eps) * cost of u = K x on the
          linearised model, terminal set the largest ellipsoid on which J_f(f(x, Kx)) + g <= J_f(x),
          |Kx| <= u_max and |w+| <= 8 hold on a fine sample (numerical certificate, not a proof).
  NMPC    certified parallel NMPC with N steps (CasADi + IPOPT, multiple shooting): one problem per
          unit and per upright copy th = 2 pi k, initial guess from the greedy policy on J_ref.
          coverage = share of x0 where some problem is feasible (IPOPT is local, so this is a lower
          bound); closed loop runs NMPC until x enters a terminal set, then the unit's LQR law.
"""
import argparse
import multiprocessing as mp
import time

import casadi as ca
import numpy as np
from scipy.linalg import solve_discrete_are, solve_discrete_lyapunov

DT, W_MAX, QW, RU = 0.05, 8.0, 0.1, 0.001
PEN = 1e3            # value-iteration penalty for violating |w+| <= 8
STEPS = 400          # closed-loop length (20 s)


def wrap(th):
    return (th + np.pi) % (2 * np.pi) - np.pi


def step(th, w, u):
    w2 = w + (15.0 * np.sin(th) + 3.0 * u) * DT
    return th + w2 * DT, w2


def cost(th, w, u):
    return 2 * (1 - np.cos(th)) + QW * w ** 2 + RU * u ** 2


# ----------------------------------------------------------------------------- grid value iteration
class Grid:
    def __init__(self, n_th, n_w, n_u, u_max):
        self.th = np.linspace(-np.pi, np.pi, n_th, endpoint=False)   # n_th even: th = 0 is a node
        self.w = np.linspace(-W_MAX, W_MAX, n_w)                      # n_w odd: w = 0 is a node
        self.u = np.linspace(-u_max, u_max, n_u)                      # n_u odd: u = 0 is a level
        self.n_th, self.n_w = n_th, n_w
        TH, W, U = np.meshgrid(self.th, self.w, self.u, indexing="ij")
        th2, w2 = step(TH, W, U)
        self.stage = (cost(TH, W, U) + PEN * (np.abs(w2) > W_MAX)).reshape(-1, n_u)
        self.idx, self.wts = self.corners(th2, w2)
        self.idx = self.idx.reshape(-1, n_u, 4)
        self.wts = self.wts.reshape(-1, n_u, 4)

    def corners(self, th, w):
        sh = np.shape(th)
        a = (wrap(th) + np.pi) / (2 * np.pi) * self.n_th
        i0 = np.floor(a).astype(int) % self.n_th
        ta = a - np.floor(a)
        b = np.clip((w + W_MAX) / (2 * W_MAX) * (self.n_w - 1), 0, self.n_w - 1 - 1e-9)
        j0 = np.floor(b).astype(int)
        tb = b - j0
        i1 = (i0 + 1) % self.n_th
        idx = np.stack([i0 * self.n_w + j0, i1 * self.n_w + j0, i0 * self.n_w + j0 + 1, i1 * self.n_w + j0 + 1], -1)
        wts = np.stack([(1 - ta) * (1 - tb), ta * (1 - tb), (1 - ta) * tb, ta * tb], -1)
        return idx.reshape(*sh, 4), wts.reshape(*sh, 4)

    def value_iteration(self, tol=1e-6, max_iter=20000):
        J = np.zeros(self.n_th * self.n_w)
        for k in range(max_iter):
            Jn = (self.stage + (self.wts * J[self.idx]).sum(-1)).min(axis=1)
            if np.max(np.abs(Jn - J)) < tol:
                return Jn, k + 1
            J = Jn
        raise RuntimeError("value iteration did not converge")

    def greedy(self, J, th, w):
        """Greedy action on J at continuous states (arrays)."""
        TH, U = th[:, None], self.u[None, :]
        th2, w2 = step(TH, w[:, None], U)
        idx, wts = self.corners(th2, w2)
        q = cost(TH, w[:, None], U) + PEN * (np.abs(w2) > W_MAX) + (wts * J[idx]).sum(-1)
        return self.u[np.argmin(q, axis=1)]

    def rollout(self, J, th0, w0, steps):
        th, w = np.array(th0, float), np.array(w0, float)
        total, viol = np.zeros_like(th), np.zeros_like(th, dtype=bool)
        traj = []
        for _ in range(steps):
            u = self.greedy(J, th, w)
            traj.append((th.copy(), w.copy(), u.copy()))
            total += cost(th, w, u)
            th, w = step(th, w, u)
            viol |= np.abs(w) > W_MAX + 1e-9
        done = (np.abs(wrap(th)) < 0.05) & (np.abs(w) < 0.05)
        return total, done, viol, traj


# ----------------------------------------------------------------------------- certified LQR units
def lqr_units(u_max, r_designs=(0.001, 0.1, 1.0), eps=0.1):
    A = np.array([[1 + 15 * DT ** 2, DT], [15 * DT, 1.0]])
    B = np.array([[3 * DT ** 2], [3 * DT]])
    Q = np.diag([1.0, QW])
    s = np.linspace(-1.5, 1.5, 801)
    TH, W = [a.ravel() for a in np.meshgrid(s, 4 * s, indexing="ij")]
    X = np.stack([TH, W], 1)
    units = []
    for rd in r_designs:
        Pd = solve_discrete_are(A, B, Q, np.array([[rd]]))
        K = -np.linalg.solve(rd + B.T @ Pd @ B, B.T @ Pd @ A).ravel()
        P = (1 + eps) * solve_discrete_lyapunov((A + B @ K[None, :]).T, Q + RU * np.outer(K, K))
        V = np.einsum("ni,ij,nj->n", X, P, X)
        u = X @ K
        th2, w2 = step(TH, W, u)
        X2 = np.stack([th2, w2], 1)
        V2 = np.einsum("ni,ij,nj->n", X2, P, X2)
        edge = (np.abs(TH) >= s[-1]) | (np.abs(W) >= 4 * s[-1])
        ok = (np.abs(u) <= u_max) & (np.abs(w2) <= W_MAX) & (V2 + cost(TH, W, u) <= V) & ~edge
        units.append({"r": rd, "K": K, "P": P, "c": float(np.min(V[~ok]) * (1 - 1e-6))})
    return units


def in_terminal(th, w, unit):
    x = np.array([wrap(th), w])
    return x @ unit["P"] @ x <= unit["c"]


# ----------------------------------------------------------------------------- certified NMPC
class NMPC:
    """N-step problem for one unit; parameters x0 (2) and the upright copy th_t = 2 pi k."""

    def __init__(self, N, unit, u_max):
        X, U = ca.SX.sym("X", 2, N + 1), ca.SX.sym("U", N)
        p = ca.SX.sym("p", 3)
        obj, g, lbg, ubg = 0, [X[:, 0] - p[:2]], [0, 0], [0, 0]
        for k in range(N):
            th, w = X[0, k], X[1, k]
            w2 = w + (15 * ca.sin(th) + 3 * U[k]) * DT
            g += [X[:, k + 1] - ca.vertcat(th + w2 * DT, w2)]
            lbg += [0, 0]
            ubg += [0, 0]
            obj += 2 * (1 - ca.cos(th)) + QW * w ** 2 + RU * U[k] ** 2
        e = X[:, N] - ca.vertcat(p[2], 0)
        P = ca.DM(unit["P"])
        obj += ca.mtimes([e.T, P, e])
        g += [ca.mtimes([e.T, P, e])]
        lbg += [-ca.inf]
        ubg += [unit["c"]]
        z = ca.vertcat(ca.reshape(X, -1, 1), U)
        opts = {"ipopt.print_level": 0, "print_time": 0, "ipopt.max_iter": 300, "ipopt.tol": 1e-8}
        self.solver = ca.nlpsol("s", "ipopt", {"x": z, "f": obj, "g": ca.vertcat(*g), "p": p}, opts)
        nX = 2 * (N + 1)
        lbz = -np.inf * np.ones(nX + N)
        ubz = np.inf * np.ones(nX + N)
        lbz[1:nX:2], ubz[1:nX:2] = -W_MAX, W_MAX       # |w_k| <= 8
        lbz[nX:], ubz[nX:] = -u_max, u_max
        self.lbz, self.ubz, self.lbg, self.ubg, self.N, self.nX = lbz, ubz, lbg, ubg, N, nX

    def solve(self, x0, th_t, guess):
        t = time.perf_counter()
        r = self.solver(x0=guess, p=np.r_[x0, th_t], lbx=self.lbz, ubx=self.ubz, lbg=self.lbg, ubg=self.ubg)
        _W.setdefault("times", []).append(time.perf_counter() - t)
        ok = self.solver.stats()["success"]
        z = np.array(r["x"]).ravel()
        if ok:   # double-check feasibility: IPOPT can report success at tolerance-level violations
            gv = np.array(r["g"]).ravel()
            ok = np.all(np.abs(gv[:-1]) < 1e-6) and gv[-1] <= self.ubg[-1] + 1e-6
        return (float(r["f"]) if ok else np.inf), z

    def guess_from(self, th, w, u, th_t):
        """Initial guess from a trajectory (arrays of length >= N+1 / N)."""
        z = np.empty(self.nX + self.N)
        z[0:self.nX:2], z[1:self.nX:2] = th[:self.N + 1], w[:self.N + 1]
        z[self.nX:] = u[:self.N]
        return z


_W = {}


def _init(N, units, u_max, grid_args, J, naive):
    _W["mpcs"] = [NMPC(N, un, u_max) for un in units]
    _W["units"], _W["u_max"], _W["naive"] = units, u_max, naive
    _W["grid"], _W["J"] = Grid(*grid_args), J


def guess_traj(th0, w0, n):
    g, J = _W["grid"], _W["J"]
    th, w, us = [th0], [w0], []
    for _ in range(n):
        u = float(g.greedy(J, np.array([th[-1]]), np.array([w[-1]]))[0])
        t2, w2 = step(th[-1], w[-1], u)
        th.append(t2)
        w.append(float(np.clip(w2, -W_MAX, W_MAX)))
        us.append(u)
    return np.array(th), np.array(w), np.array(us)


def run_x0(x0):
    """Coverage at x0 and closed-loop cost of the certified parallel NMPC."""
    mpcs, units = _W["mpcs"], _W["units"]
    N = mpcs[0].N
    th, w = float(x0[0]), float(x0[1])
    if _W["naive"]:      # no learned information: hold x0, zero torque, nearest upright copy
        gth, gw, gu = np.full(N + 1, th), np.full(N + 1, w), np.zeros(N)
        k0 = 0
    else:
        gth, gw, gu = guess_traj(th, w, N)
        k0 = int(np.round(gth[-1] / (2 * np.pi)))
    best = (np.inf, None, None, None)
    for i, m in enumerate(mpcs):
        for k in (k0 - 1, k0, k0 + 1):
            # shift the guess so that it ends near th = 2 pi k
            v, z = m.solve([th, w], 2 * np.pi * k, m.guess_from(gth + 2 * np.pi * (k - k0), gw, gu, 2 * np.pi * k))
            if v < best[0]:
                best = (v, i, k, z)
    if not np.isfinite(best[0]):
        return False, np.inf, 0, []
    _W["times"] = []                        # time only the warm-started closed-loop solves
    _, i, k, z = best
    total, t = 0.0, 0
    while t < STEPS:
        un = units[i]
        if in_terminal(th, w, un):          # inside the certified set: the unit's LQR law from here on
            for _ in range(STEPS - t):
                x = np.array([wrap(th), w])
                u = float(np.clip(x @ un["K"], -_W["u_max"], _W["u_max"]))
                total += cost(th, w, u)
                th, w = step(th, w, u)
            break
        m = mpcs[i]
        u = z[m.nX]
        total += cost(th, w, u)
        th, w = step(th, w, u)
        t += 1
        # warm start: shift the previous plan and append the LQR step
        Xp = z[:m.nX].reshape(-1, 2)
        Up = z[m.nX:]
        xe = Xp[-1] - [2 * np.pi * k, 0]
        ue = float(np.clip(xe @ units[i]["K"], -_W["u_max"], _W["u_max"]))
        te, we = step(Xp[-1, 0], Xp[-1, 1], ue)
        gX = np.vstack([Xp[1:], [te, we]])
        gU = np.r_[Up[1:], ue]
        cand = []
        for j, mj in enumerate(mpcs):
            vj, zj = mj.solve([th, w], 2 * np.pi * k, np.r_[gX.ravel(), gU])
            cand.append((vj, j, zj))
        v, i, z = min(cand, key=lambda c: c[0])
        if not np.isfinite(v):
            return True, np.nan, t, _W["times"]   # recursive feasibility lost: should not happen for certified units
    return True, total, t, _W["times"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--umax", type=float, default=2.0)
    ap.add_argument("--grid", type=int, nargs=3, default=[300, 241, 41], help="n_th n_w n_u for value iteration")
    ap.add_argument("--horizons", type=int, nargs="+", default=[10, 20, 40, 60, 80])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--naive-guess", action="store_true",
                    help="initialise IPOPT at x0 with zero torque instead of the greedy trajectory on J_ref")
    args = ap.parse_args()
    t0 = time.time()
    print(f"pendulum swing-up, u_max = {args.umax}, value-iteration grid {args.grid}, "
          f"IPOPT initial guess: {'x0 held, u = 0' if args.naive_guess else 'greedy trajectory on J_ref'}", flush=True)

    # reference: value iteration at two resolutions, checked by greedy rollouts
    th0 = np.linspace(-np.pi, np.pi, 12, endpoint=False)
    w0 = np.array([-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0])
    TH0, W0 = [a.ravel() for a in np.meshgrid(th0, w0, indexing="ij")]
    coarse = [args.grid[0] // 4 * 2, args.grid[1] // 2 | 1, args.grid[2] // 2 | 1]   # half resolution
    nz = ~((TH0 == 0) & (W0 == 0))          # the upright x0 has J* = 0 and is left out of all ratios
    refs = {}
    for name, gargs in (("coarse", coarse), ("fine", args.grid)):
        grid = Grid(*gargs, args.umax)
        J, it = grid.value_iteration()
        tot, done, viol, _ = grid.rollout(J, TH0, W0, STEPS)
        idx, wts = grid.corners(TH0, W0)
        Jx = (wts * J[idx]).sum(-1)
        refs[name] = (tot, J, grid)
        r = tot[nz] / Jx[nz]
        print(f"  value iteration {name} {gargs}: {it} iterations; greedy rollout reaches upright from "
              f"{done.mean():.0%} of x0, speed-limit violations {viol.mean():.0%}; "
              f"rollout cost / table value: mean {np.mean(r):.4f}, min {np.min(r):.4f}, "
              f"max {np.max(r):.4f}  ({time.time() - t0:.0f} s)", flush=True)
    tot_c, tot_f = refs["coarse"][0], refs["fine"][0]
    r = tot_c[nz] / tot_f[nz]
    print(f"  greedy rollout cost, coarse / fine grid: mean {np.mean(r):.4f}, min {np.min(r):.4f}, "
          f"max {np.max(r):.4f}", flush=True)
    jref = np.minimum(tot_c, tot_f)          # best achievable cost found, an upper bound on J*
    _, J, grid = refs["fine"]

    units = lqr_units(args.umax)
    for un in units:
        print(f"  LQR unit, design R = {un['r']}: K = {np.round(un['K'], 2)}, certified level c = {un['c']:.3f}",
              flush=True)

    print(f"\n  {'N':>3} {'horizon':>8} | {'coverage':>8} {'|w0|<=2':>8} | gap to J_ref on covered x0: "
          f"{'mean':>7} {'median':>7} {'max':>7}  share > 5%", flush=True)
    for N in args.horizons:
        with mp.Pool(args.workers, initializer=_init, initargs=(N, units, args.umax, args.grid + [args.umax], J, args.naive_guess)) as pool:
            res = pool.map(run_x0, list(zip(TH0, W0)))
        cov = np.array([r[0] for r in res])
        cl = np.array([r[1] for r in res])
        lost = np.sum(np.isnan(cl))
        ok = cov & np.isfinite(cl) & nz
        gap = cl[ok] / jref[ok] - 1
        slow = np.abs(W0) <= 2
        times = np.concatenate([r[3] for r in res if len(r[3])]) if any(len(r[3]) for r in res) else np.array([np.nan])
        stats = (f"{np.mean(gap):7.2%} {np.median(gap):7.2%} {np.max(gap):7.2%}  {np.mean(gap > 0.05):6.1%}"
                 if ok.any() else "   (none covered)")
        print(f"  {N:>3} {N * DT:>6.2f} s | {cov[nz].mean():8.1%} {cov[slow & nz].mean():8.1%} | {'':27}{stats}"
              + f"  | IPOPT per unit solve {1e3 * np.median(times):.1f} ms (median, warm start)"
              + (f"  [{lost} closed loops lost feasibility]" if lost else "") + f"  ({time.time() - t0:.0f} s)",
              flush=True)


if __name__ == "__main__":
    main()
