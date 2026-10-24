"""Sanity check: sample-set pruning in (linear, nominal) learning MPC.

Mass-spring-damper from RALMPC main.m with true parameters, no disturbance.
Compares three sample-set policies after every iteration:
  full    : keep every stored state (what the LMPC theory assumes)
  xhull   : keep only vertices of the convex hull of the x-projection
            (+ the initial trajectory), mimicking get_updateSampleSet in solve_RALMPC.m
  lifted  : drop a point only if it is dominated in the lifted space (x, J), i.e.
            some convex combination of other points has the same x and no larger J
"""
import numpy as np
import cvxpy as cp
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

Ts, m, k_star, c_star = 0.1, 1.0, 0.5, 0.3
A = np.array([[1.0, Ts], [-Ts * k_star / m, 1 - Ts * c_star / m]])
B = np.array([[0.0], [Ts / m]])
Q = np.diag([1.0, 1e-2])
R = np.array([[1e-1]])
x_lo, x_hi = np.array([-0.2, -5.0]), np.array([4.1, 5.0])
u_max = 15.0
x0 = np.array([4.0, 0.0])
N = 3
ITERS = 8
TOL = 1e-4


def stage(x, u):
    return float(x @ Q @ x + u @ R @ u)


def initial_trajectory():
    # sluggish saturated state feedback that respects x1 >= -0.2
    K0 = np.array([[0.5, 2.0]])
    xs, us = [x0.copy()], []
    x = x0.copy()
    for _ in range(600):
        u = np.clip(-K0 @ x, -u_max, u_max)
        x = A @ x + B @ u
        us.append(u)
        xs.append(x.copy())
        if np.linalg.norm(x) < 1e-6:
            break
    xs, us = np.array(xs), np.array(us)
    assert np.all(xs >= x_lo - 1e-9) and np.all(xs <= x_hi + 1e-9), "initial trajectory infeasible"
    return xs, us


def cost_to_go(xs, us):
    J = np.zeros(len(xs))
    for t in range(len(us) - 1, -1, -1):
        J[t] = stage(xs[t], us[t]) + J[t + 1]
    J[-1] = 0.0
    # last state is ~origin; force exact origin with zero cost so the hull contains 0
    return J


def lmpc_step(x, Z, Jz):
    n_pts = Z.shape[0]
    X = cp.Variable((N + 1, 2))
    U = cp.Variable((N, 1))
    lam = cp.Variable(n_pts, nonneg=True)
    cons = [X[0] == x, cp.sum(lam) == 1, X[N] == Z.T @ lam]
    obj = 0
    for k in range(N):
        cons += [X[k + 1] == A @ X[k] + B @ U[k],
                 X[k] >= x_lo, X[k] <= x_hi, cp.abs(U[k]) <= u_max]
        obj += cp.quad_form(X[k], Q) + cp.quad_form(U[k], R)
    obj += Jz @ lam
    prob = cp.Problem(cp.Minimize(obj), cons)
    prob.solve(solver=cp.CLARABEL)
    if prob.status not in ("optimal", "optimal_inaccurate"):
        return None, None
    return U.value[0], prob.value


def prune_lifted(Z, Jz):
    keep = np.ones(len(Z), dtype=bool)
    order = np.argsort(-Jz)  # try to remove expensive points first
    for i in order:
        if np.linalg.norm(Z[i]) < 1e-9:
            continue  # always keep the origin
        idx = np.where(keep)[0]
        idx = idx[idx != i]
        # min sum mu_j J_j  s.t.  sum mu_j z_j = z_i, sum mu_j = 1, mu >= 0
        A_eq = np.vstack([Z[idx].T, np.ones(len(idx))])
        b_eq = np.concatenate([Z[i], [1.0]])
        res = linprog(Jz[idx], A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
        if res.status == 0 and res.fun <= Jz[i] + 1e-9:
            keep[i] = False
    return Z[keep], Jz[keep]


def prune_xhull(Z, Jz, Z0, J0):
    hull = ConvexHull(Z)
    Zk, Jk = Z[hull.vertices], Jz[hull.vertices]
    allZ = np.vstack([Zk, Z0])
    allJ = np.concatenate([Jk, J0])
    _, uniq = np.unique(np.round(allZ, 12), axis=0, return_index=True)
    return allZ[np.sort(uniq)], allJ[np.sort(uniq)]


def run(policy):
    xs0, us0 = initial_trajectory()
    J0 = cost_to_go(xs0, us0)
    Z0, Jz0 = xs0.copy(), J0.copy()
    Z, Jz = Z0.copy(), Jz0.copy()
    log = []
    for it in range(ITERS):
        x = x0.copy()
        xs, us, Vs = [x.copy()], [], []
        decrease_violations = 0
        while np.linalg.norm(x) > TOL and len(us) < 200:
            u, V = lmpc_step(x, Z, Jz)
            if u is None:
                log.append(dict(it=it, status="INFEASIBLE", t=len(us)))
                return log
            Vs.append(V)
            us.append(u)
            x = A @ x + B @ u
            xs.append(x.copy())
        xs, us = np.array(xs), np.array(us)
        for t in range(len(Vs) - 1):
            # LMPC Lyapunov property: V(x_{t+1}) <= V(x_t) - l(x_t, u_t)
            if Vs[t + 1] > Vs[t] - stage(xs[t], us[t]) + 1e-6 * max(1.0, Vs[t]):
                decrease_violations += 1
        Jcl = sum(stage(xs[t], us[t]) for t in range(len(us)))
        Jtraj = cost_to_go(xs, us)
        Z = np.vstack([Z, xs])
        Jz = np.concatenate([Jz, Jtraj])
        if policy == "lifted":
            Z, Jz = prune_lifted(Z, Jz)
        elif policy == "xhull":
            Z, Jz = prune_xhull(Z, Jz, Z0, Jz0)
        log.append(dict(it=it, status="ok", J=Jcl, V0=Vs[0], steps=len(us),
                        pts=len(Z), viol=decrease_violations))
    return log


if __name__ == "__main__":
    results = {p: run(p) for p in ("full", "lifted", "xhull")}
    print(f"N={N}, iterations={ITERS}")
    print(f"{'it':>3} | {'full: J (pts, viol)':>26} | {'lifted: J (pts, viol)':>26} | {'xhull: J (pts, viol)':>26}")
    for it in range(ITERS):
        row = []
        for p in ("full", "lifted", "xhull"):
            r = results[p][it] if it < len(results[p]) else None
            if r is None:
                row.append("-")
            elif r["status"] != "ok":
                row.append(f"INFEASIBLE at t={r['t']}")
            else:
                row.append(f"{r['J']:.6f} ({r['pts']:4d}, {r['viol']})")
        print(f"{it:>3} | {row[0]:>26} | {row[1]:>26} | {row[2]:>26}")
