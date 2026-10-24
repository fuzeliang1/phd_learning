"""Python port of the numerical examples of Parallel MPC (Li, Karapetyan, Schmid, Lygeros, Johansson,
Martensson; IEEE TAC 71(2), 2026; arXiv 2309.14560v3, Section IV).

The MATLAB originals (MPT3 + YALMIP + LMI Lab) are in pr-mpc-main. File mapping:
  paper Example 4.1 (LQ with constraints, timing)  <- example_4_3.m, example_4_3_computational_complexity.m
  paper Example 4.2 (PWA dynamics, Table I)         <- example_4_6.m, functions/*.m
  paper Example 4.3 (switching system, Table II)    <- example_4_7.m, functions/vi_low_bound.m

Everything is recomputed with numpy/scipy/cvxpy: invariant sets by the backward iteration of
functions/max_reach.m and MPT's invariantSet, Lyapunov matrices by dlyap or by the min-trace LMI of
functions/Lyap_fun.m (solved to optimality here; LMI Lab's mincx stops at relative accuracy 1e-2),
lookahead problems as QPs/SOCPs with Clarabel. The one matrix that cannot be recomputed, the LMI gain
of Example 4.2 (LMI Lab's feasp returns an arbitrary feasible point), is taken from the paper's
Appendix B.

Known differences between the paper text and the code, and how they are handled, are listed in
README.md. Example 4.1 is run in both versions.
"""
import argparse

import cvxpy as cp
import numpy as np
from scipy.linalg import solve_discrete_are, solve_discrete_lyapunov
from scipy.optimize import linprog

XMAX, UMAX = 5.0, 1.0
BOX_H = np.vstack([np.eye(2), -np.eye(2)])
BOX_h = XMAX * np.ones(4)
Q2 = np.eye(2)


# ----------------------------------------------------------------------------- polytopes
def lp_max(c, H, h):
    res = linprog(-c, A_ub=H, b_ub=h, bounds=[(None, None)] * H.shape[1], method="highs")
    if res.status == 3:
        return np.inf
    if res.status != 0:
        raise RuntimeError(res.message)
    return -res.fun


def prune(H, h, tol=1e-9):
    """Normalise rows and drop redundant inequalities (functions/slim_constraint.m)."""
    nrm = np.linalg.norm(H, axis=1)
    keep = nrm > 1e-12
    if np.any(~keep & (h < 0)):
        raise ValueError("empty polytope")
    H, h = H[keep] / nrm[keep, None], h[keep] / nrm[keep]
    rows = list(range(len(h)))
    for i in reversed(range(len(h))):
        rest = [j for j in rows if j != i]
        if rest and lp_max(H[i], H[rest], h[rest]) <= h[i] + tol:
            rows.remove(i)
    return H[rows], h[rows]


def max_invariant(H, h, maps, tol=1e-7, max_iter=200):
    """Largest subset of {H x <= h} kept invariant by every linear map in `maps`.

    Same fixed point as functions/max_reach.m (both closed-loop maps applied everywhere) and, for a
    single map, as MPT's LTISystem.invariantSet().
    """
    H, h = prune(H, h)
    for _ in range(max_iter):
        pre = np.vstack([H @ M for M in maps])
        pre_h = np.tile(h, len(maps))
        if all(lp_max(r, H, h) <= b + tol for r, b in zip(pre, pre_h)):
            return H, h
        H, h = prune(np.vstack([H, pre]), np.concatenate([h, pre_h]))
    raise RuntimeError("invariant set iteration did not converge")


# ----------------------------------------------------------------------------- control helpers
def dlqr(A, B, Q, R):
    """Gain L with u = L x (MATLAB dlqr/idare return K with u = -K x)."""
    B = B.reshape(-1, 1)
    P = solve_discrete_are(A, B, Q, np.atleast_2d(R))
    K = np.linalg.solve(np.atleast_2d(R) + B.T @ P @ B, B.T @ P @ A)
    return -K.ravel(), P


def dlyap_cost(Acl, Qcl):
    """X with Acl' X Acl - X + Qcl = 0, i.e. MATLAB dlyap(Acl', Qcl)."""
    return solve_discrete_lyapunov(Acl.T, Qcl)


def common_lyap(Acls, Qcls):
    """min trace P s.t. [[P, P Ac], [Ac' P, P - Qc]] >= 0 for every mode (functions/Lyap_fun.m)."""
    P = cp.Variable((2, 2), symmetric=True)
    cons = [P >> 0]
    for Ac, Qc in zip(Acls, Qcls):
        M = cp.bmat([[P, P @ Ac], [Ac.T @ P, P - Qc]])
        cons.append((M + M.T) / 2 >> 0)
    cp.Problem(cp.Minimize(cp.trace(P)), cons).solve(solver=cp.CLARABEL)
    return (P.value + P.value.T) / 2


def ellipsoid_level(P, L):
    """Largest a with {x' P x <= a} inside |x|_inf <= 5, |L x| <= 1 (find_ellipsoid in example_4_3.m)."""
    Pinv = np.linalg.inv(P)
    rows = [(r, XMAX) for r in BOX_H] + [(L, UMAX), (-L, UMAX)]
    return min(b ** 2 / (r @ Pinv @ r) for r, b in rows)


# ----------------------------------------------------------------------------- lookahead problems
class Lookahead:
    """One lookahead minimisation with a fixed mode sequence.

    min  sum_{k<N} x_k' Q x_k + R u_k^2 + x_N' P x_N
    s.t. x_{k+1} = A_k x_k + B_k u_k,  |u_k| <= 1,  |x_k|_inf <= 5 (k = 0..N),
         terminal set (polytope or ellipsoid), optional region constraints a' x_k <= b.
    """

    def __init__(self, dyn, R, P, term=None, ellipse=None, regions=(), solver=None):
        N = len(dyn)
        self.x0 = cp.Parameter(2)
        x, u = cp.Variable((N + 1, 2)), cp.Variable(N)
        cons = [x[0] == self.x0, cp.abs(u) <= UMAX, cp.abs(x) <= XMAX]
        cons += [x[k + 1] == A @ x[k] + B * u[k] for k, (A, B) in enumerate(dyn)]
        if term is not None:
            cons.append(term[0] @ x[N] <= term[1])
        if ellipse is not None:
            cons.append(cp.quad_form(x[N], ellipse[0]) <= ellipse[1])
        cons += [a @ x[k] <= b for k, a, b in regions]
        obj = sum(cp.quad_form(x[k], Q2) + R * cp.square(u[k]) for k in range(N))
        if P is not None:
            obj += cp.quad_form(x[N], (P + P.T) / 2)
        self.prob = cp.Problem(cp.Minimize(obj), cons)
        self.x, self.u = x, u
        # Clarabel first; near the origin (|x| ~ 1e-5) it can stall, then fall back
        self.solvers = [solver] if solver else [cp.CLARABEL] + ([cp.GUROBI] if ellipse else [cp.PIQP, cp.GUROBI])
        self.solve_time = 0.0

    def __call__(self, x0):
        self.x0.value = np.asarray(x0, dtype=float)
        for s in self.solvers:
            try:
                self.prob.solve(solver=s, **({"Threads": 1} if s == cp.GUROBI else {}))
            except cp.error.SolverError:
                continue
            if self.prob.status in ("optimal", "infeasible", "infeasible_or_unbounded"):
                break
        self.solve_time = self.prob.solver_stats.solve_time or 0.0
        # the objective is bounded below by 0, so Gurobi's INF_OR_UNBD means infeasible
        if self.prob.status in ("infeasible", "infeasible_inaccurate", "infeasible_or_unbounded"):
            return np.inf, None
        if self.prob.status not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(f"lookahead solve failed: {self.prob.status}")
        return self.prob.value, self.u.value.copy()


def best(problems, x0):
    """Smallest value over a list of lookahead problems: (value, first control, index)."""
    out = [p(x0) for p in problems]
    i = int(np.argmin([v for v, _ in out]))
    v, u = out[i]
    return v, (None if u is None else u[0]), i


def fmt(v, nd=1):
    return "inf" if not np.isfinite(v) else f"{v:.{nd}f}"


# ============================================================================= Example 4.2 (PWA)
A1_PWA = np.array([[0.35, -0.6062], [0.6062, 0.35]])   # paper text rounds 0.6062 to 0.61
A2_PWA = np.array([[0.35, 0.6062], [-0.6062, 0.35]])
B_PWA = np.array([0.0, 1.0])
R_PWA = 0.4
X0_PWA = [(-4.0, 4.6), (5.0, 3.4), (-4.5, 2.7), (-5.0, -5.0)]
TABLE1 = {  # paper Table I: J~1, J~2, J~3, J~4, selected unit, T^l Jbar, J_mu~
    (-4.0, 4.6): (53.4, np.inf, 53.4, 54.3, 3, 53.4, 52.2),
    (5.0, 3.4): (53.1, np.inf, np.inf, 54.7, 1, 53.1, 52.0),
    (-4.5, 2.7): (39.2, 39.2, 39.1, 40.4, 3, 39.1, 38.4),
    (-5.0, -5.0): (np.inf, np.inf, 91.1, 85.9, 4, 85.9, 81.4),
}
# Table I numbers its columns like example_4_6.m (LQR first); Appendix B's text calls the Lazar unit
# "unit 1" and the LQR unit "unit 2". Units are kept by name here and printed in Table I's order.
UNITS_PWA = ["LQR", "Lazar", "LMI", "zero"]
# Appendix B values: the LMI gain is used as printed (it cannot be recomputed); K and S are an
# alternative set of terminal ingredients for the check below
APPB_L_LMI = (-np.array([0.5112, 0.2963]), np.array([0.5112, -0.2963]))
APPB_K = {"LQR": np.diag([1.4786, 1.9605]), "LMI": np.diag([1.7847, 2.1375]), "zero": np.diag([1.9631, 1.9631])}
APPB_S = {
    "LQR": np.array([[-0.4452, -0.3354], [0.4452, 0.3354], [0.4452, -0.3354], [-0.4452, 0.3354]]),
    "LMI": np.array([[0.15, -0.33], [-0.15, 0.33], [-0.15, -0.33], [0.15, 0.33],
                     [-0.51, -0.30], [0.51, 0.30], [0.51, -0.30], [-0.51, 0.30]]),
}


def pwa_mode(x):
    return 0 if x[0] >= 0 else 1


def pwa_step(x, u):
    return (A1_PWA if x[0] >= 0 else A2_PWA) @ x + B_PWA * u


def pwa_units():
    """Terminal ingredients of the four units: mode-dependent gains u = L^m x, min-trace common
    Lyapunov matrix (functions/Lyap_fun.m) and invariant set (functions/max_reach.m)."""
    gains = {
        "LQR": (dlqr(A1_PWA, B_PWA, Q2, R_PWA)[0], dlqr(A2_PWA, B_PWA, Q2, R_PWA)[0]),
        "Lazar": (-np.array([0.611, 0.3572]), np.array([0.611, -0.3572])),
        "LMI": APPB_L_LMI,
        "zero": (np.zeros(2), np.zeros(2)),
    }
    units = {}
    for name, (L1, L2) in gains.items():
        Ac1, Ac2 = A1_PWA + np.outer(B_PWA, L1), A2_PWA + np.outer(B_PWA, L2)
        P = common_lyap([Ac1, Ac2], [Q2 + R_PWA * np.outer(L1, L1), Q2 + R_PWA * np.outer(L2, L2)])
        # functions/safe_set.m imposes both modes' input bounds everywhere, max_reach.m both maps
        H0 = np.vstack([L1, -L1, L2, -L2, BOX_H])
        h0 = np.concatenate([UMAX * np.ones(4), BOX_h])
        S = max_invariant(H0, h0, [Ac1, Ac2])
        units[name] = {"L": (L1, L2), "P": P, "S": S}
    return units


def pwa_lookaheads(P, S, ell=2):
    """All mode sequences for one unit; the mode of x_0 is fixed by x_0, later ones are chosen."""
    probs = []
    for seq in np.ndindex(*(2,) * ell):
        dyn = [((A1_PWA, A2_PWA)[m], B_PWA) for m in seq]
        regions = [(k, np.array([-1.0, 0.0]) if m == 0 else np.array([1.0, 0.0]), 0.0)
                   for k, m in enumerate(seq) if k > 0]
        probs.append((seq[0], Lookahead(dyn, R_PWA, P, term=S, regions=regions)))
    return probs


def pwa_unit_value(probs, x):
    cand = [p for m0, p in probs if m0 == pwa_mode(x) or x[0] == 0]
    return best(cand, x)[:2]


def pwa_closed_loop(policy, x0, T=50):
    """Cost as in closed_loop_cost_calculator of example_4_6.m: T-1 steps plus x_{T-1}' Q x_{T-1}."""
    x, J = np.array(x0, dtype=float), 0.0
    for _ in range(T - 1):
        v, u = policy(x)
        if not np.isfinite(v):
            return np.inf
        J += x @ x + R_PWA * u ** 2
        x = pwa_step(x, u)
    return J + x @ x


def pwa_lower_bound(x0, N):
    """(T^N 0)(x0) <= J*(x0) for the PWA system, as an MIQP in mixed logical dynamical form.

    Written in gurobipy directly: through cvxpy the quadratic objective becomes second-order cones,
    which exceeds the size limit of the pip-installed Gurobi licence for N > 10.
    """
    import gurobipy as gp
    M = 20.0
    m = gp.Model()
    m.Params.OutputFlag, m.Params.Threads, m.Params.MIPGap = 0, 1, 1e-9
    x = m.addVars(N + 1, 2, lb=-XMAX, ub=XMAX)
    u = m.addVars(N, lb=-UMAX, ub=UMAX)
    d = m.addVars(N, vtype=gp.GRB.BINARY)          # d_k = 1: x_k[0] >= 0, mode A1
    m.addConstrs(x[0, i] == x0[i] for i in range(2))
    for k in range(N):
        m.addConstr(x[k, 0] >= -M * (1 - d[k]))
        m.addConstr(x[k, 0] <= M * d[k])
        for A, on in ((A1_PWA, 1 - d[k]), (A2_PWA, d[k])):
            for i in range(2):
                nxt = A[i, 0] * x[k, 0] + A[i, 1] * x[k, 1] + B_PWA[i] * u[k]
                m.addConstr(x[k + 1, i] - nxt <= M * on)
                m.addConstr(nxt - x[k + 1, i] <= M * on)
    m.setObjective(gp.quicksum(x[k, 0] * x[k, 0] + x[k, 1] * x[k, 1] + R_PWA * u[k] * u[k] for k in range(N)))
    m.optimize()
    if m.Status != gp.GRB.OPTIMAL:
        raise RuntimeError(f"lower bound MIQP status {m.Status}")
    return m.ObjVal


def same(a, b):
    return (not np.isfinite(a) and not np.isfinite(b)) or (np.isfinite(a) and np.isfinite(b) and round(a, 1) == b)


def example_4_2():
    print("=" * 100)
    print("Example 4.2 (PWA dynamics), paper Table I")
    units = pwa_units()
    for name, d in units.items():
        H, h = d["S"]
        print(f"  {name:>5}: L^1 = {np.round(d['L'][0], 4)}, L^2 = {np.round(d['L'][1], 4)}, "
              f"P = diag{tuple(np.round(np.diag(d['P']), 4).tolist())} (off-diagonal {abs(d['P'][0, 1]):.0e}), "
              f"S has {len(h)} facets")
        if name in APPB_K:
            print(f"         Appendix B: K = diag{tuple(np.diag(APPB_K[name]).tolist())}")
        if name in APPB_S:
            print(f"         S written as P x <= 1: {np.round(H / h[:, None], 4).tolist()}")

    variants = {"recomputed": {n: (d["P"], d["S"]) for n, d in units.items()}}
    appb = {}
    for n in UNITS_PWA:
        S = (APPB_S[n], np.ones(len(APPB_S[n]))) if n in APPB_S else units[n]["S"]
        appb[n] = (APPB_K.get(n, units[n]["P"]), S)
    variants["Appendix B K and S (Lazar unit recomputed, not printed)"] = appb

    lbs = {}
    for name, ingr in variants.items():
        probs = {n: pwa_lookaheads(P, S) for n, (P, S) in ingr.items()}

        def parallel(x):
            vals = [pwa_unit_value(probs[n], x) for n in UNITS_PWA]
            j = int(np.argmin([v for v, _ in vals]))
            return vals[j]

        print(f"\n  terminal ingredients: {name}")
        print("  columns in Table I order: 1 = LQR, 2 = Lazar, 3 = LMI, 4 = zero; the last column is the")
        print("  lower bound (T^15 0)(x0) <= J*(x0); gaps are relative to it")
        print(f"  {'x0':>13} | {'J~1':>6} {'J~2':>6} {'J~3':>6} {'J~4':>6} | {'i~':>2} | {'T^l Jbar':>8} | "
              f"{'J_mu~':>8} | single-unit closed loops J_mu~i | Table I | {'LB':>7} gap(J_mu~) gap(best single)")
        for x0 in X0_PWA:
            Jt = [pwa_unit_value(probs[n], np.array(x0))[0] for n in UNITS_PWA]
            isel = int(np.argmin(Jt)) + 1
            Jmu = pwa_closed_loop(parallel, x0)
            singles = [pwa_closed_loop(lambda x, n=n: pwa_unit_value(probs[n], x), x0) for n in UNITS_PWA]
            ours = Jt + [isel, min(Jt), Jmu]
            hits = sum(same(a, b) if k != 4 else a == b for k, (a, b) in enumerate(zip(ours, TABLE1[x0])))
            lb = lbs.setdefault(x0, pwa_lower_bound(np.array(x0), 15))
            print(f"  {str(x0):>13} | " + " ".join(f"{fmt(v):>6}" for v in Jt) +
                  f" | {isel:>2} | {fmt(min(Jt)):>8} | {fmt(Jmu, 3):>8} | " +
                  " ".join(f"{fmt(v, 3):>7}" for v in singles) + f" | {hits}/7     | {lb:7.3f} "
                  f"{(Jmu - lb) / Jmu:9.1e}  {(min(singles) - lb) / min(singles):9.1e}")


# ============================================================================= Example 4.3 (switching)
A_SW = (np.array([[2.0, 1.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [0.0, 0.5]]))
B_SW = (np.array([1.0, 1.0]), np.array([1.0, 2.0]))
X0_SW = [(-4.0, 4.6), (1.2, 1.5), (-3.5, 2.0), (-1.5, -0.5)]
TABLE2 = {  # paper Table II: T^l Jbar, J_mu~, reported bound on (J_mu~ - T^8 Jbar_0) / J_mu~
    (-4.0, 4.6): (65.8, 65.8, 1e-5), (1.2, 1.5): (89.2, 86.6, 1e-3),
    (-3.5, 2.0): (123.3, 113.5, 1e-4), (-1.5, -0.5): (34.6, 34.6, 1e-4),
}


def example_4_3(T=80, n_lb=8):
    print("=" * 100)
    print("Example 4.3 (switching system), simplified parallel MPC, paper Table II")
    term = []
    for d in range(2):
        L, P = dlqr(A_SW[d], B_SW[d], Q2, 1.0)
        S = max_invariant(np.vstack([BOX_H, L, -L]), np.concatenate([BOX_h, [UMAX, UMAX]]),
                          [A_SW[d] + np.outer(B_SW[d], L)])
        term.append((P, S))
        print(f"  mode {d + 1}: L = {np.round(L, 4)}, K = {np.round(P, 4).tolist()}, S has {len(S[1])} facets")
    # example_4_7.m: unit (first mode, later modes, terminal ingredients of the later mode)
    seqs = [(0, 0), (1, 0), (1, 1), (0, 1)]
    units = []
    for d0, d in seqs:
        modes = [d0] + [d] * 4
        units.append((d0, Lookahead([(A_SW[m], B_SW[m]) for m in modes], 1.0, term[d][0], term=term[d][1])))
    lb_probs = [Lookahead([(A_SW[m], B_SW[m]) for m in modes], 1.0, None) for modes in np.ndindex(*(2,) * n_lb)]

    print(f"  {'x0':>12} | {'T^l Jbar':>8} {'J_mu~':>8} {'LB=T^8 0':>9} {'gap':>8} | paper: T^l Jbar, J_mu~, gap")
    for x0 in X0_SW:
        x, J = np.array(x0), 0.0
        for j in range(T + 1):
            vals = [p(x) for _, p in units]
            i = int(np.argmin([v for v, _ in vals]))
            if j == 0:
                v0 = vals[i][0]
            v = vals[i][1][0]
            J += x @ x + v ** 2
            x = A_SW[units[i][0]] @ x + B_SW[units[i][0]] * v
        J += x @ term[0][0] @ x        # example_4_7.m adds x_end' P_Ly{1} x_end
        lb = min(p(np.array(x0))[0] for p in lb_probs)
        gap = (J - lb) / J
        ref = TABLE2[x0]
        print(f"  {str(x0):>12} | {v0:8.1f} {J:8.1f} {lb:9.2f} {gap:8.1e} | {ref[0]}, {ref[1]}, <{ref[2]:.0e}")


# ============================================================================= Example 4.1 (LQ)
A_LQ = np.array([[1.0, 1.0], [0.0, 1.0]])
B_LQ = np.array([1.0, 0.5])
X0_LQ = [(-5.0, 2.7), (2.3, -0.6)]
CONFIGS_LQ = {   # lookahead lengths and the unit whose terminal set is an ellipsoid
    "paper text (l = 2,2,3,3; ellipsoid on unit 4)": ((2, 2, 3, 3), 3),
    "code (l = 3,3,3,3; ellipsoid on unit 3)": ((3, 3, 3, 3), 2),
}


def lq_units():
    gains = [dlqr(A_LQ, B_LQ, Q2, 1.0)[0], -np.array([0.1, 1.2]), -np.array([0.2, 0.7]), -np.array([0.3, 0.8])]
    units = []
    for L in gains:
        Acl = A_LQ + np.outer(B_LQ, L)
        P = dlyap_cost(Acl, Q2 + np.outer(L, L))
        S = max_invariant(np.vstack([BOX_H, L, -L]), np.concatenate([BOX_h, [UMAX, UMAX]]), [Acl])
        units.append({"L": L, "P": P, "S": S, "E": (P, ellipsoid_level(P, L))})
    return units


def lq_lookaheads(units, ells, e_idx, solver=None):
    out = []
    for i, (d, ell) in enumerate(zip(units, ells)):
        dyn = [(A_LQ, B_LQ)] * ell
        if i == e_idx:
            out.append(Lookahead(dyn, 1.0, d["P"], ellipse=d["E"], solver=solver))
        else:
            out.append(Lookahead(dyn, 1.0, d["P"], term=d["S"], solver=solver))
    return out


class Alternative:
    """Left side of eq. (14): one lookahead of length max(l_i) with binary unit selection (big-M).

    This is what example_4_3_computational_complexity.m intends; there the loop
    `for i = 1:g, p(g) = binvar(1); end` creates only p(4), so p(1..3) are the constant 0 and the
    'naive' controller always uses unit 4.
    """

    def __init__(self, units, ells, e_idx):
        L0, L1 = min(ells), max(ells)
        self.x0 = cp.Parameter(2)
        x, u = cp.Variable((L1 + 1, 2)), cp.Variable(L1)
        z, w = cp.Variable(len(units), boolean=True), cp.Variable(len(units), nonneg=True)
        cons = [x[0] == self.x0, cp.abs(u) <= UMAX, cp.sum(z) == 1]
        cons += [x[k + 1] == A_LQ @ x[k] + B_LQ * u[k] for k in range(L1)]
        cons += [cp.abs(x[k]) <= XMAX for k in range(L0 + 1)]
        short = sum(z[i] for i, l in enumerate(ells) if l < L1) if L0 < L1 else 0
        cons += [cp.abs(x[k]) <= XMAX + 100 * short for k in range(L0 + 1, L1 + 1)]
        for i, (d, ell) in enumerate(zip(units, ells)):
            tail = sum(cp.sum_squares(x[k]) + cp.square(u[k]) for k in range(L0, ell)) + cp.quad_form(x[ell], d["P"])
            big = 51 * (ell - L0) + 50 * np.linalg.eigvalsh(d["P"]).max() + 1
            cons.append(tail <= w[i] + big * (1 - z[i]))
            if i == e_idx:
                Pe, a = d["E"]
                cons.append(cp.quad_form(x[ell], Pe) <= a + (50 * np.linalg.eigvalsh(Pe).max() - a) * (1 - z[i]))
            else:
                H, h = d["S"]
                bigH = np.maximum(XMAX * np.abs(H).sum(axis=1) - h, 0) + 100 * (ell > L0)
                cons.append(H @ x[ell] <= h + cp.multiply(bigH, 1 - z[i]))
        obj = sum(cp.sum_squares(x[k]) + cp.square(u[k]) for k in range(L0)) + cp.sum(w)
        self.prob, self.u, self.z = cp.Problem(cp.Minimize(obj), cons), u, z

    def __call__(self, x0):
        self.x0.value = np.asarray(x0, dtype=float)
        # IntFeasTol: with the default 1e-5, z_i = 1 - 1e-5 times big-M ~ 250 lets the epigraph
        # undercut the true terminal cost by ~1e-3 near the origin
        self.prob.solve(solver=cp.GUROBI, Threads=1, MIPGap=1e-9, IntFeasTol=1e-9)
        t = self.prob.solver_stats.solve_time
        if self.prob.status in ("infeasible", "infeasible_or_unbounded"):
            return np.inf, None, t
        return self.prob.value, self.u.value[0], t


def lq_closed_loop(policy, x0, T=50):
    x, J, xs = np.array(x0, dtype=float), 0.0, []
    for _ in range(T - 1):
        v, u = policy(x)
        if not np.isfinite(v):
            return np.inf, xs
        xs.append(x.copy())
        J += x @ x + u ** 2
        x = A_LQ @ x + B_LQ * u
    return J + x @ x, xs


def lq_lower_bound(x0, N=60):
    return Lookahead([(A_LQ, B_LQ)] * N, 1.0, None)(x0)[0]


def example_4_1(trials):
    print("=" * 100)
    print("Example 4.1 (LQ with constraints): eq. (14) equivalence and computation time")
    units = lq_units()
    for i, d in enumerate(units, 1):
        print(f"  unit {i}: L = {np.round(d['L'], 4)}, P = {np.round(d['P'], 3).tolist()}, "
              f"S has {len(d['S'][1])} facets, ellipsoid level {d['E'][1]:.3f}")
    lbs = {x0: lq_lower_bound(np.array(x0)) for x0 in X0_LQ}
    for cname, (ells, e_idx) in CONFIGS_LQ.items():
        print(f"\n  configuration: {cname}")
        probs = lq_lookaheads(units, ells, e_idx)

        def parallel(x):
            v, u, _ = best(probs, x)
            return v, u

        print(f"  {'x0':>12} | {'J~1':>6} {'J~2':>6} {'J~3':>6} {'J~4':>6} | i~ | {'J_mu~':>7} | "
              f"single-unit J_mu~i            | {'LB':>7} gap(J_mu~)")
        states = []
        for x0 in X0_LQ:
            Jt = [p(np.array(x0))[0] for p in probs]
            Jmu, xs = lq_closed_loop(parallel, x0)
            states += xs
            singles = [lq_closed_loop(lambda x, p=p: best([p], x)[:2], x0)[0] for p in probs]
            lb = lbs[x0]
            print(f"  {str(x0):>12} | " + " ".join(f"{fmt(v):>6}" for v in Jt) +
                  f" | {int(np.argmin(Jt)) + 1:>2} | {Jmu:7.3f} | " + " ".join(f"{fmt(v, 3):>7}" for v in singles) +
                  f" | {lb:7.3f} {(Jmu - lb) / Jmu:9.1e}")

        # eq. (14): min over the units equals the single problem with binary variables
        alt = Alternative(units, ells, e_idx)
        gprobs = lq_lookaheads(units, ells, e_idx, solver=cp.GUROBI)
        diffs, t_max, t_mean, t_alt = [], [], [], []
        for _ in range(trials):
            for x in states:
                vals = [p(x)[0] for p in gprobs]
                ts = [p.solve_time for p in gprobs]
                va, _, ta = alt(x)
                diffs.append(abs(min(vals) - va) / (1 + abs(va)))
                t_max.append(max(ts))
                t_mean.append(np.mean(ts))
                t_alt.append(ta)
        print(f"  eq. (14) check on {len(states)} closed-loop states: max |min_i J~i - V_alt| / (1 + V_alt) = "
              f"{max(diffs):.1e}")
        print(f"  Gurobi solve time per state, {trials} trials, 1 thread each (ms): "
              f"parallel max over units {1e3 * np.mean(t_max):.2f}, mean over units {1e3 * np.mean(t_mean):.2f}, "
              f"binary formulation {1e3 * np.mean(t_alt):.2f}; ratio max/alt {np.mean(t_max) / np.mean(t_alt):.0%}, "
              f"mean/alt {np.mean(t_mean) / np.mean(t_alt):.0%}")


# ============================================================================= main
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--example", default="all", choices=["4.1", "4.2", "4.3", "all"])
    ap.add_argument("--trials", type=int, default=10, help="timing repetitions in Example 4.1")
    args = ap.parse_args()
    if args.example in ("4.1", "all"):
        example_4_1(args.trials)
    if args.example in ("4.2", "all"):
        example_4_2()
    if args.example in ("4.3", "all"):
        example_4_3()
