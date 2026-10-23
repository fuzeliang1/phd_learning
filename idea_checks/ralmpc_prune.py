"""Robust check of sample-set pruning in RALMPC (Petrenz, Koehler, Borrelli, arXiv 2504.11261).

This is an independent Python implementation written from the paper's equations (22)-(29) and
Algorithm 1, using the mass-spring-damper data of main.m. It does not reuse the MATLAB code, whose
cost scaling and sample-set bookkeeping differ from the paper (see README.md).

Conventions follow the paper: B_p is the unit hypercube, Theta_0 = theta_bar0 + eta0 * B_p with
theta_bar0 = 0 and eta0 = 1, so k in [0.5, 1.5] and c in [0.1, 0.3]. True parameters k*=0.5, c*=0.3
correspond to theta* = (1, -1).

Sample-set strategies compared after every iteration:
  full    : keep every stored point (what the theory assumes)
  lifted  : keep only extreme points of the epigraph in the lifted space (x, s, J)
  xhull   : keep vertices of the convex hull of the (x1, x2) projection plus the initial
            trajectory, mimicking get_updateSampleSet in solve_RALMPC.m
Checks recorded at every time step:
  feasibility of the MPC problem, constraint satisfaction of the real closed loop,
  tube containment H (x_{t+1} - x_{1|t}) <= s_{1|t}  (sanity check of the implementation),
  cost decrease of Theorem 1, eq. (45):  V_{t+1} - V_t <= -l(x_t, u_t) + l_max,s
"""
import argparse
import time

import cvxpy as cp
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import ConvexHull, HalfspaceIntersection

# ----------------------------------------------------------------------------- system (main.m)
Ts, mass, k_nom, c_nom = 0.1, 1.0, 1.0, 0.2
A0 = np.array([[1.0, Ts], [-Ts * k_nom / mass, 1 - Ts * c_nom / mass]])
A1 = 0.5 * k_nom / mass * np.array([[0.0, 0.0], [Ts, 0.0]])
A2 = 0.5 * c_nom / mass * np.array([[0.0, 0.0], [0.0, Ts]])
B0 = np.array([[0.0], [Ts / mass]])
THETA_STAR = np.array([1.0, -1.0])
A_STAR = A0 + THETA_STAR[0] * A1 + THETA_STAR[1] * A2
B_STAR = B0.copy()
W_MAX = 0.2 / mass * Ts                      # |d_2| <= 0.02, d_1 = 0
D_VERTS = np.array([[0.0, W_MAX], [0.0, -W_MAX]])
Q = np.diag([1.0, 1e-2])
R = np.array([[1e-1]])
F = np.array([[1 / 4.1, 0], [-1 / 0.2, 0], [0, 1 / 5], [0, -1 / 5], [0, 0], [0, 0]])
G = np.array([[0], [0], [0], [0], [1 / 15], [-1 / 15]])
X0 = np.array([4.0, 0.0])
E_VERTS = np.array([[1, 1], [1, -1], [-1, 1], [-1, -1]], dtype=float)  # vertices of unit box B_p
ETA0 = 1.0
THETA_BAR0 = np.zeros(2)


def A_of(theta):
    return A0 + theta[0] * A1 + theta[1] * A2


def stage(x, u):
    return float(x @ Q @ x + u @ R @ u)


# ----------------------------------------------------------------------------- offline design
def design_K_P(rho=0.75):
    """LMI of get_P_K.m (Assumption 2 plus rho-contraction), evaluated at the true vertices of Theta_0."""
    n, m = 2, 1
    X = cp.Variable((n, n), symmetric=True)
    Y = cp.Variable((m, n))
    Qh, Rh = np.sqrt(Q), np.sqrt(R)
    cons = []
    for e in E_VERTS:
        A = A_of(THETA_BAR0 + ETA0 * e)
        AXBY = A @ X + B0 @ Y
        M1 = cp.bmat([[X, AXBY.T, X @ Qh, Y.T @ Rh],
                      [AXBY, X, np.zeros((n, n)), np.zeros((n, m))],
                      [Qh @ X, np.zeros((n, n)), np.eye(n), np.zeros((n, m))],
                      [Rh @ Y, np.zeros((m, n)), np.zeros((m, n)), np.eye(m)]])
        M2 = cp.bmat([[rho * X, AXBY.T], [AXBY, rho * X]])
        cons += [(M1 + M1.T) / 2 >> 0, (M2 + M2.T) / 2 >> 0]
    cp.Problem(cp.Maximize(cp.log_det(X)), cons).solve(solver=cp.CLARABEL)
    P = np.linalg.inv(X.value)
    K = Y.value @ P
    for e in E_VERTS:
        Acl = A_of(ETA0 * e) + B0 @ K
        assert np.max(np.linalg.eigvalsh(Acl.T @ P @ Acl + Q + K.T @ R @ K - P)) < 1e-6
    return K, P


def remove_redundant(Hm):
    keep = []
    for i in range(len(Hm)):
        idx = keep + list(range(i + 1, len(Hm)))  # rows kept so far plus rows not yet examined
        if not idx:
            keep.append(i)
            continue
        res = linprog(-Hm[i], A_ub=Hm[idx], b_ub=np.ones(len(idx)), bounds=(None, None), method="highs")
        if res.status != 0 or -res.fun > 1 + 1e-9:
            keep.append(i)
    return Hm[keep]


# parameter points used to shape P in get_polytopeP.m: rows of the H-representation of Theta_0,
# i.e. the axis points (+-1, 0), (0, +-1). Only the shape of P depends on this choice; the tube
# constants rho, L_B and the online tube use the true vertices E_VERTS, as required by Lemma 2.
SHAPE_POINTS = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]], dtype=float)


def design_polytope(K, lam=0.75, max_iter=100, shape_points=SHAPE_POINTS):
    """Maximal lam-contractive polytope for the symmetrised constraints, built as in get_polytopeP.m.

    Any polytope satisfying Assumption 3 is admissible; this reproduces the design choice of the
    MATLAB code (lam = 0.75 at the axis points), which keeps the tube small enough for x0 = (4, 0).
    """
    Fs, Gs = F.copy(), G.copy()
    for j in range(Fs.shape[1]):
        nz = Fs[:, j] != 0
        Fs[nz, j] = np.sign(Fs[nz, j]) * np.max(np.abs(Fs[nz, j]))
    for j in range(Gs.shape[1]):
        nz = Gs[:, j] != 0
        Gs[nz, j] = np.sign(Gs[nz, j]) * np.max(np.abs(Gs[nz, j]))
    V = Fs + Gs @ K
    V = remove_redundant(V[np.linalg.norm(V, axis=1) > 0])
    frontier = V.copy()
    AKs = [(A_of(ETA0 * e) + B0 @ K) / lam for e in shape_points]
    for _ in range(max_iter):
        add = []
        for row in np.vstack([frontier @ AK for AK in AKs]):
            res = linprog(-row, A_ub=V, b_ub=np.ones(len(V)), bounds=(None, None), method="highs")
            if -res.fun > 1 + 1e-9:
                add.append(row)
        if not add:
            break
        V = remove_redundant(np.vstack([V, np.array(add)]))
        # next frontier: newly added rows that survived the redundancy removal
        frontier = np.array([r for r in add if np.min(np.linalg.norm(V - r, axis=1)) < 1e-12])
        if len(frontier) == 0:
            break
    else:
        raise RuntimeError("contractive polytope not finitely determined; increase lam")
    H = V
    verts = HalfspaceIntersection(np.hstack([H, -np.ones((len(H), 1))]), np.zeros(2)).intersections
    return H, verts


class Design:
    def __init__(self):
        self.K, self.P = design_K_P()
        self.H, self.PV = design_polytope(self.K)
        self.rho0 = self.rho_of(THETA_BAR0)
        # L_B (10): max over vertices e_l of B_p and x in P of H_i D(x, Kx) e_l
        self.L_B = max(float(np.max(self.H @ ((e[0] * A1 + e[1] * A2) @ self.PV.T))) for e in E_VERTS)
        self.d_bar = float(np.max(self.H @ D_VERTS.T))
        self.c = np.max((F + G @ self.K) @ self.PV.T, axis=1)
        self.W = self.uncertainty_pieces()
        assert self.rho0 + ETA0 * self.L_B < 1, "Assumption 3 violated"
        self.s_steady = self.d_bar / (1 - (self.rho0 + ETA0 * self.L_B))
        self.L_cost = self.lipschitz_cost()
        self.l_max_s = self.L_cost * self.s_steady

    def uncertainty_pieces(self):
        """Non-redundant rows g of w_eta(x)/eta = max_{i,l} H_i D(x, u) e_l = max_g g x.

        Here B_1 = B_2 = 0, so D does not depend on u and every piece is linear in x. A piece is
        redundant for the max over all x iff its coefficient vector lies in the convex hull of the
        others, so keeping the vertices of that hull is an exact reformulation (112 rows -> a few)."""
        g = np.array([self.H[i] @ (e[0] * A1 + e[1] * A2) for e in E_VERTS for i in range(len(self.H))])
        hull = ConvexHull(g)
        return g[np.sort(hull.vertices)]

    def rho_of(self, theta_bar):
        Acl = A_of(theta_bar) + B0 @ self.K
        return float(np.max(self.H @ Acl @ self.PV.T))

    def lipschitz_cost(self):
        """Constant of (19) via the mean value theorem: max over vertices of Z (in (x, v)) and of P u -P."""
        FK = F + G @ self.K
        Hz = np.hstack([FK, G])
        res = linprog(np.zeros(3), A_ub=Hz, b_ub=np.ones(len(Hz)) * 0.999, bounds=(None, None), method="highs")
        zv = HalfspaceIntersection(np.hstack([Hz, -np.ones((len(Hz), 1))]), res.x).intersections
        best = 0.0
        for xv in zv:
            x, v = xv[:2], xv[2:]
            grad = 2 * Q @ x + 2 * self.K.T @ R @ (self.K @ x + v)
            for p in np.vstack([self.PV, -self.PV]):
                best = max(best, abs(float(grad @ p)))
        return best

    def l_max(self, x, v, s):
        return stage(x, self.K @ x + v) + self.L_cost * s


# ----------------------------------------------------------------------------- online problem
class RALMPCProblem:
    """Problem (22) for a fixed sample set; theta_bar, eta, rho and the state are parameters."""

    LAM_REG = 1e-7  # tiny ridge on lambda: the objective is linear in lambda, so without it the optimal
    #                 lambda is massively non-unique and interior-point solvers lose accuracy

    def __init__(self, dsg, N, Z, J, d_margin=0.0, c_margin=0.0):
        """d_margin inflates d_bar and c_margin tightens the constraints; both are used only for the
        initial solution, so that the stored tube has slack and the online problems have an interior."""
        self.dsg, self.N = dsg, N
        n_pts = Z.shape[0]
        self.x0 = cp.Parameter(2)
        self.Acl = cp.Parameter((2, 2))
        self.eta = cp.Parameter(nonneg=True)
        self.rho = cp.Parameter(nonneg=True)
        self.x = cp.Variable((N + 1, 2))
        self.v = cp.Variable((N, 1))
        self.s = cp.Variable(N + 1)
        self.lam = cp.Variable(n_pts, nonneg=True)
        sig = cp.Variable()
        H, K = dsg.H, dsg.K
        cons = [self.x[0] == self.x0, self.s[0] == 0]
        obj = 0
        for k in range(N):
            xk, vk, sk = self.x[k], self.v[k], self.s[k]
            uk = K @ xk + vk
            cons += [self.x[k + 1] == self.Acl @ xk + B0 @ vk]
            # s_{k+1} >= rho s_k + d_bar + eta L_B s_k + eta max_{i,l} H_i D(x_k, u_k) e_l
            Dterm = dsg.W @ xk  # exact reduction of max_{i,l} H_i D(x_k, u_k) e_l, see uncertainty_pieces
            cons += [self.s[k + 1] >= self.rho * sk + dsg.d_bar * (1 + d_margin)
                     + self.eta * dsg.L_B * sk + self.eta * Dterm]
            cons += [(F + G @ K) @ xk + G @ vk + dsg.c * sk <= 1 - c_margin]
            obj += cp.quad_form(xk, Q) + cp.quad_form(uk, R) + dsg.L_cost * sk
        xt = Z[:, :2].T @ self.lam
        st = Z[:, 2] @ self.lam
        cons += [cp.sum(self.lam) == 1, H @ (self.x[N] - xt) <= sig, sig <= st - self.s[N]]
        obj += J @ self.lam + self.LAM_REG * cp.sum_squares(self.lam)
        self.prob = cp.Problem(cp.Minimize(obj), cons)
        assert self.prob.is_dcp(dpp=True)

    # Clarabel first; an "inaccurate" Clarabel result is accepted only if its constraint violation is
    # below FEAS_TOL. On failure fall back to HiGHS (active-set QP) and PIQP (proximal interior point),
    # which cope better with the degenerate sample sets that accumulate near the origin, then SCS and
    # OSQP. The rank used is reported as "fallback".
    SOLVERS = [(cp.CLARABEL, dict(max_iter=400)),
               (cp.HIGHS, {}),
               (cp.PIQP, {}),
               (cp.SCS, dict(eps=1e-9, max_iters=200000)),
               (cp.OSQP, dict(eps_abs=1e-9, eps_rel=1e-9, polish=True, max_iter=400000))]
    FEAS_TOL = 1e-7

    def max_violation(self):
        return max(float(np.max(np.atleast_1d(c.violation()))) for c in self.prob.constraints)

    def solve(self, x, theta_bar, eta, rho):
        self.x0.value = x
        self.Acl.value = A_of(theta_bar) + B0 @ self.dsg.K
        self.eta.value, self.rho.value = eta, rho
        for rank, (solver, opts) in enumerate(self.SOLVERS):
            try:
                self.prob.solve(solver=solver, **opts)
            except cp.SolverError:
                continue
            if self.prob.status == "optimal" or (
                    self.prob.status == "optimal_inaccurate" and self.max_violation() < self.FEAS_TOL):
                break
        else:
            return None
        return dict(V=self.prob.value, x=self.x.value.copy(), v=self.v.value.copy(),
                    s=self.s.value.copy(), fallback=rank, Qterm=float(self.prob.value - sum(
                        self.dsg.l_max(self.x.value[k], self.v.value[k], self.s.value[k]) for k in range(self.N))))


# ----------------------------------------------------------------------------- set membership (Alg. 1)
class Hypercube:
    def __init__(self, M=10):
        self.M = M
        self.theta_bar, self.eta = THETA_BAR0.copy(), ETA0
        self.window = []

    def update(self, x_prev, u_prev, x_now):
        r = x_now - A0 @ x_prev - B0 @ u_prev
        Dm = np.column_stack([A1 @ x_prev, A2 @ x_prev])
        # x_now - A_theta x - B u in D  <=>  d = r - Dm theta with d_1 = 0, |d_2| <= W_MAX
        Ad = np.array([[0, 1], [0, -1], [1, 0], [-1, 0]], dtype=float)
        bd = np.array([W_MAX, W_MAX, 0, 0])
        self.window = (self.window + [(-Ad @ Dm, bd - Ad @ r)])[-self.M:]
        A_ub = [np.vstack([np.eye(2), -np.eye(2)])]
        b_ub = [np.concatenate([self.theta_bar + self.eta, -(self.theta_bar - self.eta)])]
        for Aw, bw in self.window:
            A_ub.append(Aw)
            b_ub.append(bw)
        A_ub, b_ub = np.vstack(A_ub), np.concatenate(b_ub) + 1e-12
        lo, hi = np.zeros(2), np.zeros(2)
        for i in range(2):
            e = np.eye(2)[i]
            lo[i] = linprog(e, A_ub=A_ub, b_ub=b_ub, bounds=(None, None), method="highs").fun
            hi[i] = -linprog(-e, A_ub=A_ub, b_ub=b_ub, bounds=(None, None), method="highs").fun
        eta_new = 0.5 * float(np.max(hi - lo))
        tb_new = 0.5 * (lo + hi)
        rad = self.eta - eta_new  # projection keeps theta_bar_t - theta_bar_{t-1} in (eta_{t-1}-eta_t) B_p
        tb_new = np.clip(tb_new, self.theta_bar - rad, self.theta_bar + rad)
        self.theta_bar, self.eta = tb_new, eta_new


# ----------------------------------------------------------------------------- initial solution (Ass. 4)
def steady_chain(dsg, s_bar, tol=1e-12):
    """Points (0, s_j) with s_{j+1} = (rho0 + eta0 L_B) s_j + d_bar decaying from s_bar to s_steady.

    With x = 0 and v = 0 every level s >= s_steady is robustly invariant, but a tube that has grown
    above s_steady only reaches it asymptotically, so the terminal point of Assumption 4 is replaced by
    this chain. Costs follow (27): J_j = sum_{i >= j} (l_max(0, 0, s_i) - l_max,s).
    """
    a = dsg.rho0 + ETA0 * dsg.L_B
    s = [s_bar]
    while s[-1] - dsg.s_steady > tol:
        s.append(a * s[-1] + dsg.d_bar)
    s[-1] = dsg.s_steady
    s = np.array(s)
    costs = np.cumsum((dsg.L_cost * (s - dsg.s_steady))[::-1])[::-1]
    return np.column_stack([np.zeros((len(s), 2)), s]), costs


def initial_solution(dsg, Nbar=40):
    """Assumption 4: a feasible long-horizon tube trajectory from X0, plus the steady-state chain."""
    s_bar = 0.9 / float(np.max(dsg.c))  # largest level that keeps a margin to the tightened constraints
    Zc, Jc = steady_chain(dsg, s_bar)
    sol = RALMPCProblem(dsg, Nbar, Zc, Jc, d_margin=0.01).solve(X0, THETA_BAR0, ETA0, dsg.rho0)
    assert sol is not None, "no initial solution, increase Nbar"
    pts, costs = [], []
    Jk = sol["Qterm"]
    for k in range(Nbar - 1, -1, -1):
        Jk = dsg.l_max(sol["x"][k], sol["v"][k], sol["s"][k]) - dsg.l_max_s + Jk
        pts.append([*sol["x"][k], sol["s"][k]])
        costs.append(Jk)
    Z0 = np.vstack([np.array(pts), Zc])
    J0 = np.concatenate([np.array(costs), Jc])
    return Z0, J0  # last row is the steady-state point (0, 0, s_steady)


# ----------------------------------------------------------------------------- pruning
def dominated(Zc, Jc, z, Jz, tol=1e-9):
    A_eq = np.vstack([Zc.T, np.ones(len(Zc))])
    b_eq = np.concatenate([z, [1.0]])
    res = linprog(Jc, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method="highs")
    return res.status == 0 and res.fun <= Jz + tol * max(1.0, abs(Jz))


def prune_lifted(Z, J, protect):
    """Keep extreme points of the epigraph of the lower convex envelope of {(z_i, J_i)}."""
    pts = np.hstack([Z, J[:, None]])
    cand = set(protect)
    try:
        scale = np.maximum(np.ptp(pts, axis=0), 1e-9)
        hull = ConvexHull(pts / scale)
        lower = hull.equations[:, 3] < -1e-12
        cand |= set(np.unique(hull.simplices[lower]).tolist())
    except Exception:
        cand = set(range(len(Z)))
    # verify every dropped point is dominated by the kept ones; add it back otherwise.
    # one pass suffices: adding candidates can only make the remaining points easier to dominate
    for i in range(len(Z)):
        if i in cand:
            continue
        idx = sorted(cand)
        if not dominated(Z[idx], J[idx], Z[i], J[i]):
            cand.add(i)
    cand = sorted(cand)
    return Z[cand], J[cand]


def dedup(Z, J, dec=10):
    """Merge numerically identical points, keeping the cheapest copy (lossless: same z, lower J)."""
    order = np.argsort(J, kind="stable")
    _, first = np.unique(np.round(Z[order], dec), axis=0, return_index=True)
    keep = np.sort(order[first])
    return Z[keep], J[keep]


def prune_xhull(Z, J, Z0, J0):
    hull = ConvexHull(Z[:, :2])
    allZ = np.vstack([Z[hull.vertices], Z0])
    allJ = np.concatenate([J[hull.vertices], J0])
    _, idx = np.unique(np.round(allZ, 12), axis=0, return_index=True)
    idx = np.sort(idx)
    return allZ[idx], allJ[idx]


# ----------------------------------------------------------------------------- experiment
def print_row(r):
    print(f"{r['h']:2d} {r['pts']:5d} {r['cost']:9.4f} {r['infeasible']:5d} {r['cons_viol']:4d} "
          f"{r['tube_viol']:4d} {r['dec_viol']:6d}   {r['dec_max']:.2e} {r['eta']:.4f} {1e3 * r['t_solve']:8.1f} "
          f"{r['fallback']:6d}", flush=True)


def run(dsg, strategy, N, iters, T, adaptive, disturbance, seed=0, verbose=False):
    rng = np.random.default_rng(seed)
    Z0, J0 = initial_solution(dsg)
    Z, J = Z0.copy(), J0.copy()
    hc = Hypercube()
    rows = []
    for h in range(iters):
        prob = RALMPCProblem(dsg, N, Z, J)
        x = X0.copy()
        V_prev = l_prev = None
        log = dict(h=h, pts=len(Z), infeasible=0, cons_viol=0, tube_viol=0, dec_viol=0, dec_max=0.0, cost=0.0,
                   fallback=0)
        new_Z, new_J = [], []
        t_solve = 0.0
        x_prev = u_prev = None
        for t in range(T):
            if adaptive and t > 0:
                hc.update(x_prev, u_prev, x)
            rho = dsg.rho_of(hc.theta_bar) if adaptive else dsg.rho0
            t0 = time.perf_counter()
            sol = prob.solve(x, hc.theta_bar, hc.eta, rho)
            t_solve += time.perf_counter() - t0
            if sol is None:
                log["infeasible"] += 1
                break
            log["fallback"] += int(sol["fallback"] > 0)
            u = dsg.K @ x + sol["v"][0]
            if np.any(F @ x + G @ u > 1 + 1e-9):
                log["cons_viol"] += 1
            if V_prev is not None:
                excess = sol["V"] - V_prev + l_prev - dsg.l_max_s
                if excess > 1e-6 * max(1.0, abs(V_prev)):
                    log["dec_viol"] += 1
                    log["dec_max"] = max(log["dec_max"], excess)
            # worst-case cost-to-go along the prediction, eq. (28)
            Jk = sol["Qterm"]
            for k in range(N - 1, -1, -1):
                Jk = dsg.l_max(sol["x"][k], sol["v"][k], sol["s"][k]) - dsg.l_max_s + Jk
                new_Z.append([*sol["x"][k], sol["s"][k]])
                new_J.append(Jk)
            d = D_VERTS[0] * 0.5 if disturbance == "constant" else np.array([0.0, rng.uniform(-W_MAX, W_MAX)])
            x_next = A_STAR @ x + B_STAR @ u + d
            if np.any(dsg.H @ (x_next - sol["x"][1]) > sol["s"][1] + 1e-7):
                log["tube_viol"] += 1
            log["cost"] += stage(x, u)
            V_prev, l_prev = sol["V"], stage(x, u)
            x_prev, u_prev, x = x, u, x_next
        log["t_solve"] = t_solve / max(1, t + 1)
        log["eta"] = hc.eta
        rows.append(log)
        if verbose:
            print_row(log)
        if log["infeasible"]:
            break
        Z = np.vstack([Z, np.array(new_Z)])
        J = np.concatenate([J, np.array(new_J)])
        Z, J = dedup(Z, J)
        if strategy == "lifted":
            protect = [int(np.argmin(np.linalg.norm(Z - Z0[-1], axis=1)))]
            Z, J = prune_lifted(Z, J, protect)
        elif strategy == "xhull":
            Z, J = prune_xhull(Z, J, Z0, J0)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=8)
    ap.add_argument("--iters", type=int, default=8)
    ap.add_argument("--T", type=int, default=60)
    ap.add_argument("--adaptive", action="store_true")
    ap.add_argument("--disturbance", choices=["constant", "random"], default="constant")
    ap.add_argument("--strategies", default="full,lifted,xhull")
    args = ap.parse_args()
    dsg = Design()
    print(f"K={np.round(dsg.K, 4)}, |P vertices|={len(dsg.PV)}, rho0={dsg.rho0:.4f}, L_B={dsg.L_B:.4f}, "
          f"d_bar={dsg.d_bar:.4f}, s_steady={dsg.s_steady:.4f}, L_cost={dsg.L_cost:.3f}, l_max_s={dsg.l_max_s:.4f}")
    for strat in args.strategies.split(","):
        print(f"\n== {strat}  (N={args.N}, adaptive={args.adaptive}, disturbance={args.disturbance})")
        print(" h   pts   cost      infeas cons tube dec_viol dec_max    eta    solve_ms fallback", flush=True)
        run(dsg, strat, args.N, args.iters, args.T, args.adaptive, args.disturbance, verbose=True)
