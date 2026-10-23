"""Toy check of 'certified' parallel MPC with learned terminal costs.

Scalar example 1.1 of Li et al., Parallel MPC (TAC): x+ = 2x + u, |u| <= 1,
g(x,u) = x^2 + u^2 for |x| <= 2, infinite otherwise.

Units (lookahead length, terminal function):
  A: l=2, J = indicator{x=0}               (certified, from the paper)
  B: l=1, J = 13/3 x^2 on |x| <= 2/3        (certified, from the paper)
  C: l=2, J = learned value with a spurious 'hole' near x=c (not in the L-region there)
  D: l=2, J = learned value close to J* (small multiplicative error)

Policies compared:
  paper   : parallel MPC over {A, B}
  naive   : parallel MPC over {A, B, C, D}, pick the smallest predicted cost
  cert    : same units, but a learned unit is admissible only if its terminal function passes
            the one-step test  T J(x_N) <= J(x_N) - (1 - alpha) g(x_N, u_N*)  at the predicted
            terminal state; the shifted previous plan is always kept as a candidate.
"""
import numpy as np

U = np.linspace(-1, 1, 401)
XMAX = 2.0
ALPHA = 0.5
BREAK_TOL = 1e-6  # decrease violations below this are round-off near the origin, not breaks


def g(x, u):
    x = np.asarray(x, dtype=float)
    return np.where(np.abs(x) <= XMAX, x ** 2 + u ** 2, np.inf)


def f(x, u):
    return 2 * x + u


# ---------- reference J* by value iteration on a grid ----------
XG = np.linspace(-1.2, 1.2, 2401)


def value_iteration(iters=200):
    # start from unit B's terminal function (an L-region upper bound), so T^k J decreases to J*
    J = np.where(np.abs(XG) <= 2 / 3, 13 / 3 * XG ** 2, 1e3)
    UF = np.linspace(-1, 1, 2001)  # fine grid keeps 2x+u on the state grid
    X, Uu = np.meshgrid(XG, UF, indexing="ij")
    Xn = f(X, Uu)
    stage = g(X, Uu)
    for _ in range(iters):
        Jn = np.interp(Xn, XG, J, left=np.inf, right=np.inf)
        J = np.min(stage + Jn, axis=1)
    return J


JSTAR = value_iteration()


def jstar(x):
    return np.interp(x, XG, JSTAR, left=np.inf, right=np.inf)


C_HOLE = 0.45


def J_A(x):
    return np.where(np.abs(x) < 3e-3, 0.0, np.inf)


def J_B(x):
    return np.where(np.abs(x) <= 2 / 3, 13 / 3 * np.asarray(x) ** 2, np.inf)


def J_C(x):  # learned value with a spurious minimum near x = C_HOLE
    return jstar(x) * (1 - 0.97 * np.exp(-((np.asarray(x) - C_HOLE) / 0.08) ** 2))


def J_D(x):  # good learned value: J* with a smooth 15% error
    return jstar(x) * (1 + 0.15 * np.sin(5 * np.asarray(x)))


UNITS = {"A": (2, J_A, True), "B": (1, J_B, True), "C": (2, J_C, False), "D": (2, J_D, False)}


def solve_unit(x, l, J):
    """Brute-force l-step lookahead (l in {1,2}); returns value, first input, terminal state, plan."""
    if l == 1:
        xn = f(x, U)
        v = g(x, U) + J(xn)
        k = int(np.argmin(v))
        return v[k], U[k], xn[k], [U[k]]
    U1, U2 = np.meshgrid(U, U, indexing="ij")
    x1 = f(x, U1)
    x2 = f(x1, U2)
    v = g(x, U1) + g(x1, U2) + J(x2)
    k = np.unravel_index(int(np.argmin(v)), v.shape)
    return v[k], U1[k], x2[k], [U1[k], U2[k]]


def one_step(xN, J):
    v = g(xN, U) + J(f(xN, U))
    k = int(np.argmin(v))
    return v[k], U[k]


def terminal_test(J, xN, stage_now):
    """Relaxed L-region test at the predicted terminal state: T J(xN) - J(xN) <= (1-alpha) g(x_t, u_t)."""
    TJ, uN = one_step(xN, J)
    return bool(TJ - J(xN) <= (1 - ALPHA) * stage_now + 1e-9), uN


def simulate(x0, units, certified, T=40, return_stats=False, lock=True, log=None):
    """Closed-loop cost from x0 (plus J* of the final state).

    lock: after a certificate break, use offline-certified units only.
    log:  if a list, one (x_t, unit, v_t, g_t) tuple is appended per step.
    """
    x, total, prev = x0, 0.0, None  # prev = (name, terminal function, remaining plan)
    locked, breaks, V_prev, g_prev = False, 0, None, None
    for _ in range(T):
        cands = []
        for name in units:
            l, J, offline_certified = UNITS[name]
            if certified and locked and not offline_certified:
                continue
            v, u0, xN, plan = solve_unit(x, l, J)
            if not np.isfinite(v):
                continue
            if certified and not offline_certified and not terminal_test(J, xN, float(g(x, u0)))[0]:
                continue
            cands.append((v, u0, name, J, plan))
        if certified and prev is not None and not locked and not UNITS[prev[0]][2]:
            # shifted previous plan of a learned unit, extended by its one-step minimizer;
            # it is admissible only if its new terminal state passes the same test
            name_p, Jp, plan_p = prev
            xs = x
            for u in plan_p[1:]:
                xs = f(xs, u)
            _, uN = one_step(xs, Jp)
            rest = plan_p[1:] + [uN]
            xs, v = x, 0.0
            for u in rest:
                v += g(xs, u)
                xs = f(xs, u)
            v += Jp(xs)
            if np.isfinite(v) and terminal_test(Jp, xs, float(g(x, rest[0])))[0]:
                cands.append((v, rest[0], name_p, Jp, rest))
        if not cands:
            return (np.inf, breaks, locked) if return_stats else np.inf
        v, u, name, J, plan = min(cands, key=lambda c: c[0])
        if certified and V_prev is not None and v > V_prev - ALPHA * g_prev + BREAK_TOL:
            breaks += 1      # certificate chain broken: from now on use offline-certified units only
            locked = lock
        V_prev, g_prev = float(v), float(g(x, u))
        if log is not None:
            log.append((float(x), name, V_prev, g_prev))
        prev = (name, J, plan)
        total += float(g(x, u))
        x = f(x, u)
        if not np.isfinite(total):
            return (np.inf, breaks, locked) if return_stats else np.inf
    out = total + float(jstar(x))  # tail cost from where we stopped
    return (out, breaks, locked) if return_stats else out


if __name__ == "__main__":
    x0s = np.round(np.linspace(-0.9, 0.9, 19), 3)
    rows = []
    for x0 in x0s:
        rows.append((x0, float(jstar(x0)),
                     simulate(x0, ["A", "B"], certified=False),
                     simulate(x0, ["A", "B", "C", "D"], certified=False),
                     simulate(x0, ["A", "B", "C", "D"], certified=True)))
    print(f"{'x0':>6} {'J*':>9} {'paper{A,B}':>11} {'naive{ABCD}':>12} {'cert{ABCD}':>11}")
    for r in rows:
        print(f"{r[0]:6.2f} {r[1]:9.4f} {r[2]:11.4f} {r[3]:12.4f} {r[4]:11.4f}")
    arr = np.array([r[1:] for r in rows])
    fin = np.all(np.isfinite(arr), axis=1)
    fin = fin & (arr[:, 0] > 1e-9)
    gap = lambda k: np.mean(arr[fin, k] / arr[fin, 0] - 1) * 100
    print(f"\nmean suboptimality over x0 where all finite: paper {gap(1):.2f}%, naive {gap(2):.2f}%, cert {gap(3):.2f}%")
    print(f"worst case ratio J/J*: paper {np.max(arr[fin,1]/arr[fin,0]):.3f}, naive {np.max(arr[fin,2]/arr[fin,0]):.3f}, cert {np.max(arr[fin,3]/arr[fin,0]):.3f}")
    stats = [simulate(x0, ["A", "B", "C", "D"], certified=True, return_stats=True) for x0 in x0s]
    print("certificate breaks per x0 (cert):", [s[1] for s in stats])
    feas = lambda k: int(np.sum(np.isfinite(arr[:, k])))
    print(f"feasible initial states out of {len(x0s)}: paper {feas(1)}, naive {feas(2)}, cert {feas(3)}")
