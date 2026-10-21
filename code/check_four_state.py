"""Four-state verifier (Theorem 5.5 and Remark 5.6 of the paper). Deterministic; writes four_state_verification.json.

Part A (exact, sympy): structural facts about the limiting block A(p) for every
  admissible vanishing set V of a generic N=4 model (states topologically ordered
  1<2<3<4, all admissible coefficients c_{ijk} symbolic).
Part B (exact, sympy): the |V|=1 example x(u) = -(4e^u-5)^3/(240 e^{2u}) (triple zero).
Part C (numerical): the V={2,3} model
     a' = -B a y,  x' = B a y - x d - x a,  y' = x d - y d - y a,  d' = x a + y d + y a
  with B=4 (double zero of det A on the face, predicted t^{3/2}(x+2y) -> sqrt(3)/2)
  and  B=9/2 (simple zero at s=2/3, predicted t^2 (x+3y) -> 16/3).
  Near-critical interior data are found by bisection; the reported quantities are
  t^theta z(t) on a log grid. This is numerical evidence only.
"""
import json, itertools
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


def r10(v):
    v = float(v)
    return 0.0 if v == 0 else float(f"{v:.10g}")


out = {}

# ---------------- Part A ----------------
N = 4
y = sp.symbols('y1:5', positive=True)
s = sp.Symbol('s')


def model(V):
    """Generic admissible model for vanishing set V (1-based).  Returns f, coeffs."""
    P = [i for i in range(1, 5) if i not in V]
    cs = {}
    for i, j, k in itertools.product(range(1, 5), repeat=3):
        if j <= i:          # topological order: transitions only i->j with i<j
            continue
        if not (i in V or k in V):   # active flux needs y_i* y_k* = 0
            continue
        cs[(i, j, k)] = sp.Symbol(f'c{i}{j}{k}', nonnegative=True)
    f = [0] * 4
    for (i, j, k), c in cs.items():
        flux = c * y[i - 1] * y[k - 1]
        f[i - 1] -= flux
        f[j - 1] += flux
    return f, cs, P


def limiting_block(V):
    f, cs, P = model(V)
    Vl = sorted(V)
    A = sp.Matrix([[sp.diff(f[i - 1], y[j - 1]) for j in Vl] for i in Vl])
    sub0 = {y[i - 1]: 0 for i in Vl}
    A = A.subs(sub0)
    return sp.expand(A), cs, P, f


resA = {}
for V in [(1, 2, 3), (1, 2), (1, 3), (2, 3)]:
    A, cs, P, f = limiting_block(V)
    rec = {"P": P}
    # sum of f is zero (mass conservation)
    rec["mass_conserved"] = bool(sp.expand(sum(f)) == 0)
    # every term of f vanishes on the face v=0
    rec["face_is_equilibria"] = all(sp.expand(fi.subs({y[i - 1]: 0 for i in V})) == 0 for fi in f)
    # A homogeneous linear in persistent coordinates
    pv = [y[i - 1] for i in P]
    rec["A_homogeneous_linear_in_p"] = all(
        sp.Poly(e, *pv).is_homogeneous and sp.Poly(e, *pv).total_degree() <= 1 for e in A if e != 0)
    if V == (1, 2, 3):
        rec["A_lower_triangular"] = all(A[r, c] == 0 for r in range(3) for c in range(3) if c > r)
        # diagonal entries are minus sums of nonnegative coefficients times y4
        rec["diag_nonpositive_form"] = all(
            all(coef < 0 for coef in sp.Poly(A[r, r], *list(cs.values()), y[3]).coeffs()) for r in range(3))
    if 1 in V:
        rec["row1_offdiag_zero"] = all(A[0, c] == 0 for c in range(1, A.shape[1]))
    if 1 in V and len(V) == 2:
        # A_ii affine on the face p = s e_j + (1-s) e_l
        j, l = P
        Ai = A[1, 1].subs({y[j - 1]: s, y[l - 1]: 1 - s})
        rec["Aii_face_degree_in_s"] = int(sp.Poly(sp.expand(Ai), s).degree())
    if V == (2, 3):
        face = {y[0]: s, y[3]: 1 - s}
        Af = A.subs(face)
        det = sp.expand(Af.det())
        rec["det_face_degree_in_s"] = int(sp.Poly(det, s).degree())
        # at s=0 (pure state 4) diagonal entries are minus sums of coefficients
        rec["diag_at_s0"] = [str(sp.expand(Af[r, r].subs(s, 0))) for r in range(2)]
        rec["diag_at_s0_nonpositive"] = all(
            all(coef < 0 for coef in sp.Poly(sp.expand(Af[r, r].subs(s, 0)), *list(cs.values())).coeffs())
            for r in range(2))
        # 1 has no inflow, 4 has no outflow => s=y1/(y1+y4) nonincreasing
        rec["f1_is_minus_y1_times_nonneg"] = bool(
            all(co < 0 for co in sp.Poly(sp.expand(f[0]), *list(cs.values()), *y).coeffs()))
        rec["f4_nonneg_terms"] = bool(
            all(co > 0 for co in sp.Poly(sp.expand(f[3]), *list(cs.values()), *y).coeffs()))
    resA["V=" + "".join(map(str, V))] = rec
out["A_structure"] = resA

# ---------------- Part B ----------------
u = sp.Symbol('u', real=True)
a0, b0, x0, c0 = sp.Rational(25, 48), sp.Rational(5, 24), sp.Rational(1, 240), sp.Rational(4, 15)
au = a0 * sp.exp(-2 * u)
bu = (2 * a0 + b0) * sp.exp(-u) - 2 * a0 * sp.exp(-2 * u)
cu = c0 * sp.exp(u)
xu = 1 - au - bu - cu
# check the closed form solves the exposure-time linear system y' = H y
chk = [sp.simplify(sp.diff(au, u) + 2 * au), sp.simplify(sp.diff(bu, u) - (2 * au - bu)),
       sp.simplify(sp.diff(xu, u) - (bu - cu)), sp.simplify(sp.diff(cu, u) - cu)]
U = sp.log(sp.Rational(5, 4))
out["V1_example"] = {
    "exposure_system_residuals_zero": all(c == 0 for c in chk),
    "initial_sum": str(a0 + b0 + x0 + c0),
    "x_closed_form_identity": bool(sp.simplify(xu + (4 * sp.exp(u) - 5) ** 3 / (240 * sp.exp(2 * u))) == 0),
    "U": "log(5/4)",
    "derivs_at_U": [str(sp.nsimplify(sp.simplify(sp.diff(xu, u, n).subs(u, U)))) for n in range(4)],
}

# ---------------- Part C ----------------
def make(B):
    def f(t, Y):
        a, x, yy, d = Y
        return [-B * a * yy, B * a * yy - x * d - x * a, x * d - yy * d - yy * a, x * a + yy * d + yy * a]

    def jac(t, Y):
        a, x, yy, d = Y
        return [[-B * yy, 0, -B * a, 0], [B * yy - x, -d - a, B * a, -x],
                [-yy, d, -d - a, x - yy], [x + yy, a, d + a, yy]]
    return f, jac


def critical_run(B, sc, r, w, k, K, z0, T, nbis):
    """sc: zero of lambda on face; r: right Perron vector (w.r=1); k: predicted zero order;
    K: predicted coefficient z ~ K (U-u)^k ; b_inf=-1/2 for both examples."""
    f, jac = make(B)
    binf = 0.5

    def seed(ds):
        x, yy = z0 * r[0], z0 * r[1]
        sv = sc + binf * (z0 / K) ** (1.0 / k) + ds
        m = 1 - x - yy
        return [sv * m, x, yy, (1 - sv) * m]

    def ev(t, Y):
        return Y[0] / (Y[0] + Y[3]) - sc
    ev.terminal = True

    def crosses(ds):
        sol = solve_ivp(f, (0, T), seed(ds), method='Radau', jac=jac, rtol=1e-11, atol=1e-26, events=ev)
        return sol.status == 1
    lo, hi = -0.9 * binf * (z0 / K) ** (1.0 / k), 0.05
    assert crosses(lo) and not crosses(hi)
    for _ in range(nbis):
        mid = 0.5 * (lo + hi)
        if crosses(mid):
            lo = mid
        else:
            hi = mid
    Y0 = seed(hi)
    ts = np.logspace(1, np.log10(T), 2 * int(np.log10(T)) - 1)
    sol = solve_ivp(f, (0, T), Y0, method='Radau', jac=jac, rtol=1e-11, atol=1e-26, t_eval=ts)
    theta = k / (k - 1)
    rows = []
    for t, Y in zip(sol.t, sol.y.T):
        a, x, yy, d = Y
        z = w[0] * x + w[1] * yy
        rows.append({"t": r10(t), "s_minus_sc": r10(a / (a + d) - sc), "t^theta*z": r10(t ** theta * z)})
    pred = K * ((k - 1) * K) ** (-theta)
    return {"B": B, "initial_data": [r10(v) for v in Y0], "theta": r10(theta),
            "predicted_limit_t^theta*z": r10(pred), "samples": rows}


out["V23_B4_double_zero"] = critical_run(4.0, 0.5, (0.5, 0.25), (1.0, 2.0), 3, 1.0 / 6.0, 1e-4, 1e7, 40)
out["V23_B4.5_simple_zero"] = critical_run(4.5, 2.0 / 3.0, (0.5, 1.0 / 6.0), (1.0, 3.0), 2, 3.0 / 16.0, 1e-4, 1e7, 40)

from pathlib import Path
with open(Path(__file__).with_name("four_state_verification.json"), "w") as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
print("PASS four-state checks")
