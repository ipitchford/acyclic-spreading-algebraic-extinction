"""Exact and high-precision checks for the derivation of the counterexample,
the SWIR example, the support and rate criteria, the slowing law and
monotonicity.

Writes code/extensions_verification.json.  Deterministic: no timings and no
version stamps; floats are printed with a fixed number of significant digits.
"""
import json
from pathlib import Path
import sympy as sp
import mpmath as mp

mp.mp.dps = 60
OUT = Path(__file__).with_name("extensions_verification.json")
R = {"checks": {}}


def fl(v, n=12):
    return mp.nstr(mp.mpf(sp.N(v, n + 10)) if not isinstance(v, mp.mpf) else v, n)


def ok(name, cond):
    R["checks"][name] = bool(cond)
    assert cond, name


q, rho, t, u = sp.symbols("q rho t u", real=True)

# ---------------------------------------------------------------- Section 3
# Derivation of the counterexample: impose a triple root on the exposure cubic.
s0, w0, r0 = sp.symbols("s0 w0 r0", positive=True)
Q = s0 * q**3 - (w0 + 2 * s0) * q**2 + q - r0          # rates (2,1,1)
sol = sp.solve(sp.Poly(Q - s0 * (q - rho) ** 3, q).coeffs(), [s0, w0, r0], dict=True)[0]
ok("triple_root_data_formula",
   sp.simplify(sol[s0] - 1 / (3 * rho**2)) == 0
   and sp.simplify(sol[w0] - (3 * rho - 2) / (3 * rho**2)) == 0
   and sp.simplify(sol[r0] - rho / 3) == 0)
base = {k: sp.nsimplify(vv.subs(rho, sp.Rational(4, 5))) for k, vv in sol.items()}
x0 = 1 - base[s0] - base[w0] - base[r0]
ok("base_initial_data", (base[s0], base[w0], x0, base[r0]) ==
   (sp.Rational(25, 48), sp.Rational(5, 24), sp.Rational(1, 240), sp.Rational(4, 15)))
th = (16 + sp.Rational(2, 3) * t) ** sp.Rational(-1, 2)
ze = 1 + th
sf, wf, xf, rf = ze**2 / 3, ze * (3 - 2 * ze) / 3, th**3 / (3 * ze), 1 / (3 * ze)
res = [sp.diff(sf, t) + 2 * sf * xf, sp.diff(wf, t) - 2 * sf * xf + wf * xf,
       sp.diff(xf, t) - (wf - rf) * xf, sp.diff(rf, t) - rf * xf, sf + wf + xf + rf - 1]
ok("base_closed_form_residuals", all(sp.simplify(e) == 0 for e in res))
ok("base_tail_constant", sp.simplify(sp.limit(t**sp.Rational(3, 2) * xf, t, sp.oo)
                                     - sp.sqrt(sp.Rational(3, 8))) == 0)
Z = sp.symbols("zeta", positive=True)
Reff = Z**2 * (3 - 2 * Z)                       # w/r along the trajectory
ok("Reff_identity", sp.expand(1 - Reff - (Z - 1) ** 2 * (2 * Z + 1)) == 0)
ok("Reff_asymptotic", sp.limit(t * (1 - Reff.subs(Z, ze)), t, sp.oo) == sp.Rational(9, 2))
ok("no_peak_base", sp.simplify((Z * (3 - 2 * Z) / 3 - 1 / (3 * Z))
                                + (Z - 1) ** 2 * (2 * Z + 1) / (3 * Z)) == 0)
ok("no_peak_triple_family",
   sp.simplify(sol[w0] - sol[r0] + (rho - 1) ** 2 * (rho + 2) / (3 * rho**2)) == 0)
R["base"] = {"initial": ["25/48", "5/24", "1/240", "4/15"], "rho": "4/5",
             "tail_constant": "sqrt(3/8)", "limit_t_times_one_minus_Reff": "9/2"}

# ---------------------------------------------------------------- double zero
hh, rr = sp.symbols("h r_star", positive=True)
a0d = hh / rr**2
w0d = (3 * hh + 1) / (2 * rr) - 2 * hh / rr**2
r0d = rr * (1 - hh) / 2
Qd = a0d * q**3 - (w0d + 2 * a0d) * q**2 + q - r0d
etad = rr * (1 - hh) / (2 * hh)
ok("double_factorisation", sp.simplify(Qd - a0d * (q - rr) ** 2 * (q - etad)) == 0)
J2 = sp.Matrix([a0d, w0d, r0d]).jacobian([hh, rr])
minor = sp.factor(J2[[0, 2], :].det())
ok("double_family_rank2", sp.simplify(minor) != 0)
R["double_family_jacobian_minor"] = str(minor)
J3 = sp.Matrix([sol[s0], sol[w0], sol[r0]]).diff(rho)
ok("triple_family_rank1", sp.simplify(J3[0]) != 0)
ex = {hh: sp.Rational(1, 2), rr: sp.Rational(9, 10)}
d_ex = [sp.nsimplify(e.subs(ex)) for e in (a0d, w0d, r0d)]
ok("double_example_data", d_ex == [sp.Rational(50, 81), sp.Rational(25, 162), sp.Rational(9, 40)]
   and 1 - sum(d_ex) == sp.Rational(11, 3240))
ok("double_example_poly", sp.expand(Qd.subs(ex) - (10 * q - 9) ** 2 * (20 * q - 9) / 3240) == 0)
R["double_example"] = {"initial": ["50/81", "25/162", "11/3240", "9/40"],
                       "endpoint": ["1/2", "1/4", "0", "1/4"],
                       "endpoint_after_jump": ["1/8", "3/8", "0", "1/2"],
                       "tail": "4 t^-2"}

# ------------------------------------------------- general-rate chain endpoint
al, be, ga = sp.symbols("alpha beta gamma", positive=True)
A_, B_, C_ = sp.symbols("A B C", positive=True)
Xu = be * B_ - ga * C_
Xuu = al * be * A_ - be**2 * B_ - ga**2 * C_
Xuuu = -al * be * (al + be) * A_ + be**3 * B_ - ga**3 * C_
tri = sp.solve([Xu, Xuu, A_ + B_ + C_ - 1], [A_, B_, C_], dict=True)[0]
ok("triple_endpoint_general", sp.simplify(Xuuu.subs(tri) + al * be * ga) == 0)

# ---------------------------------------------------------------- hierarchy
d = sp.symbols("d", positive=True)
hier = {}
for m in range(2, 22, 2):
    rr_ = sp.Rational(1, m + 1)
    z = 1 + d
    w1 = rr_ * z * sum((1 + l) * (-d) ** l for l in range(m))
    hier[str(m)] = sp.simplify(w1 - rr_ / z + rr_ * d**m * (m * d + m + 1) / z) == 0
ok("no_peak_hierarchy_m2_to_20", all(hier.values()))
# odd m fails: active fraction of the formula is negative
m = 3
rr_ = sp.Rational(1, m + 1)
xodd = sp.nsimplify(rr_ * (-(1 - (1 + sp.Rational(1, 2 * m)))) ** (m + 1) / (1 + sp.Rational(1, 2 * m)))
z0 = 1 + sp.Rational(1, 2 * m)
x_m3 = 1 - sum(rr_ * (1 - (1 - z0) ** j) for j in range(1, m + 1)) - rr_ / z0
ok("odd_m3_negative_active_fraction", x_m3 < 0)
R["odd_m3_active_fraction_at_t0"] = str(x_m3)

# ---------------------------------------------------------------- SWIR
nu, et = sp.Integer(4), sp.Integer(4)
us = sp.Rational(1, 10)
Ss, Ws = 1 / nu, 1 / et
Rs = 1 - Ss - Ws
S0 = Ss * sp.exp(nu * us)
W0 = Ws * sp.exp(et * us) - nu * S0 * us
R0 = Rs - us
I0 = sp.simplify(1 - S0 - W0 - R0)
ok("swir_I0_closed", sp.simplify(I0 - (3 - 2 * sp.exp(sp.Rational(2, 5))) / 5) == 0)
ok("swir_W0_closed", sp.simplify(W0 - 3 * sp.exp(sp.Rational(2, 5)) / 20) == 0)
Su = S0 * sp.exp(-nu * u)
Wu = sp.exp(-et * u) * (W0 + nu * S0 * u)
Ru = R0 + u
Iu = 1 - Su - Wu - Ru
ok("swir_exposure_eqs", sp.simplify(sp.diff(Wu, u) - (nu * Su - et * Wu)) == 0
   and sp.simplify(sp.diff(Su, u) + nu * Su) == 0)
ders = [sp.simplify(sp.diff(Iu, u, k).subs(u, us)) for k in range(4)]
ok("swir_triple_zero", ders[0] == 0 and ders[1] == 0 and ders[2] == 0
   and sp.simplify(ders[3] + nu * et) == 0)
ok("swir_positive_data", all(float(e) > 0 for e in (S0, W0, I0, R0)))
If = sp.lambdify(u, Iu, "mpmath")
pred = mp.sqrt(3) / 8
swir_tail = []
for k in (3, 4, 5):
    delta = mp.mpf(10) ** (-k)
    ue = mp.mpf(1) / 10 - delta
    tt = mp.quad(lambda x: 1 / If(x), [0, ue / 2, ue - 10 * delta, ue])
    swir_tail.append({"u_star_minus_u": f"1e-{k}", "t": fl(tt, 8),
                      "t^(3/2)*I/(sqrt(3)/8)": fl(tt**1.5 * If(ue) / pred, 8)})
ok("swir_tail_converges", abs(float(swir_tail[-1]["t^(3/2)*I/(sqrt(3)/8)"]) - 1) < 5e-3)
R["swir"] = {"rates": {"S->W partner I": 4, "W->I partner I": 4,
                       "I->R partners S,W,I,R": 1},
             "initial": ["exp(2/5)/4", "3exp(2/5)/20", "(3-2exp(2/5))/5", "2/5"],
             "initial_numeric": [fl(e, 10) for e in (S0, W0, I0, R0)],
             "endpoint": ["1/4", "1/4", "0", "1/2"], "tail": "sqrt(3)/8 t^-3/2",
             "tail_check": swir_tail}
dep = {"S": {"I"}, "W": {"I"}, "I": {"S", "W", "I", "R"}, "R": set()}
persistent = {"S", "W", "R"}


def simple_cycles(g):
    found = set()
    def dfs(start, node, path):
        for n in sorted(g[node]):
            if n == start:
                found.add(tuple(path))
            elif n not in path and n > start:
                dfs(start, n, path + [n])
    for s_ in sorted(g):
        dfs(s_, s_, [s_])
    return sorted(found)


cycles = simple_cycles(dep)
ok("swir_G2", all(any(dep[v_] & persistent for v_ in c) for c in cycles))
R["swir"]["dependency_cycles"] = ["->".join(c) for c in cycles]

# ---------------------------------------------------------------- support and rate criteria
def block_report(name, c, P, V):
    rep = {}
    for i in V:
        auto = {k: c.get((k, i, i), 0) for k in P}
        hik = {k: sum(val for (a, b, kk), val in c.items() if a == i and kk == k) for k in P}
        rep[i] = {"autocatalytic_inflow_from_persistent": {k: str(vv) for k, vv in auto.items() if vv},
                  "dominance_holds": bool(all(auto[k] <= hik[k] for k in P)
                                          and any(auto[k] < hik[k] for k in P))}
    R.setdefault("criterion_examples", {})[name] = rep
    return rep


cb = {("S", "W", "X"): 2, ("W", "X", "X"): 1, ("X", "R", "R"): 1}
ok("base_violates_dominance", block_report("base_chain", cb, ["S", "W", "R"], ["X"])["X"]["dominance_holds"] is False)
csw = {("S", "W", "I"): 4, ("W", "I", "I"): 4, ("I", "R", "S"): 1, ("I", "R", "W"): 1,
       ("I", "R", "I"): 1, ("I", "R", "R"): 1}
ok("swir_violates_dominance", block_report("swir_eta4", csw, ["S", "W", "R"], ["I"])["I"]["dominance_holds"] is False)
csw1 = dict(csw)
csw1[("W", "I", "I")] = sp.Rational(1, 2)
ok("swir_eta_half_dominance", block_report("swir_eta_half", csw1, ["S", "W", "R"], ["I"])["I"]["dominance_holds"])

# ---------------------------------------------------------------- slowing law
def slow_check(name, F, G, U, k, ak, epsilons, branch):
    out = []
    g = G(U)
    for e in epsilons:
        e = mp.mpf(e)
        f = lambda x: F(x) + e * G(x)
        scale = abs(e * g / ak) ** (mp.mpf(1) / k)
        root = mp.findroot(f, U + branch(e) * scale)
        rate = -mp.diff(f, root)
        law = k * ak ** (mp.mpf(1) / k) * abs(e * g) ** (mp.mpf(k - 1) / k)
        out.append({"epsilon": mp.nstr(e, 3), "rate": fl(rate, 10),
                    "rate_over_law": fl(rate / law, 10),
                    "shift_over_scale": fl((root - U) / scale, 10)})
    R.setdefault("slowing_law", {})[name] = {"order_k": k, "a_k": fl(ak, 12),
                                             "g": fl(g, 12), "results": out}
    return out


a0b = mp.mpf(25) / 48
rb_ = mp.mpf(4) / 5
Fb = lambda x: a0b * (mp.e ** (-x) - rb_) ** 3 / mp.e ** (-x)
Gb = lambda x: mp.e ** x
Ub = mp.log(mp.mpf(5) / 4)
eps_b = ["1e-6", "-1e-6", "1e-12", "-1e-12"]
o1 = slow_check("base_triple", Fb, Gb, Ub, 3, mp.mpf(1) / 3, eps_b, lambda e: 1 if e > 0 else -1)
ok("slowing_triple", all(abs(float(r["rate_over_law"]) - 1) < 1e-3 for r in o1[2:]))
for r_, e in zip(o1, eps_b):
    exact = 3 * a0b ** (mp.mpf(1) / 3) * abs(mp.mpf(e)) ** (mp.mpf(2) / 3)
    ok("slowing_triple_matches_exact_" + e, abs(mp.mpf(r_["rate"]) / exact - 1) < 1e-8)
a0x = mp.mpf(50) / 81
rx = mp.mpf(9) / 10
etx = mp.mpf(9) / 20
Fd = lambda x: a0x * (mp.e ** (-x) - rx) ** 2 * (mp.e ** (-x) - etx) / mp.e ** (-x)
Ud = mp.log(mp.mpf(10) / 9)
o2 = slow_check("double_example", Fd, Gb, Ud, 2, a0x * rx * (rx - etx), ["-1e-8", "-1e-12"], lambda e: -1)
ok("slowing_double", all(abs(float(r["rate_over_law"]) - 1) < 1e-3 for r in o2))
fpos = lambda x: Fd(x) + mp.mpf("1e-8") * Gb(x)
ok("double_no_root_near_U_for_positive_eps",
   all(fpos(Ud + dlt) > 0 for dlt in (mp.mpf("-1e-3"), 0, mp.mpf("1e-3"))))
jroot = mp.findroot(fpos, mp.log(mp.mpf(20) / 9) - mp.mpf("1e-3"))
R["slowing_law"]["double_example"]["positive_eps_1e-8_first_root_u"] = fl(jroot, 10)
R["slowing_law"]["double_example"]["positive_eps_1e-8_final_r"] = fl((mp.mpf(9) / 40 - mp.mpf("1e-8")) * mp.e ** jroot, 10)
m = 4
rm = mp.mpf(1) / (m + 1)
Um = mp.log(1 + mp.mpf(1) / (2 * m))
Fh = lambda x: rm * (mp.e ** (Um - x) - 1) ** (m + 1) * mp.e ** (-(Um - x))
o3 = slow_check("hierarchy_m4", Fh, Gb, Um, m + 1, rm, ["1e-10", "-1e-10", "1e-20", "-1e-20"],
                lambda e: 1 if e > 0 else -1)
ok("slowing_order5", all(abs(float(r["rate_over_law"]) - 1) < 1e-3 for r in o3[2:]))

OUT.write_text(json.dumps(R, indent=2, sort_keys=True) + "\n")
print("PASS", sum(R["checks"].values()), "of", len(R["checks"]), "checks")
