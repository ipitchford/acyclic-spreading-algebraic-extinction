#!/usr/bin/env python3
"""Semantic negative controls for the base counterexample (Theorem 3.1).

Standard library only (fractions). The baseline checks must pass and every
corrupted variant must be rejected; the script exits non-zero otherwise.
Assertions are not used, so the controls also run under ``python -O``.

Baseline, with theta(t) = (16 + 2t/3)^(-1/2), zeta = 1 + theta:
  s = zeta^2/3, w = zeta(3 - 2 zeta)/3, x = theta^3/(3 zeta), r = 1/(3 zeta)
solves s' = -2sx, w' = 2sx - wx, x' = (w - r)x, r' = rx
with initial data (25/48, 5/24, 1/240, 4/15). Using theta' = -theta^3/3 the
ODE residuals are rational functions of theta and are checked exactly at
rational theta values. Graph acyclicity, the R_eff identity and the tail
constant t^(3/2) x -> sqrt(3/8) are checked too.
"""
from fractions import Fraction as F
import sys

RATES = {"SW": F(2), "WX": F(1), "XR": F(1)}
INIT = (F(25, 48), F(5, 24), F(1, 240), F(4, 15))


def traj(th):
    z = 1 + th
    return (z * z / 3, z * (3 - 2 * z) / 3, th ** 3 / (3 * z), 1 / (3 * z))


def dtraj(th):
    # derivatives with respect to t, using d theta/dt = -theta^3/3
    z = 1 + th
    dth = -th ** 3 / 3
    ds = 2 * z / 3 * dth
    dw = (3 - 4 * z) / 3 * dth
    dx = (3 * th ** 2 * z - th ** 3) / (3 * z * z) * dth
    dr = -1 / (3 * z * z) * dth
    return ds, dw, dx, dr


def residual(th, rates, traj_fn=traj):
    s, w, x, r = traj_fn(th)
    ds, dw, dx, dr = dtraj(th)
    a, b, c = rates["SW"], rates["WX"], rates["XR"]
    rhs = (-a * s * x, a * s * x - b * w * x, (b * w - c * r) * x, c * r * x)
    return [l - rr for l, rr in zip((ds, dw, dx, dr), rhs)]


def solves(rates=RATES, init=INIT, traj_fn=traj):
    th0 = F(1, 4)  # theta(0) = 16^(-1/2)
    if tuple(traj_fn(th0)) != tuple(init):
        return False
    for th in (F(1, 4), F(1, 5), F(1, 7), F(1, 10), F(1, 100), F(3, 1000)):
        if any(v != 0 for v in residual(th, rates, traj_fn)):
            return False
    return True


def acyclic(edges):
    nodes = {u for e in edges for u in e}
    indeg = {n: 0 for n in nodes}
    for _, v in edges:
        indeg[v] += 1
    stack = [n for n in nodes if indeg[n] == 0]
    seen = 0
    while stack:
        n = stack.pop(); seen += 1
        for u, v in edges:
            if u == n:
                indeg[v] -= 1
                if indeg[v] == 0:
                    stack.append(v)
    return seen == len(nodes) and all(u != v for u, v in edges)


TRANSITION = [("S", "W"), ("W", "X"), ("X", "R")]
DEPENDENCY = [("S", "X"), ("W", "X"), ("X", "R")]  # donor -> catalyst of its outflow


def reff_ok(th, claimed=lambda th: None):
    s, w, x, r = traj(th)
    reff = w / r  # production / removal rate of spreaders
    t = (1 / (th * th) - 16) * F(3, 2)
    return reff, t


def tail_constant_ok(power_num, power_den, target_sq):
    # t^(p) x(t) -> c with c^2 = target_sq, checked as a limit along theta -> 0
    vals = []
    for th in (F(1, 10 ** 3), F(1, 10 ** 4), F(1, 10 ** 5)):
        t = (1 / (th * th) - 16) * F(3, 2)
        x = traj(th)[2]
        # compare (t^p x)^(2 den) with target^(den) without irrational powers
        lhs = (t ** power_num) * (x ** (2 * power_den))
        vals.append(abs(lhs / (target_sq ** power_den) - 1))
    return vals[-1] < F(1, 10 ** 3) and vals[0] > vals[-1]


def main():
    ok = True
    def check(name, cond, expect):
        nonlocal ok
        good = (cond == expect)
        print(("PASS " if good else "FAIL ") + name + (" (rejected)" if not expect and good else ""))
        ok = ok and good

    # baseline
    check("baseline trajectory solves the ODE exactly with the stated data", solves(), True)
    check("baseline transition graph acyclic", acyclic(TRANSITION), True)
    check("baseline dependency graph acyclic", acyclic(DEPENDENCY), True)
    # (t^(3/2) x)^2 = t^3 x^2 -> 3/8
    check("baseline tail t^(3/2) x -> sqrt(3/8)", tail_constant_ok(3, 1, F(3, 8)), True)
    # 1 - R_eff ~ 9/(2t)
    th = F(1, 10 ** 4)
    reff, t = reff_ok(th)
    check("baseline 1 - R_eff ~ 9/(2t)", abs((1 - reff) * 2 * t / 9 - 1) < F(1, 50), True)

    # corruptions: each must be rejected
    check("corrupted rate S->W = 3", solves(dict(RATES, SW=F(3))), False)
    check("corrupted rate X->R = 2", solves(dict(RATES, XR=F(2))), False)
    check("corrupted initial x0 = 1/120", solves(init=(INIT[0], INIT[1], F(1, 120), INIT[3])), False)
    check("corrupted trajectory x = theta^2/(3 zeta)",
          solves(traj_fn=lambda th: (traj(th)[0], traj(th)[1], th ** 2 / (3 * (1 + th)), traj(th)[3])), False)
    check("dependency graph with an added cycle X->W", acyclic(DEPENDENCY + [("X", "W")]), False)
    check("transition graph with a self-loop", acyclic(TRANSITION + [("X", "X")]), False)
    # wrong exponent: t^2 x would need a finite nonzero limit; (t^2 x)^2 = t^4 x^2 diverges
    check("claimed t^(-2) tail with constant sqrt(3/8)", tail_constant_ok(4, 1, F(3, 8)), False)
    check("claimed t^(-3/2) tail with constant 1/2", tail_constant_ok(3, 1, F(1, 4)), False)

    print("ALL CONTROLS PASSED" if ok else "NEGATIVE CONTROLS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
