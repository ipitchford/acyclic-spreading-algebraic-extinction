#!/usr/bin/env python3
"""Independent exact-arithmetic check of the even-m acyclic chain family.

Uses only the Python standard library.  Coordinates are represented by
polynomial numerators N_i(d) over the common denominator z=1+d.  The checker
cross-multiplies the ORIGINAL reaction ODEs and compares polynomial
coefficients exactly.  It does not use numerical ODE integration.
"""

from fractions import Fraction as F
from math import comb
import json
from pathlib import Path


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for j, value in enumerate(a):
        out[j] += value
    for j, value in enumerate(b):
        out[j] += value
    return trim(out)


def scale(a, c):
    return trim([c * value for value in a])


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for j, aj in enumerate(a):
        for k, bk in enumerate(b):
            out[j + k] += aj * bk
    return trim(out)


def derivative(a):
    return trim([j * a[j] for j in range(1, len(a))] or [F(0)])


def evaluate(a, d):
    value = F(0)
    for coefficient in reversed(a):
        value = value * d + coefficient
    return value


def topological_check(edges, order):
    position = {node: j for j, node in enumerate(order)}
    return all(position[a] < position[b] for a, b in edges)


def verify(m):
    assert m >= 2 and m % 2 == 0
    r = F(1, m + 1)
    d0 = F(1, 2 * m)
    z = [F(1), F(1)]
    labels = [f"q{k}" for k in range(1, m + 1)] + ["x", "R"]
    x_index, r_index = m, m + 1
    numerators = []
    positivity_bounds = {}

    for k in range(1, m + 1):
        power = [F(comb(k + 1, j)) for j in range(k + 2)]
        alternating = [r * comb(k + ell, k) * (-1) ** ell
                       for ell in range(m - k + 1)]
        numerators.append(mul(power, alternating))
        if k < m:
            ratios = [F(k + ell + 1, ell + 1) * d0
                      for ell in range(m - k)]
            assert all(F(0) < ratio <= F(1, 2) for ratio in ratios)
            lower = 1 - (k + 1) * d0
            assert lower >= F(1, 2)
            positivity_bounds[f"q{k}"] = str(lower)
        else:
            positivity_bounds[f"q{k}"] = "1"

    numerators.append([F(0)] * (m + 1) + [r])
    numerators.append([r])

    # Independent mass identity, after multiplying by z.
    total = [F(0)]
    for numerator in numerators:
        total = add(total, numerator)
    assert total == z

    # Reconstruct donor, recipient, catalyst, coefficient from the reaction list.
    reactions = []
    for k in range(1, m + 1):
        donor = k - 1
        recipient = k - 2 if k > 1 else x_index
        reactions.append((donor, recipient, x_index, F(k)))
    reactions.append((x_index, r_index, r_index, F(1)))

    rhs = [[F(0)] for _ in labels]
    for donor, recipient, catalyst, coefficient in reactions:
        flux = scale(mul(numerators[donor], numerators[catalyst]), coefficient)
        rhs[donor] = add(rhs[donor], scale(flux, -1))
        rhs[recipient] = add(rhs[recipient], flux)

    # y_i=N_i/z and d'=-r d^(m+1).  Multiply y_i'=reaction RHS by z^2.
    d_prime = [F(0)] * (m + 1) + [-r]
    for numerator, target in zip(numerators, rhs):
        quotient_derivative_numerator = add(mul(z, derivative(numerator)),
                                            scale(numerator, -1))
        lhs = mul(quotient_derivative_numerator, d_prime)
        assert lhs == target

    initial = [evaluate(numerator, d0) / (1 + d0) for numerator in numerators]
    limit = [numerator[0] for numerator in numerators]
    assert all(F(0) < value < F(1) for value in initial)
    assert sum(initial) == 1
    assert limit == [r] * m + [F(0), r]
    assert sum(limit) == 1

    transition_edges = {(donor, recipient) for donor, recipient, _, _ in reactions}
    dependency_edges = {(donor, catalyst) for donor, _, catalyst, _ in reactions}
    transition_order = list(range(m - 1, -1, -1)) + [x_index, r_index]
    dependency_order = list(range(m)) + [x_index, r_index]
    assert topological_check(transition_edges, transition_order)
    assert topological_check(dependency_edges, dependency_order)
    assert dependency_edges == {(k, x_index) for k in range(m)} | {(x_index, r_index)}

    # x~C t^(-(m+1)/m), C^m=r^m/(m*r)^(m+1).
    coefficient_power = r ** m / (m * r) ** (m + 1)
    assert coefficient_power == F(m + 1, m ** (m + 1))
    assert numerators[x_index][m + 1] == r
    assert all(value == 0 for value in numerators[x_index][:m + 1])

    return {
        "m": m,
        "states": m + 2,
        "result": "PASS",
        "initial": dict(zip(labels, map(str, initial))),
        "limit": dict(zip(labels, map(str, limit))),
        "exponent": str(F(m + 1, m)),
        "leading_coefficient_to_power_m": str(coefficient_power),
        "integrated_exposure": f"log({1 + d0})",
        "alternating_sum_lower_bounds": positivity_bounds,
        "transition_order": [labels[j] for j in transition_order],
        "dependency_order": [labels[j] for j in dependency_order],
        "checked": ["positive rational initial fractions", "total mass polynomial identity",
                    "all original reaction ODEs as polynomial identities", "endpoint fractions",
                    "transition DAG", "dependency DAG", "alternating-sum positivity bound",
                    "zero multiplicity", "asymptotic coefficient power"]
    }


if __name__ == "__main__":
    output = {"method": "standard-library Fraction polynomial arithmetic",
              "cases": [verify(m) for m in range(2, 21, 2)]}
    first = output["cases"][0]
    assert first["initial"] == {"q1": "5/24", "q2": "25/48", "x": "1/240", "R": "4/15"}
    assert first["leading_coefficient_to_power_m"] == "3/8"
    destination = Path(__file__).with_name("chain_verification.json")
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(f"PASS: {len(output['cases'])} even-m cases (m=2,...,20); exact ODE and graph identities.")
    print(destination)
