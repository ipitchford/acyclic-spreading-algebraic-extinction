# Algebraic extinction in acyclic spreading models

Explicit counterexamples, genericity and support-checkable repairs.
Author: Anonymous. 21 October 2026. Version 1.0.0-candidate. DOI: https://doi.org/10.5281/zenodo.23548109.

## In plain English

Compartment models of epidemics and rumours divide a population into
groups, such as "unaware", "heard once", "spreading" and "stopped". People
move between groups when they meet someone. In models like SIR the number
of active spreaders eventually falls by a fixed fraction per unit time,
which is exponential fade-out. Chan, Crommelin and Mandjes
(arXiv:2511.02437v1) conjectured that this always happens when two graphs
describing the model are acyclic or nearly so.

The conjecture was already refuted in an earlier release from the same
publisher, Evidence Press:
*Algebraic extinction in acyclic compartment models* (Anonymous, 20 October
2026, https://doi.org/10.5281/zenodo.23501636). That release gave a
four-group example that fades like `18/t^2`, proved that four groups are
needed, and gave a spectral repair. This paper is a successor to it.

The reason the conjecture fails is that an outbreak can end exactly
at the threshold where each spreader replaces itself (effective
reproduction number 1). The activity then fades like a power of time.
This paper adds the following:

- a four-group model in which *both* graphs are acyclic and the
  trajectory is explicit, fading like `t^(-3/2)`;
- a proof that no four-group model fades more slowly than `t^(-3/2)`;
- that power-law fade-out needs finely tuned starting conditions in models
  driven by one active group;
- that exponential decay returns at a predictable rate after a perturbation;
- that in the families studied the slow tails never follow an outbreak peak;
- simple conditions on the coefficient supports or the rates that rule the
  mechanism out (the support conditions need the source's identifiability
  condition to be read off the graphs);
- a simulation of what happens in finite populations.

A SWIR-type two-step contagion submodel also gives a `t^(-3/2)`
counterexample. The paper's §9 says exactly which results overlap with
the earlier release.

## What is proved, computed, simulated or open

| Result | Where in the paper | Status |
|---|---|---|
| A no-direct-infection SWIR/GEP submodel admits a tricritical-type triple zero, hence a `t^(-3/2)` counterexample | §3.4 | proof; exact check |
| Four-state model with both graphs acyclic, explicit trajectory, `x ~ sqrt(3/8) t^(-3/2)` | Theorem 3.1 | proof; exact check |
| `t^(-2)` example of codimension one | §3.3, eq. (9) | proof; exact check |
| Three compartments never suffice (same statement as the prior release, Thm 5) | Theorem 4.1 | proof |
| Tails `t^(-1-1/m)` for every even `m`, both graphs acyclic; slower at each even N than the prior family `t^(-m/(m-1))` | Theorem 5.1 | proof; exact check `m=2..20` |
| Sharp zero-order bound for `y' = x H y` only | Proposition 5.3 | proof |
| Slowest tail at `N=4` with one vanishing state is `t^(-3/2)` | Corollary 5.4 | proof |
| Any four-compartment model: exponential, `t^(-2)` or `t^(-3/2)`; nothing slower | Theorem 5.7 | proof; numerical examples (Remark 5.8) |
| Complete classification of the four-state chain | Theorem 6.1 | proof |
| Algebraic extinction confined to a null set (class `y'=xHy`); codimensions 1 and 2 in the chain | Theorem 6.2 | proof; exact Jacobians |
| General-order slowing law `rate ∝ |ε|^((k-1)/k)` | Theorem 6.3 | proof; high-precision check |
| No outbreak peak before a slow tail (chain, SWIR and hierarchy families only) | Proposition 6.5 | proof; exact identities |
| Log-rate equals the spectral abscissa of the limiting Metzler block (same result as the prior release, Thm 7 and Cor. 8) | Theorem 7.1 | proof (needs the limit) |
| Support-checkable and rate-checkable conditions for exponential decay | Corollaries 7.3–7.4 | proof |
| Finite populations: extinction at about `1.5 N^0.42` | §8 | exploratory simulation |

"Proof" means a complete written argument in the paper. No proof
assistant has checked these proofs. "Exact" means a symbolic or
exact-rational computation. "Numerical" means a floating-point check, and
"exploratory" means a simulation without proof.

## Files

- `paper/spreading_decay.pdf`, `paper/spreading_decay.tex`: the manuscript.
- `code/check_chain_family.py` → `code/chain_verification.json`: exact
  rational check of the hierarchy, `m = 2, 4, …, 20` (standard library only;
  internal names `q_k, R` are the paper's `w_k, r`).
- `code/check_extensions.py` → `code/extensions_verification.json`:
  symbolic and 60-digit checks of the derivation, the `R_eff` identity,
  the double-zero family and its Jacobian, the SWIR example, the
  monotonicity identities, the dominance examples and the slowing law.
- `code/check_four_state.py` → `code/four_state_verification.json`: checks
  for Theorem 5.7 and Remark 5.8 (about one minute).
- `code/numerical_checks_and_figures.py` → `evidence/*.csv`,
  `evidence/numerical_verification.json`, `figures/*`: integration of the
  original equations and the figures.
- `code/finite_population_gillespie.py` → `code/finite_population_gillespie.json`:
  the exploratory simulation, with fixed seeds.
- `evidence/target.yaml`: claims, status labels and source pinpoints.
- `sources/provenance.json`: URLs and SHA-256 hashes of the source versions
  inspected. Third-party papers are not redistributed.
- `MANIFEST`, `SHA256SUMS`: the file list and hashes.

## Replay

Requires Python 3.13 with `numpy`, `scipy`, `sympy`, `mpmath`, `matplotlib`
(see `requirements.txt`).

```bash
python verify_all.py                    # everything, about 90 seconds
python verify_all.py --skip-simulation  # without the stochastic simulation
```

`verify_all.py` copies the package to a temporary directory and reruns
every script there. It then compares every regenerated file with the
archived copy byte for byte; the package itself is not modified. The outputs
contain no timings or version stamps, and floating-point values are
rounded. The exact outputs are platform-independent. Floating-point
outputs and figures are byte-identical with the pinned library versions;
on other platforms the last digits may differ.

Rebuild the paper from `paper/`, with the figures one level up in
`figures/`:

```bash
pdflatex spreading_decay.tex && pdflatex spreading_decay.tex
```

## Source of the problem

Chan, K. M. D., Crommelin, D. T., & Mandjes, M. R. H. (2025). *Asymptotic
behavior for a general class of spreading models*. arXiv:2511.02437v1.
https://arxiv.org/abs/2511.02437v1. The paper quotes the hypotheses
verbatim with page pinpoints. Only version 1 exists on arXiv. A 2026
Leiden dissertation by the first author could not be accessed.
