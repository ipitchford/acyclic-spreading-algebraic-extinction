# Agent navigation

## Identity

- Release slug: `acyclic-spreading-algebraic-extinction`; version `1.0.0-candidate`; DOI 10.5281/zenodo.23548109 (concept DOI 10.5281/zenodo.23548108)
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md))
- **Successor release.** Parent: *Algebraic extinction in acyclic compartment models*, doi:10.5281/zenodo.23501636 (https://evidencepress.org/releases/algebraic-extinction-compartments/), which first refuted the same conjecture. Overlap: paper §9 "Relation to the prior release".
- Manuscript: [PDF](paper/spreading_decay.pdf), [LaTeX source](paper/spreading_decay.tex); human summary and status table: [README.md](README.md); file list: [MANIFEST](MANIFEST)
- Machine-readable claims: [CLAIMS.json](CLAIMS.json); frozen target with claims, exclusions and source pinpoints: [evidence/target.yaml](evidence/target.yaml)
- Assurance limits: [ASSURANCE.md](ASSURANCE.md); internal review summary: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)
- Provenance: [PROVENANCE.md](PROVENANCE.md); prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json); source versions and hashes: [sources/provenance.json](sources/provenance.json)
- Citation: [CITATION.cff](CITATION.cff), [Zenodo metadata](.zenodo.json); licences: [LICENSES.md](LICENSES.md)

## Scientific map

Object: conservative quadratic spreading models y_j' = Σ β^j_{ik} y_i y_k − … of Chan, Crommelin and Mandjes (arXiv:2511.02437v1), with transition graph and dependency graph; Conjecture 1 (p. 25): under (G1) acyclic transitions and (G2) the dependency-cycle condition, every vanishing compartment decays exponentially.

- **Base counterexample** (Theorem 3.1, proof; exact): s' = −2sx, w' = 2sx − wx, x' = (w − r)x, r' = rx; data (25/48, 5/24, 1/240, 4/15); ϑ = (16 + 2t/3)^(−1/2), x = ϑ³/(3(1+ϑ)) ~ √(3/8) t^(−3/2); both graphs acyclic; R_eff = w/r → 1, 1 − R_eff ~ 9/(2t).
- **SWIR-type submodel** (§3.4, proof; exact): no-direct-infection SWIR/GEP submodel (Janssen–Müller–Stenull 2004; Choi–Lee–Kahng) with recovery as in source eq. (8) is in the class and has a tricritical-type triple zero, tail (√3/8) t^(−3/2). Not identified with the tricritical point analysed by Janssen et al.
- **Bounds** (Theorems 4.1, 5.1, 5.7): N ≤ 3 impossible (as in the parent); tails t^(−1−1/m) for even m; any N = 4 model: exponential, t^(−2) or t^(−3/2).
- **Structure** (Theorems 6.1–6.3, Proposition 6.5): chain classification; algebraic extinction on a null set of data for y' = xHy (codimensions 1 and 2 in the chain); slowing law rate ∝ |ε|^((k−1)/k); no outbreak peak precedes a slow tail in the chain, SWIR and hierarchy families.
- **Repairs** (Theorem 7.1, same as parent; Corollaries 7.3–7.4): log-rate = s(J) of the limiting Metzler block; support-checkable and rate-checkable conditions for exponential decay.
- **Finite populations** (§8, exploratory): median extinction ≈ 1.48 N^0.42.

Do not read: the conjecture's refutation as new (the parent refuted it first); the critical data as generic; the simulation as a theorem; any statement as an epidemic forecast.

## Replay

Python 3.13 with numpy, scipy, sympy, mpmath, matplotlib ([requirements.txt](requirements.txt); [ENVIRONMENT.txt](ENVIRONMENT.txt)):

```sh
python3 -B -O negative_controls.py         # exact controls; corruptions must be rejected (standard library)
python3 -B verify_all.py                   # rerun every script in a temporary copy; byte-for-byte comparison
python3 -B verify_all.py --skip-simulation # without the stochastic simulation
shasum -a 256 -c SHA256SUMS
```

Scripts: [code/check_chain_family.py](code/check_chain_family.py), [code/check_extensions.py](code/check_extensions.py), [code/check_four_state.py](code/check_four_state.py), [code/numerical_checks_and_figures.py](code/numerical_checks_and_figures.py), [code/finite_population_gillespie.py](code/finite_population_gillespie.py); outputs in `code/*.json`, `evidence/`, `figures/`.

## Useful next research question

A stochastic theory of the near-tricritical endgame in finite populations, or the slowest tail and genericity beyond the single-catalyst class (paper §9, open problems 1–4).
