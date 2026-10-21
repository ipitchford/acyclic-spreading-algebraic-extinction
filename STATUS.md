# Status

Version: 1.0.0-candidate (21 October 2026). DOI: [10.5281/zenodo.23548109](https://doi.org/10.5281/zenodo.23548109). Scholarly creator: Anonymous.

**Unrefereed research candidate.** Not peer reviewed by a journal or by named human specialists. Checked by internal AI-assisted review in several roles, an AI-assisted adversarial review by a different model (OpenAI gpt-5.6-sol) and a verification re-review ([REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)). That is not peer review, independent reproduction or formal verification.

**Successor release.** Conjecture 1 of Chan, Crommelin and Mandjes (arXiv:2511.02437v1, p. 25) was first refuted by an earlier Evidence Press release from the same publisher: Anonymous (2026), *Algebraic extinction in acyclic compartment models*, version 1.0.0-candidate, doi:[10.5281/zenodo.23501636](https://doi.org/10.5281/zenodo.23501636). That release gave a four-compartment 18/t² example with a cyclic dependency graph, proved four compartments necessary, gave a family with exponents m/(m−1), and proved the limiting-Jacobian log-rate theorem. Section 9 of the paper states the overlap item by item. No priority is claimed.

* **Proved in writing:** Theorem 3.1 (four-state model, both graphs acyclic, explicit trajectory, x ~ √(3/8) t^(−3/2)); the t^(−2) example (§3.3); the no-direct-infection SWIR/GEP submodel counterexample with a tricritical-type triple zero (§3.4); Theorem 4.1 (three compartments never suffice; same statement as the prior release); Theorem 5.1 (tails t^(−1−1/m), every even m); Proposition 5.3 and Corollary 5.4; Theorem 5.7 (any four-compartment model: exponential, t^(−2) or t^(−3/2)); Theorem 6.1 (classification of the chain); Theorem 6.2 (genericity in the single-catalyst class); Theorem 6.3 (slowing law); Proposition 6.5 (no outbreak peak in the chain, SWIR and hierarchy families); Theorem 7.1 (same statement as the prior release); Corollaries 7.3–7.4 (support- and rate-checkable conditions).
* **Exact computation:** exact rational and symbolic checks (`code/check_chain_family.py`, `code/check_extensions.py`).
* **Numerical:** 60-digit and floating-point checks, integration of the original equations, Remark 5.8 examples.
* **Exploratory:** finite-population Gillespie simulation, median extinction ≈ 1.5 N^0.42 (§8).
* **Open:** slowest tail for N ≥ 5 outside the single-catalyst class; genericity beyond it; a peaked algebraic tail; a stochastic theory of the endgame; a siphon formulation (§9).

No proof assistant was used. The novelty statement is bounded by the searches summarised in [RESEARCH_GATES.json](RESEARCH_GATES.json).
