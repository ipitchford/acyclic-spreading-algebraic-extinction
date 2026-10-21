# Internal review summary (public, rating-free)

This release was checked by internal AI-assisted review in several roles (editor, methodology, domain, perspective and devil's-advocate), an AI-assisted adversarial review by a different model (OpenAI gpt-5.6-sol) with a follow-up audit, and a verification re-review. None of this is external peer review, independent reproduction or formal verification. Original reports are kept privately; this summary keeps the substantive objections and responses.

No reviewer mathematical claim was found to be wrong. Every item below was re-verified by the package scripts or by hand.

## Main objections and responses

1. **Exact hypotheses.** Reviewers asked for verbatim quotations of the source's hypotheses and an item-by-item check. Response: §2.2 quotes the source (arXiv:2511.02437v1) with page pinpoints; Table 3 checks the counterexample against each hypothesis; the dependency-edge orientation and the persistence reading were confirmed against the source's Definitions 7–8; the counterexamples do not depend on the orientation convention. The source's Corollary 1 is quoted with its identifiability condition.
2. **SWIR / generalised epidemic process.** The domain review pointed out that a mean-field SWIR-type model at a tricritical-type point already violates the conjecture under the source's recovery encoding. Response: this is stated plainly in §1 and §3.4 and credited to Janssen–Müller–Stenull (2004) and Choi–Lee–Kahng; the reviewer's construction was re-verified exactly (triple zero, tail (√3/8) t^(−3/2)). After the adversarial audit the wording was narrowed to "a no-direct-infection SWIR/GEP submodel admits a tricritical-type triple zero", explicitly not identified with the point analysed by Janssen et al. The novelty was recentred on what is new relative to SWIR.
3. **Earlier Evidence Press release.** Before publication it was found that an earlier release from the same publisher (doi:10.5281/zenodo.23501636) had already refuted the conjecture (18/t², cyclic dependency graph), proved four compartments necessary, built a family with exponents m/(m−1) and proved the limiting-Jacobian log-rate theorem. Response: the paper is framed as a successor; the abstract, §1, Table 1, Table 7 and §9 credit it; minimality and Theorem 7.1 are marked as the same statements; no priority wording remains. The cited details were checked against the archived PDF.
4. **Honest framing of tuning.** The critical data are tuned (codimension two). Response: stated after Theorem 3.1; the t^(−2) example (codimension one) has its own section; genericity is proved in Theorem 6.2.
5. **Overstated repair.** "Exact spectral replacement" was renamed; Theorem 7.1 is described as an a posteriori criterion using standard Metzler/Perron–Frobenius tools, silent at s(J) = 0. "Graph-checkable" became "support-checkable", because the labels are read off the graphs only under the source's identifiability condition. The rate-only negative result was restricted to criteria that ignore the initial data.
6. **Proof details.** Unstated steps were written out (Lemma 2.2 successor clause, Proposition 5.3 hypothesis, Theorem 7.1 diagonal placement, the double-family argument, the minimality clause, triple-case attainability, the Theorem 5.7 case {2,3} estimate, the Theorem 6.2(c) containment and stratification, the Theorem 6.3 endpoint expansion).
7. **Scope of no-peak.** Restricted throughout to the chain, SWIR and hierarchy families.
8. **Literature.** 32 references checked by DOI, arXiv or library record; incorrect attributions (SWIR to Chung et al.; titles of Anderson 2011 and Choi–Lee–Kahng) were corrected. The Pólya–Szegő problem number could not be verified, so the lemma is proved in the text.
9. **Readability.** The paper was rewritten reader-first: the problem in words, the R_eff picture, a contributions table with status labels, a notation table, a diagram of both graphs and a derivation of the counterexample instead of a verification by differentiation.

## Not done / open

- Siphon restatement of the graph conditions (open problem).
- No contact with the source authors before release.
- No proof-assistant formalisation.
- The Leiden dissertation of the source's first author could not be accessed.
