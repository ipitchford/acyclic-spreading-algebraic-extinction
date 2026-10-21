# Assurance

| Dimension | State | Note |
|---|---|---|
| Written proofs | present | Complete written proofs of every theorem; not machine-checked. The general statements do not depend on any computation |
| Exact checks | passed | `code/check_chain_family.py` (exact rationals, hierarchy m = 2, 4, …, 20) and the symbolic parts of `code/check_extensions.py` |
| High-precision / numerical checks | passed | 60-digit checks in `code/check_extensions.py`; `code/check_four_state.py`; `code/numerical_checks_and_figures.py` (integration of the original equations) |
| Exploratory simulation | labelled | `code/finite_population_gillespie.py`, fixed seeds; no theorem |
| Producer replay | passed | `verify_all.py` reruns every script in a temporary copy and compares all 12 outputs byte for byte. Clean full replays of the final package were run in two separate copies before release; the released scripts and outputs were compared with them by SHA-256 |
| Negative controls | passed | `negative_controls.py` (standard library, exact): corrupted rates, initial data, trajectory, graphs (added cycle or self-loop) and tail exponent/constant are all rejected |
| Independent rerun | not assessed | No unaffiliated rerun |
| Formal verification | not assessed | No proof assistant |
| Specialist review | not assessed | AI-assisted internal and adversarial review is not external specialist review |
| Editorial peer review | not assessed | None |
| Prior art | bounded | Earlier Evidence Press release credited as first refutation; dated bounded search in [RESEARCH_GATES.json](RESEARCH_GATES.json); no priority claim |

Trust boundary: the general statements rest on the written proofs. Exact checks cover the stated finite instances; floating-point outputs are reproducible with the pinned library versions. The simulation is exploratory. The model class is deterministic and no data are analysed; the critical initial data are tuned.
