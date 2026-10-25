# Corrections

Errata to release 1.0.0-candidate (DOI 10.5281/zenodo.23548109). The archived PDF,
`paper/spreading_decay.tex` on this branch, `SHA256SUMS`, the Zenodo deposit and the tag
`v1.0.0-candidate` are unchanged.

## 2026-10-25 — Table 4, time scale for m = 6

Table 4 (`paper/spreading_decay.tex` line 745) gives $T_6\approx3.48\times10^{7}$.
With $T_m=(m+1)(2m)^m/m$, $T_6=7\cdot12^6/6=3{,}483{,}648$, so the correct entry is
$T_6\approx3.48\times10^{6}$.

## 2026-10-25 — §3.2, dimension of the simplex of initial data

§3.2 ("How tuned is this?", line 540) says the $t^{-3/2}$ trajectories form a curve "in the
two-dimensional simplex of initial data". It is the **three-dimensional** simplex of initial
data; the stated codimension two of the curve is correct for that simplex.

Neither statement is used in any proof; the theorem, certificates and evidence state are
unaffected. Also recorded at <https://evidencepress.org/releases/acyclic-spreading-algebraic-extinction/>.
