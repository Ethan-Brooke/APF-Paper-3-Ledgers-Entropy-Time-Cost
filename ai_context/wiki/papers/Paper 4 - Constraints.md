---
type: paper
domain: apf
layer: 4-constraints
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Paper 4 - Constraints
## "Gauge Group, Field Content, Masses, and Mixing"

**Status:** Updated PDF available (Paper_4_CONSTRAINTS_updated.pdf)

**Core claim:** The admissible gauge structure forces 45 fermions in 3 generations with exact mass ratios and mixing parameters. No mass hierarchy problem — the masses are determined by [[Enforcement Budget]] topology and [[Gram Matrix]] properties.

**Key theorems:**
- Theorem_R (R1/R2/R3 checks) — Representation forcing in SU(3)×SU(2)×U(1)
- L_count (`check_L_count`) — Exactly 61 effective degrees of freedom
- T_field (`check_field_content`) — 45 fermions forced; no exotic matter
- [[Gram Matrix]] properties — Determines mass/mixing structure
- T27c (`check_T27c`, x = 1/2) — Capacity-binding condition for three generations
- T_CKM (`check_T_CKM`) — Quark mixing angles
- T_PMNS (`check_T_PMNS`) — Neutrino mixing angles

**Content:**
- Gauge representation theory and generation counting
- Why exactly three families; fourth generation forbidden
- Fermionic mass ratios from Gram spectrum
- [[Non-Closure Theorem]] prevents additional particles
- Mixing matrix derivation: CKM and PMNS
- Comparison to experimental values

**Implementation:** `gauge.py` (29 checks) + `generations__8_.py` (86 checks)

**Status notes:**
- v6.7 complete and stable
- PDF updated; matches codebase

## See also
- [[Gram Matrix]]
- [[Predictions Catalog]]
- [[Paper 2 - Structure]]
