---
type: paper
domain: apf
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-17
sources:
  - Paper_13_v8_1 (1).tex
  - Paper_13_v8_3-PLEC.pdf
  - Paper_13_v8_2.pdf
  - Edit Proposals - Paper 13 Descriptive Framing.md
---

# Paper 13 - Minimal Admissibility Core
## "The Master Reference"

**Status:** v8.3-PLEC (Phase 2.8a PLEC propagation, April 2026)

## Changelog

- **v8.3-PLEC (2026-04-17):** Phase 2.8a surgical PLEC propagation. Header/fancyhead/title retagged v8.3-PLEC; epigraph (already the canonical PLEC lay statement) annotated with a framing note identifying PLEC and its four structurally necessary components (A1, MD, A2, BW); abstract expanded with a PLEC-alignment sentence pointing to where the non-A1 components enter (L_col argmin, C_total=61 count); L_col subsection gains a "note on PLEC alignment" reconciling Paper 13's "derived from A1" shorthand with Paper 1 v4.0-PLEC's treatment of A2 as structurally necessary; v8.3-PLEC row added to Appendix J version history. No theorems, numerics, or code references changed. Full v9.0 restructure around PLEC's four-component inventory deferred. Compiled clean via xelatex (three passes), 64pp.
- **v8.2 (2026-04-16):** Descriptive-framing pass applied. 18 edits (P13.0–P13.17) from [[Edit Proposals - Paper 13 Descriptive Framing]]: new abstract reader's note, new `\subsection{Descriptive versus selective readings}` framing box (§2.4), section title rewrites ("What the Axiom Forces" → "What the Axiom Admits"; "Minimality — nature is parsimonious" → "Minimality — the admissibility set is an $\arg\min$"), and fourteen running-prose rewrites converting "selects/forces/selection" to "admits/unique admissible/admissibility". Incidental fix to pre-existing `\not\xrightarrow` syntax in §A.3 (rewritten as prose). Zero formal content changes. Compiled via xelatex (paper uses fontspec); two passes clean; 63pp.
- **v8.1 (pre-April 2026):** Comprehensive self-contained exposition. 18 modules, 335 bank-registered theorems, 348 verify_all checks (counts canonicalized to codebase v6.8 on 2026-04-18; v8.1 was originally quoted as "294+ theorems, 349 checks" which turned out to be stale metadata against the actual loaded count).

**Core claim:** A comprehensive, self-contained exposition of the entire APF framework in minimal form. All 18 bank-registered modules + apf/standalone/, 335 bank-registered theorems, 348 verify_all checks (codebase v6.8, 2026-04-18), unified derivation from [[Axiom A1]]. The canonical reference for understanding the full architecture.

**Structure:**
- Layer 1 (Spine): Axiom, operability, canonical object
- Layer 2 (Structure): Non-closure, gauge uniqueness
- Layer 3 (Ledgers): Cost ledgers, irreversibility (abbreviated)
- Layer 4 (Constraints): Field content, fermion spectrum, masses
- Layer 5 (Quantum): Hilbert space, Born rule, CPTP evolution
- Layer 6 (Dynamics): Spacetime, gravity, cosmology
- Layer 7 (Action): Internalization, Lagrangian, BSM
- Cross-cutting: Predictions, open problems, reconstruction pathways

**Content:**
- Unified notation and conventions
- All core theorems with proof sketches
- Module-by-module codebase map
- 47+ predictions with experimental status
- Open problems and future directions
- Reading guide for papers 1-7

**Implementation:** All 18 modules reference v8.1

**Status notes:**
- Comprehensive and stable
- Acts as a master index
- Suitable as a standalone reference

## See also
- [[Derivation Chain]]
- [[Papers 1-7|Paper 1 - Spine]]
- [[APF Codebase]]
