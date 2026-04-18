---
type: paper
domain: apf
layer: 0-ontology
created: 2026-04-14
updated: 2026-04-17 (v3.0 rebuild)
sources: []
---

# What Physics Permits
## Paper 0 — An Ontology for Admissibility Physics

**Status:** **v3.0 (April 2026), 56pp .tex/.pdf**. Rebuilt from scratch to replace the Feb-2026 v2.0 .docx "gateway" narrative. v3.0 is the ontology paper of the APF series: the reference statement of what the framework commits to, from which all technical papers derive.

**Core claim:** The framework's ontological commitments are a four-layer hierarchy: (1) the keystone — *meaning requires enforceable distinction*; (2) finiteness of enforcement (A1); (3) admissibility structure (A1 + MD + BW); (4) realization as the $\arg\min$ over the admissible set (A2; the [[Principle of Least Enforcement Cost]] unifies all four). PLEC is stated canonically: *reality is the minimum-cost expression of distinction compatible with finite capacity*. This sentence is the framework's frozen lay handle.

**Three-book structure:**

**Book I — The Ontology** (Ch 1–5). Substantive core.
- Ch 1 The Ontological Keystone. Meaning requires enforceable distinction. Why the keystone is prior to dynamics and prior to representation. What the keystone rules out (zero-cost distinctions, distinctions with no consequence, notational fictions).
- Ch 2 The Finiteness of Enforcement (A1). A1 stated with its capacity-bound formula. What A1 does (admissible set from above). What A1 does not do (doesn't exclude zero-cost distinctions; doesn't select realization; doesn't impose cost variation; doesn't fix capacity value). Stanford-reviewer countermodel separating A1 from MD.
- Ch 3 The Admissibility Structure. A(ρ, Γ) with all four bounds. MD (positive floor) with independence countermodel. BW (non-degeneracy). Four-layer hierarchy stated canonically. Structural independence of all four components (pairwise logical independence with named countermodels).
- Ch 4 The Principle of Least Enforcement Cost (PLEC). Formal statement with canonical lay sentence. The argmin is a locator, not a process. Scope of the argmin (over A, not outside it). PLEC in the variational tradition (Fermat, Maupertuis, Hamilton, Jaynes). Non-merge discipline.
- Ch 5 Descriptive Reading. The convention stated. Why it matters. How the paper uses language (sentence-level rewrites, not a global flag). Applied glossary of common narrative↔descriptive translations.

**Book II — Derivational Commitments** (Ch 6–13). One chapter per paper under fixed six-section template (Scope / Commitments drawn on / Structural results / Where "forces" is a descriptor / What the paper does not claim / Status). Expands as papers mature without structural churn.
- Ch 6 Paper 1 (Enforceability of Distinction) — full.
- Ch 7 Paper 2 (Structure of Admissible Physics) — full.
- Ch 8 Paper 3 (Ledgers) — full.
- Ch 9 Paper 4 (Constraints — Field Content) — **placeholder** (pending .tex rewrite; current codebase content cited).
- Ch 10 Paper 5 (Quantum Structure) — **placeholder** (remaining scope under review post-Paper 1 v4.0 absorption).
- Ch 11 Paper 6 (Dynamics and Geometry) — **placeholder** (pending .tex rewrite).
- Ch 12 Paper 7 (Action, Internalization, and the Lagrangian) — full.
- Ch 13 Paper 8 (Correlation Space, reserved) — **placeholder** (reserved slot, scope TBD).

**Book III — The Framework in the World** (Ch 14–16).
- Ch 14 Where Standard Physics Lives. Regime-scoping of quantum mechanics, spacetime geometry, thermodynamics, action principles, classical limit. What the framework does not change.
- Ch 15 Status, Predictions, Falsifiers. Codebase v6.7 scorecard (312 theorems, 39 tested, 36/39 within 5%, 3.83% mean error). Representative predictions table including Ω_Λ = 42/61, Ω_m = 19/61, sin²θ_W = 3/13, Cabibbo. Falsifiers organized by layer attacked. What would not refute the framework. Open problems.
- Ch 16 How to Read the Papers. Map of the series. Reading orders by audience. Status labels. How to critique productively — five specific attack surfaces. Closes with the canonical sentence.

**Canonical sentence placements (load-bearing):**
- Title-page epigraph
- Book I Ch 4.2 (inside the PLEC canonical-statement tcolorbox)
- Book III Ch 16.6 (closing)

**What changed vs. v2.0:**
- v2.0 was a "gateway, not a textbook" five-act narrative (Preface + Act I The Problem + Act II The Reframe + Act III The Map + Act IV Where Standard Physics Lives + Act V The Reading Guide + Appendix Reference Table).
- v3.0 is the ontology paper, not a gateway. Four-layer hierarchy stated explicitly. PLEC named (v2.0 never used the word). Canonical sentence placed (v2.0 did not carry it). Paper 7 at v2.0 scope (spectral action / SM Lagrangian), not v1.0 scope (minimum action only). Paper 1 correctly described as carrying the full quantum reconstruction (v2.0 had said "does NOT derive any specific physics," now false under v4.0-PLEC). MD kept structurally separate from A1 (v2.0 elided this). Descriptive-framing discipline applied sentence-level, not via a global flag. Three-book structure replaces the five-act container so Book II can expand as technical papers mature.

**Audit driver (2026-04-17):** Ran full ontological audit against Paper 1 v4.0-PLEC + Paper 2 v5.3-PLEC + both Technical Supplements + Axiom Inventory v2.1. Found 9 blocker-level divergences (PLEC absence, keystone/A1 conflation, Paper 7 scope, Paper 1 scope, MD-as-separate, canonical sentence, five-vs-six layer inconsistency, etc.), 12 major, 4 minor. Ethan directed a full rebuild at v3.0 scope rather than a surgical port — "this is the ontology paper, words matter a lot, no shortcuts, take the time to craft precise language based on the math."

**Status notes:**
- v3.0 is shippable. Book II Chs 9–11 and 13 carry `pendingbox`-styled placeholders that will fill in as Papers 4, 5, 6, 8 mature; the structure is forward-compatible.
- pdflatex × 3 passes clean. No undefined references, no errors, 12 cosmetic overfull hboxes (description-list labels with long code refs, Paper 1 style inheritance).

**Archive:** `Old/Reference - What Physics Permits v2 Complete_pre-v3.docx` (the Feb-2026 .docx v2.0 that v3.0 replaces).

**Source:** `Papers/Paper 23 - What Physics Permits/Paper_0_WPP_v3.0.tex` + `.pdf`.

**Zenodo DOI (v2.0 deposit):** https://doi.org/10.5281/zenodo.18605692. A v3.0 deposit will follow once Ethan reviews the compiled PDF.

## See also
- [[Paper 1 - Spine]]
- [[Paper 13 - Minimal Admissibility Core]]
- [[PLEC Rollout Plan]]
- [[Derivation Chain]]
- [[Axiom A1]]
- [[Principle - Minimum-Cost Ontology]]
