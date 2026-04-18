---
type: paper
domain: apf
layer: 3-ledgers
created: 2026-04-14
updated: 2026-04-18 (lit-review remediation v3.0 → v3.1; Technical Supplement v1.0 built)
sources: []
---

# Paper 3 - Ledgers
## "Admissible Dynamics, Entropy, and the Arrow of Time"

**Status:** **v3.2 main (32pp with wayfinder) + Technical Supplement v1.1 (15pp)** (Package E lit-review remediation + wayfinder convention, both 2026-04-18). Series-wide "Operational primitives" wayfinder tcolorbox added at start of Part I mapping five primitives (ledger, CPTP map, committed capacity / cross-interface correlation, admissibility carrying capacity ACC, constraint functional Φ_c) to their definition and formalization sites. Five packages landed in the remediation: (A) 35→40-reference bibliography in main + 16→19 in Supplement; (B) PLEC/A1 consistency sweep; (C) §12 Summary + numerical-corroboration table; (D) Technical Supplement mirroring Paper 1/2/7 layout; **(E) dedicated `Relation to prior work` section in main paper** positioning each major claim against its mainstream comparator with explicit reinterpretation framing, plus comparison sub-sections in Supplement §§4, 6, 7, plus T2K F6 falsifier, plus softening of §12 Summary to acknowledge reinterpretation status. Five new references: Spohn 1978, Gibbons-Hawking 1977, Bény-Osborne 2015, Jarlskog 2005/2012, T2K 2020.

**Changelog:**
- **v3.2 main + Supplement v1.1 (2026-04-18):** Package E — Real literature-review remediation, driven by the observation that the v3.1 Package A response had conflated "bibliography" with "literature review" and shipped only the former. New §3 "Relation to prior work" in main paper (2–3pp) with seven sub-sections positioning each major claim against its mainstream literature (CPTP/Stinespring-Choi-GKSL-Kraus; entropy/Lieb-Ruskai-Lindblad-Spohn; horizon/Bekenstein-Hawking-Gibbons-Hawking-Jacobson-holography; Fisher/Braunstein-Caves-Petz-Bény-Osborne-Amari; reconstruction/Hardy-CDP-Masanes-Müller-Dakic-Brukner; arrow-of-time/Maccone-Caticha-Chen-Di Biagio-Buchholz; CP/Jarlskog-T2K). Five new references: Spohn 1978, Gibbons-Hawking 1977, Bény-Osborne 2015, Jarlskog 2005/2012, T2K 2020. Supplement mirror: new §4.4 Spohn comparison, §6 Bekenstein-Hawking-Gibbons-Hawking positioning, §7 Bény-Osborne comparison. §12 Summary softened to acknowledge reinterpretation status throughout. New falsifier F6 anchoring PMNS-phase claim to T2K data. Archive: `Old/Paper_3_LEDGERS_v3.1_pre-PackageE.*`, `Old/Paper_3_Technical_Supplement_v1.0_pre-PackageE.*`. Lit review preserved at `Reference - Paper 3 Lit Review 2 (2026-04-18).md`. Also in this session: new `audit-response` skill authored at `__APF Library/skills/audit-response/SKILL.md` to discipline future responses to external audits.
- **v3.1 + Supplement v1.0 (2026-04-18):** Lit-review remediation pass (bibliography, PLEC consistency, conclusion, Supplement v1.0). See prior Log entry. Archive: `Old/Paper_3_LEDGERS_v3.0_pre-review-edits.tex` + `.pdf`. Lit review 1 preserved at `Reference - Paper 3 Lit Review (2026-04-18).md`. Technical Supplement first appeared at this version.
- **v3.0 (2026-04-17):** Scope expansion sweeping in content from Early Papers 3/34/35/36/39 and the 2026-04 partition-function / saturation-exactness investigations. Added: (a) ledger-accounting principle paragraph in §5 Entropy — entropy as multiplicity count, robust under cost-landscape corrections; (b) new Part II closing section "Geometric Carrying Capacity" (4 subsections) — ACC extends enforcement cost to geometric interfaces, 1/2 × 1/2 horizon bookkeeping factors, unification of thermodynamic / entanglement / horizon entropy; (c) new Arrow-of-Time subsections — macroscopic arrow from loss of global admissible coordination, low-entropy = restricted phase space rebutting typicality; (d) new Part III section "Cosmogenesis and Staged Relaxation" (5 subsections) — constraint functional Φ_c = ln(d_field/d_admissible), overconstrained initial regime, staged relaxation → scale-dependent imprints, trapped-surface threshold, horizon enforceability preview. v1.2 coderef completions folded in: T_deSitter_entropy (gravity.py), L_crossing_entropy (generations.py), L_Fisher_entropy_budget (generations.py); Part I T_tensor, L_nc, L_loc coderefs explicit. Two new falsifiers F4 (ringdown memory) and F5 (structured pre-saturation imprints). Appendix C updated. Compiled clean: 23 pages, 407789 bytes. Archive: `Old/Paper_3_LEDGERS_v1.1-PLEC.tex` + `.pdf`.
- **v1.1-PLEC (2026-04-17):** Surgical PLEC propagation. Abstract + §1 intro reframed around PLEC (Principle of Least Enforcement Cost) with subordination sentence preserved. T_κ upper-bound attribution corrected: was "A1's minimality principle," now "A2, the argmin component of PLEC." Appendix B assumption inventory expanded from "A1, MD, BW" to include A2 as the fourth PLEC component, matching Papers 1 v4.0-PLEC and 2 v5.3-PLEC. No theorems, proofs, numerics, or code references changed. Compiled clean, 16 pages.
- **v1.0:** Pre-PLEC baseline (archived as `Old/Paper_3_LEDGERS_v1.0_pre-PLEC.tex/.pdf`).

**Core claim:** Admissible dynamics = CPTP maps (as a theorem, not a definition). Entropy = committed capacity, with κ = 2 derived from Cauchy uniqueness. Second law, arrow of time, and CP violation follow from a single irreversibility lemma L_irr. Horizon entropy, low-entropy initial conditions, and cosmogenesis follow from the same ledger accounting extended to geometric and cosmological interfaces.

**Key theorems (v3.0 inventory — 17 entries):**
- Part I: T_tensor, L_nc, L_loc, T_CPTP, T_seq
- Part II (Ledger): L_irr, L_irr-uniform, T_κ, T_entropy, T_zeroth, T_first, T_second_law
- Geometric Carrying Capacity: T_deSitter_entropy (gravity.py) — S_dS = 61 ln 102
- Part III (Time & Structure): T_CPT, L_Fisher, L_crossing_entropy, L_Fisher_entropy_budget
- All [P] (proved from A1 with the PLEC four-component inventory).

**Content (v3.0):**
- Part I: admissible dynamics = CPTP maps, sequential composition, no-cloning/deleting/bit-commitment
- Part II: the ledger, κ = 2 from Cauchy uniqueness, entropy as committed capacity, thermodynamic laws, second law (subsystem + cosmological S_dS = 61 ln 102 ≈ 282.12 nats), ledger-accounting principle (entropy as count, robust under cost-landscape corrections)
- Part II closing: Geometric Carrying Capacity — ACC as general interface capacity, 1/2 × 1/2 horizon bookkeeping yielding S = A/4, horizon entropy as saturation, three-entropies-one-accounting unification
- Part III: arrow of time (L_irr quantum-level + macroscopic coordination-loss complement), low-entropy = restricted phase space, T_CPT with maximal T-violation at π/4
- Cosmogenesis section: Φ_c constraint functional, overconstrained initial regime, staged relaxation, trapped-surface threshold, horizon enforceability (Paper 39 / Paper 6 preview)
- Fisher information metric, RG flow as gradient descent, channel-crossing lemmas feeding Paper 4
- CP violation: CKM holonomy + PMNS entropy-optimization two-mechanism DAG
- Five falsifiers (F1 positive-not-CP, F2 temporal Tsirelson, F3 entropy-distinguishability decoupling, F4 ringdown memory, F5 structured pre-saturation imprints)

**Implementation:** core.py (T_CPTP, T_seq, T_κ, T_entropy, L_irr, L_Fisher) + supplements.py (T_zeroth, T_first, T_second_law, T_CPT) + gravity.py (T_deSitter_entropy) + generations.py (L_crossing_entropy, L_Fisher_entropy_budget)

**Status notes:**
- v3.1 main + Supplement v1.0 both compile clean (pdflatex × 3 each, no undefined refs)
- Scope now extends from quantum interfaces to causal and cosmological interfaces
- Paper 7 cross-reference live for thermodynamic formalism (partition function, finite-T corrections)
- Bibliography now in place (35 references in main, 16 in Supplement)
- Three additional non-PLEC hypotheses (H1, H2, H3) explicitly flagged at point of use in the Supplement red-team section
- Paper 3 is now self-citably-shippable: PDF + bibliography + standalone conclusion + Technical Supplement + numerical-corroboration table all in place

## See also
- [[Enforcement Budget]]
- [[Paper 1 - Spine]]
- [[Open Problems]]
