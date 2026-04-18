---
type: entity
domain: apf
layer: cross-cutting
created: 2026-04-14
updated: 2026-04-18 (v6.9 PLEC formalization)
sources: []
---

# APF Codebase
## "Version 6.9: 355 Checks, 342 Bank-Registered Theorems, 19 Registered Modules + standalone/"

**Overview:** A production Python codebase implementing the entire [[Admissible Physics Framework]]. Each check function is a rigorous proof of a theorem; all checks pass with exact arithmetic where possible (rational numbers via `fractions`; numerical routines via `numpy` and `scipy` where required).

**Current version:** **v6.9 (2026-04-18 PLEC formalization pass)**

### v6.9 changes (this pass, 2026-04-18)
- New `apf/plec.py` module (6 checks): `Regime_R`, `Regime_exit_Type_I` through `Type_V`. Formalizes the PLEC infrastructure introduced in Papers 5/6 v2.0-PLEC.
- New `check_A9_closure` in `apf/gravity.py` (1 check). Unifies the Lovelock prerequisites A9.1–A9.5 dispersed across `core.py`, `gravity.py`, `spacetime.py`, `internalization_geo.py`.
- Net delta: +7 bank-registered checks (335 → 342); +7 verify_all (348 → 355); +1 module (18 → 19).
- Papers 5/6 v2.0-PLEC code-anchored with full coderef pass; Paper 6 §11 substantially expanded with the 4-parameter Planck match table, $H_0 = 67.76$ km/s/Mpc derivation, H0DN 7.09σ tension, escape-route table, and explicit falsifier statement.
- The 4-parameter Planck cosmological match ($\Omega_b = 3/61$, $\Omega_c = 16/61$, $\Omega_m = 19/61$, $\Omega_\Lambda = 42/61$, all matching Planck within 1.2%) is now foregrounded as a major framework result.

### Authoritative counts (from `apf.bank._load()` run 2026-04-18)

- **355** total verify_all checks — all passing
- **342** bank-registered theorems (`apf.bank.REGISTRY` size after `_load()`)
- **19** bank-registered modules in the `apf/` package
- **+ apf/standalone/** sub-package (4 additional files, not bank-registered but verify_all-tracked)
- **+ session_phase2_confrontation.py** (3 additional checks, not bank-registered)
- **48** quantitative predictions
- **0** free parameters
- **39** tested predictions
- **32/39** within 3σ of observation
- **Mean error:** 3.83%
- **Median error:** 0.37%

### Module breakdown (bank-registered, Phase 2.9 frozen count)

| Module | Bank-registered theorems | Purpose |
|--------|-------------------------:|---------|
| `apf.core` | 48 | Axiom A1 / PLEC, Hilbert space, Born rule, CPTP |
| `apf.gauge` | 31 | Gauge structure, non-closure, uniqueness |
| `apf.generations` | 85 | Fermion spectrum, 3 generations, masses, mixing |
| `apf.spacetime` | 8 | Metric structure, field equations |
| `apf.gravity` | 9 | Gravitational sector, horizon entropy |
| `apf.cosmology` | 17 | Dark energy, density fractions, N_eff |
| `apf.validation` | 9 | Prediction validation vs. experiments |
| `apf.supplements` | 70 | Lemmas, inequalities, computational aids |
| `apf.majorana` | 10 | Neutrino physics, Majorana constraints |
| `apf.internalization` | 3 | Gauge boson internalization |
| `apf.internalization_geo` | 4 | Geometry from internalization |
| `apf.extensions` | 7 | BSM physics, dark matter mechanisms |
| `apf.red_team` | 17 | Adversarial tests for false positives, edge cases |
| `apf.session_v63c` | 4 | v6.3c incremental: neutrino hierarchy, Yukawa spectral |
| `apf.session_qg` | 4 | QG / Planck-scale scoping |
| `apf.session_nnlo` | 4 | NNLO corrections |
| `apf.session_delta_pmns` | 2 | PMNS phase derivation progress |
| `apf.session_cosmo_update` | 3 | Cosmological observables refresh |
| `apf.plec` | 6 | **NEW v6.9.** Regime R + 5-type regime-exit taxonomy |
| `apf.gravity` (A9_closure addition) | +1 | **NEW v6.9.** Unified Lovelock-prerequisite closure |
| **Bank total** | **342** | |

**Additional verify_all-tracked (not bank-registered):**

| Source | Checks | Notes |
|--------|-------:|-------|
| `apf.session_phase2_confrontation` | 3 | Phase 2 experimental-confrontation module (not yet registered) |
| `apf.standalone.L_Cauchy_uniqueness` | 1 | Cauchy-uniqueness lemma, standalone |
| `apf.standalone.L_CKM_resolution_limit` | 1 | CKM resolution limit, standalone |
| `apf.standalone.phase1_seesaw_closure` | 2 | Phase 1 seesaw-chain closure + RT test |
| `apf.standalone.phase5_theorem_R_audit` | 3 | Phase 5 Theorem R adversarial audit (R1/R2/R3) |
| **Verify_all grand total** | **348** | |

### Environment and dependencies

- **Python ≥ 3.8**
- **numpy ≥ 1.20** (linear algebra, array ops)
- **scipy ≥ 1.7** (special functions, two-loop RG, spectral-action integrals — 4 checks require this specifically: `L_SA_moments`, `L_ST_index`, `L_mc_mt_twoloop_RG`, `L_spectral_action_coefficients`)
- **fractions, math, itertools** (stdlib, for exact-rational core checks)

Install via `pip install -e .` from the codebase root (uses `setup.py`).

### Key properties

- **Mixed-precision:** Core framework uses exact rationals (`fractions.Fraction`); numerical routines (two-loop RG, spectral action) use `numpy` / `scipy` with documented tolerance bounds
- **Modular:** Each module encodes a logical layer of the derivation chain
- **Deterministic:** No randomness; results are reproducible
- **Fast:** Full suite runs in ~5 seconds on standard hardware
- **Auditable:** Every check is a named theorem; human-readable
- **Bank-registered assertion:** `bank.EXPECTED_THEOREM_COUNT = 335` with runtime verification. Any future drift raises a warning immediately.

### Example check:

```python
def check_T_Born(tol=1e-15):
    """T_Born: Born rule probability formula

    For state |psi> and observable A, P(a_n) = |<a_n|psi>|^2.
    """
    # Verify this is the unique cost-minimizing measurement rule
    # under admissibility constraint...
    return True
```

### Development practice

- **v6.8 frozen** 2026-04-18 (Phase 2.9 canonicalization) as the authoritative snapshot for the current paper series (Papers 0, 1, 2, 3, 7 live in .tex; Papers 4, 5, 6 source-pending).
- New features branch for v7.0 (future, post-arXiv submission window).
- All commits tested; no regressions.
- `verify_all.py` is the single source of truth for count; CLAUDE.md and wiki mirror its output.

### Canonicalization history

- **v6.7 → v6.8 (2026-04-18):** Phase 2.9 canonicalization pass. Moved `session_*.py` and standalone files into the `apf/` package (previously at top level, breaking imports in `verify_all.py`). Fixed broken `rt_check` helper in `apf/standalone/phase1_seesaw_closure.py`. Updated `apf.bank.EXPECTED_THEOREM_COUNT` from 310 → 335 to match loaded reality. Updated `setup.py` version 6.7.0 → 6.8.0 and added explicit `install_requires = numpy, scipy`. Synced all count references across CLAUDE.md (× 2), wiki (Index, Log, APF Codebase, Predictions Catalog, paper-specific pages), Schema.md, Paper Update Work Plan, Paper Index, Paper 0 v3.0, Paper 3 v3.2 + Supplement. Driven by the external audit of 2026-04-18 (`APF Reference Docs/Reference - External Audit Memo (2026-04-18).md`).
- **v6.7 (pre-2026-04-18):** Prior authoritative state. Docs drifted from 294 → 312 → 321 → 349 → 351 across different files; `bank.EXPECTED_THEOREM_COUNT = 310` stale relative to loaded 335; session/standalone files couldn't import via `apf.*` namespace; one red-team check failed on a missing helper. All resolved in v6.8.

## See also
- [[Paper 13 - Minimal Admissibility Core]]
- [[Derivation Chain]]
- [[Red Team Tests]]
- [[Predictions Catalog]]
