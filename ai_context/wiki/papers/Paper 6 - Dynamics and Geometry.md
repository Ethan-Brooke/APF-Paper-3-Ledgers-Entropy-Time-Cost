---
type: paper
domain: apf
layer: 6-dynamics
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Paper 6 - Dynamics and Geometry
## "Spacetime, Gravity, and Cosmology"

**Status:** Updated PDF available

**Core claim:** The admissibility framework connects quantum field theory to classical spacetime geometry and cosmology. Derives Einstein's equations, cosmological density fractions (Ω_Λ, Ω_m, Ω_r), neutrino masses, and N_eff predictions from the same principle.

**Key theorems:**
- L_equip (`check_L_equip`) — Equipartition of energy into sectors
- L_dark_budget (`check_dark_budget`) — Cosmological constant forced at observed value
- L_N_eff_prediction (`check_N_eff`) — Effective number of neutrino species
- L_equation_of_state (`check_w_DE`) — Dynamical dark energy equation of state
- Neutrino mass ratios from enforcement cost hierarchies
- [[Cosmological Predictions]] — Ω_Λ, Ω_m, Ω_r, H_0, η_B

**Content:**
- Transition from quantum field theory to spacetime geometry
- Enforcement cost in the gravitational sector
- Why a positive cosmological constant (dark energy)
- Density fractions: matter, radiation, dark energy
- Neutrino physics and oscillation parameters
- Predictions vs. WMAP, Planck, DESI observations
- Open question: fully derive Riemannian geometry from A1

**Implementation:** `gravity.py` (9 checks) + `spacetime.py` (8 checks) + `cosmology.py` (17 checks)

**Status notes:**
- PDF updated; matches v6.7
- Dark energy equation-of-state under investigation with DESI

## See also
- [[Cosmological Predictions]]
- [[Predictions Catalog]]
- [[Open Problems]]
