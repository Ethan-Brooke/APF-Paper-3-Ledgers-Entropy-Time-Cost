---
type: paper
domain: apf
layer: 5-quantum
created: 2026-04-14
updated: 2026-04-14
sources: []
---

# Paper 5 - Quantum
## "Hilbert Space, the Born Rule, and Tensor Products"

**Status:** Updated PDF available

**Core claim:** [[Axiom A1]] forces quantum mechanics from first principles. The structure of admissible distinction-enforcement necessarily gives rise to complex Hilbert space, the [[Born Rule]], and entanglement via tensor products. Not assumed; derived.

**Key theorems:**
- T2 (`check_T2`) — Hilbert space forced by operability constraints
- T3 (`check_T3`) — Noncommutativity of operators is necessary
- T_Born (`check_T_Born`) — [[Born Rule]] probability formula
- T_Tsirelson (`check_T_Tsirelson`) — Tsirelson bound on correlations
- T_CPTP (`check_T_CPTP`) — Completely positive, trace-preserving maps
- [[Tensor Product Structure]] — Multi-particle Hilbert spaces

**Content:**
- Operability bound → Hilbert space structure
- Why complex numbers, not real or quaternionic
- Proof of the Born rule from enforcement cost
- Superposition and entanglement as cost-efficient encodings
- Quantum non-locality without action-at-a-distance
- CPTP maps and the arrow from pure to mixed states
- Connection to quantum information and von Neumann algebras

**Implementation:** `core.py` (48 checks)

**Status notes:**
- Fully implemented in v6.7
- PDF updated; stable

## See also
- [[Born Rule]]
- [[Axiom A1]]
- [[Paper 1 - Spine]]
