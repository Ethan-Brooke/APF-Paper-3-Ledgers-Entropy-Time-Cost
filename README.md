# Ledgers: Admissible Dynamics, Entropy, and the Arrow of Time

### Interactive Mathematical Appendix to Paper 3 of the Admissibility Physics Framework

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18604844.svg)](https://doi.org/10.5281/zenodo.18604844) [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost/blob/main/APF_Reviewer_Walkthrough.ipynb)

[Interactive Derivation DAG](https://ethan-brooke.github.io/APF-Paper-3-Ledgers-Entropy-Time-Cost/) · [Theorem Map](#theorem-mapping-table) · [Reviewers' Guide](REVIEWERS_GUIDE.md) · [Citation](#citation)

---

## Why this codebase exists

Three results from the foundation of Paper 1, all following from PLEC's four structurally necessary components: complete positivity as a theorem (T_CPTP); entropy as committed capacity (T_entropy + T_kappa with kappa=2 derived); the arrow of time from irreversibility (L_irr + T_second_law). Extends to Geometric Carrying Capacity (ACC), horizon entropy S = A/4, three-entropies unification, and staged cosmogenesis. Now includes a dedicated 'Relation to prior work' section positioning each major claim against its mainstream comparator (Stinespring/Choi/GKSL, Lieb-Ruskai/Spohn, Bekenstein-Hawking-Gibbons-Hawking, Petz/Beny-Osborne, Jarlskog/T2K).

This repository is the executable proof.

The codebase is a faithful subset of the canonical APF codebase v6.8 (frozen 2026-04-18; 348 verify_all checks, 335 bank-registered theorems across 18 modules). Each theorem in the manuscript traces to a named `check_*` function in `apf/core.py`, which can be called independently and returns a structured result.

The codebase requires Python 3.8+ and NumPy / SciPy (some numerical lemmas use them; see `pyproject.toml`).

## How to verify

Three paths, in order of increasing friction:

**1. Colab notebook — zero install.** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost/blob/main/APF_Reviewer_Walkthrough.ipynb) Every key theorem is derived inline, with annotated cells you can inspect and modify. Run all cells — the full verification takes under a minute.

**2. Browser — zero install.** Open the [Interactive Derivation DAG](https://ethan-brooke.github.io/APF-Paper-3-Ledgers-Entropy-Time-Cost/). Explore the dependency graph. Hover any node for its mathematical statement, key result, and shortest derivation chain to A1. Click **Run Checks** to watch all theorems verify in topological order.

**3. Local execution.**

```bash
git clone https://github.com/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost.git
cd APF-Paper-3-Ledgers-Entropy-Time-Cost
pip install -e .
python run_checks.py
```

Expected output:

```
      Paper 3 (Ledgers): 13 passed, 0 failed, 13 total — verified in <minutes>
```

**4. Individual inspection.**

```python
from apf.bank import get_check
r = get_check('check_T_CPTP')()
print(r['key_result'])
```

For reviewers, a [dedicated guide](REVIEWERS_GUIDE.md) walks through the logical architecture, the structural assumptions, and the anticipated objections.

---

## Theorem mapping table

This table maps every result in the manuscript to its executable verification.

| Check | Type | Summary |
|-------|------|---------|
| `check_T_CPTP` | Theorem | T_CPTP: CPTP Maps from Admissibility-Preserving Evolution. |
| `check_T_kappa` | Theorem | T_kappa: Directed Enforcement Multiplier. |
| `check_T_entropy` | Theorem | T_entropy: Von Neumann Entropy as Committed Capacity. |
| `check_T_zeroth_law` | Theorem | T_zeroth_law: Zeroth Law of Thermodynamics [P]. |
| `check_T_first_law` | Theorem | T_first_law: First Law of Thermodynamics [P]. |
| `check_T_second_law` | Theorem | T_second_law: Second Law of Thermodynamics [P]. |
| `check_T_deSitter_entropy` | Theorem | T_deSitter_entropy: de Sitter Entropy from Capacity Microstate Counting [P]. |
| `check_L_irr` | Lemma | L_irr: Irreversibility from Admissibility Physics. |
| `check_L_irr_uniform` | Lemma | L_irr_uniform: Sector-Uniform Irreversibility. |
| `check_T_CPT` | Theorem | T_CPT: CPT Invariance [P]. |
| `check_L_Fisher_measure` | Lemma | L_Fisher_measure: APF Capacity Counting Derives S = (d_eff/2) ln det G [P]. |
| `check_L_crossing_entropy` | Lemma | L_crossing_entropy: 1/α_cross = S_dS / 6 [P]. |
| `check_L_Fisher_entropy_budget` | Lemma | L_Fisher_entropy_budget: Generation Mixing Uses 17.1% of de Sitter Entropy [P]. |

All check functions reside in `apf/core.py`. Every function listed above can be called independently and returns a structured result including its logical dependencies and the mathematical content it verifies.

---

## The derivation chain

```
  Level 0: T_CPTP · T_kappa · T_entropy · T_zeroth_law · T_first_law · T_second_law · T_deSitter_entropy · L_irr · L_irr_uniform · T_CPT · L_Fisher_measure · L_crossing_entropy · L_Fisher_entropy_budget
```

The [interactive DAG](https://ethan-brooke.github.io/APF-Paper-3-Ledgers-Entropy-Time-Cost/) shows the full graph with hover details and animated verification.

---

## Repository structure

```
├── README.md                              ← you are here
├── REVIEWERS_GUIDE.md                     ← physics-first walkthrough for peer reviewers
├── apf/
│   ├── core.py                            ← 13 theorem check functions
│   ├── apf_utils.py                       ← exact arithmetic + helpers
│   └── bank.py                            ← registry and runner
├── docs/
│   └── index.html                         ← interactive derivation DAG (GitHub Pages)
├── APF_Reviewer_Walkthrough.ipynb         ← Colab notebook
├── run_checks.py                          ← convenience entry point
├── pyproject.toml                         ← package metadata
├── zenodo.json                            ← archival metadata
├── Paper_3_Ledgers_Entropy_Time_Cost_v3.2.tex                ← the paper
├── Paper_3_Ledgers_Entropy_Time_Cost_Supplement_v1.1.tex                ← Technical Supplement

└── LICENSE                                ← MIT
```

---

## What this paper derives and what it does not

**Derived:** (see Theorem mapping table above)

**Not derived here:** Specific results outside this paper's scope live in companion papers — see [`APF_Paper_Index.md`](https://github.com/Ethan-Brooke) for the full series.

---

## Citation

```bibtex
@software{apf-paper3,
  title   = {Ledgers: Admissible Dynamics, Entropy, and the Arrow of Time},
  author  = {Brooke, Ethan},
  year    = {2026},
  doi     = {10.5281/zenodo.18604844},
  url     = {https://github.com/Ethan-Brooke/APF-Paper-3-Ledgers-Entropy-Time-Cost}
}
```

## License

MIT. See [LICENSE](LICENSE).

---

*Generated by the APF `create-repo` skill on 2026-04-18. Codebase snapshot: v6.8 (frozen 2026-04-18; 348 verify_all checks, 335 bank-registered theorems, 48 quantitative predictions).*
