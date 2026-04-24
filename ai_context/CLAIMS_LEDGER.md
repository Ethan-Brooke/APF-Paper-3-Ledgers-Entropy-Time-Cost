# Claims Ledger — Paper 3

| # | Claim | Status | Proof location | Code check | Failure mode |
|---|---|---|---|---|---|
| 1 | $T_{\rm CPTP}$: admissibility forces CPTP | nontrivial | Supp §2 | `check_T_CPTP` | non-CPTP dynamics admissible |
| 2 | $T_\kappa$: cost monotone | nontrivial | Supp §3 | `check_T_kappa` | cost decreases on some admissible trajectory |
| 3 | $T_{\rm entropy}$ via Lieb-Ruskai | standard + local framing | Supp §4 | `check_T_entropy` | LR inequality violated |
| 4 | $T_{\rm second\_law}$ / arrow of time | nontrivial | Supp §5 | `check_T_second_law` | reversible admissible example |
| 5 | $S = A/4$ bookkeeping | reparametrisation | Supp §6 | `check_L_BH_cell_count` | cell-count convention wrong |
| 6 | Fisher-metric uniqueness | nontrivial | Supp §7 | `check_L_Fisher_unique` | alternative metric on admissible states |
| 7 | T_ACC formalization | nontrivial (primary receiver) | Supp §7 (v1.2) | `check_T_ACC_unification` | I1-I4 framework incomplete |
| 8 | T_interface_sector_bridge | nontrivial (primary receiver) | Supp §8 (v1.2) | `check_T_interface_sector_bridge` | 42-dim subspace non-unique |
| 9 | F4 ringdown memory falsifier | empirical | Main §F4 | — | null result at LIGO/LISA |
| 10 | F5 pre-saturation imprints | empirical | Main §F5 | — | null result at Planck |
| 11 | Ledger-accounting principle | modest | Main §v3.0 | `check_L_ledger_accounting` | cost-landscape term breaks count equality |

## Attack surface priority

Claims 1, 2, 4, 7, 8. Claims 7+8 are the primary T_ACC + bridge receivers — if they fail, Paper 8's architecture loses a key receiver.

---

*12 bank-registered checks verify this paper's subset in this repo.*
