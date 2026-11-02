# O1l — Review

Reviewer: planner (Claude), 2026-11-02. Branch `task/O1l` @ ba3b2ce (report data from 3058b25), Mac worktree `~/Documents/scopalign-manus-O1l`.

## Checker
`python checker/check_report.py tasks/O1l.md reports/O1l/report.json --prev <O1g report>` → **PASS** (0 fail, 0 warn).

## Recomputation
Light recompute: every gated median recomputed exactly from `cases.csv`, and both basin radii from `basin.csv` (per-case radius = largest t with ≥ 4/5 returns). Not recomputed: screen, LBFGS.

| gated | value |
|---|---|
| n_cases / inputs_match | 72 / 1 |
| l_star_repro_max | 1.4e-45 |
| negctrl_err_shape_A_s | 0.999 |
| errA_s_max | 0.464 (only A, γ 1.6, n 1000; every other A cell ≤ 2e-6) |
| errB_s n3000 γ 1.15 / 1.33 / 1.6 | 0.081 / 0.041 / 0.026 |
| errB_s n10000 γ 1.15 / 1.33 / 1.6 | 0.038 / 0.018 / 0.014 |
| basin_radius_A / B | 0.2 / 0.4 |

Views (γ 1.33, n 3000, r 0): no rotation gives mix 0.00 / sil 0.14; E5s gives mix 0.999 / sil 0.320; the known answer gives mix 0.998 / sil 0.320.
Deviation (declared): in-process threads (≤ 3) instead of spawned workers; local cache directories. No effect on values.
Process note: issued by Tom before O1k was formally accepted. Inputs were verified equal to O1k's (inputs_match 1; L* reproduced to 1e-45), so nothing depends on the order.

## Decision rule outcome
Validity: 72 cases, inputs match, L* reproduced, negative control 1.00 ≥ 0.30 → valid.
Rule 1: errA_s_max 0.46 > 0.10 → **search_unfixed**.

## Verdict
**ACCEPT** — outcome **search_unfixed**, as registered. The label rests on one cell: exact fields, γ 1.6, n 1000, all three replicates (err 0.78 / 0.45 / 0.46, L 0.6–5.3 vs L* 1e-30). In the other 8 exact-field cells with γ > 1, the axis-sign start recovers Q* to 1e-6.

Reading (not the outcome):
- **Fitted fields carry the rotation.** From axis starts, every fitted-field cell with γ > 1 and n ≥ 3000 is ≤ 0.10. At n = 10000 the errors shrink (0.04 / 0.02 / 0.01), so the random-start misses in O1k are fixed. At n = 1000: γ 1.15 gives 0.14, γ 1.33 gives 0.07, γ 1.6 gives 0.05 (one replicate at 0.45). Had rule 1 passed, γ_min = 1.15 → pbmc_range_ok.
- **The basin is wider with fitted fields** (0.4) than with exact ones (0.2; 0.1 in 2 of 3 replicates at γ 1.6). Starting from the axes works because it lands within that radius.
- The single failing cell needs both very elongated wells and few cells, so the sample axes are least accurate exactly where the exact-field basin is narrowest.

## Quarantined (not evidence)
- Planner seeds: at d = 30 the 2³⁰ sign screen is not feasible; a greedy or flow-based sign choice, or O1j's design (top-k block from positions plus operator signs), is the transferable version. This is synthetic, d = 10, Gaussian wells and a pure rotation; PBMC copies also differ by more than a rotation.
