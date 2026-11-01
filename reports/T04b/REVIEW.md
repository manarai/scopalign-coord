# T04b — Review

Reviewer: planner (Claude), 2026-11-01. Branch `task/T04b` @ 2d7b50b (report data from 161f333), Manus cloud sandbox.

## Checker
`python checker/check_report.py tasks/T04b.md reports/T04b/report.json --prev reports/T04c/report.json` → **PASS** (0 fail, 0 warn).

## Recomputation
Light recompute (working mode): all 13 numeric gated values recomputed from `sweep.csv` (150 rows = 30 conditions × 5 methods) with the task's definitions: identical to the report.
End-to-end spot-check: `reproduction.csv` reruns the two T04 conditions (σ 0.5, global, f 0 and 4) through the full pipeline (data, per-half PCA, all five methods including scVI). Compared with the committed `reports/T04/scores.csv`: max |Δ| 5.0e-16 (scVI 0). The embeddings were not rebuilt for the other 28 conditions.

| gated | report | recomputed |
|---|---|---|
| n_conditions | 30 | 30 |
| repro_t04_max_abs | 5.0e-16 | 5.0e-16 (vs T04 scores.csv) |
| repro_t04_scvi_abs | 0 | 0 |
| frac_rigid_removes_batch | 0.417 | 0.417 |
| frac_harmony_removes_batch | 0.667 | 0.667 |
| frac_scvi_removes_batch | 1.000 | 1.000 |
| frac_rigid_wins | 0.208 | 0.208 |
| frac_rigid_wins_global | 0.333 | 0.333 |
| frac_rigid_wins_celltype | 0.083 | 0.083 |
| median_gain | 0.158 | 0.158 |
| bio_false_max_gene_rot | 0.090 | 0.090 |
| bio_false_max_harmony | 0.041 | 0.041 |
| bio_false_max_scvi | 0.097 | 0.097 |

Figures: `fig_retention.png`, `fig_mixing.png`, `fig_wins.png` show what the task asked for (one panel per batch type × σ; win/lose/no-mix grid).

Deviation: the remaining conditions ran in two spawn workers × 3 Torch threads, per-condition checkpoints. Conditions are independent, and the T04 conditions reproduce exactly, so no gated number can move. Accepted.

## Decision rule outcome
Validity: n_conditions 30, repro 5e-16 ≤ 1e-6 → valid.
frac_rigid_wins 0.208 < 0.25 → **none**.

Wins (5 of 24): global σ 0.25 f 2, f 4; global σ 0.5 f 2, f 4; celltype σ 0.25 f 2.
How the other 19 fail, read from `sweep.csv` (description, not a gated number):
- rigid does not reach mix 0.80 at σ = 1.0 (all 8), at f = 8 (all 6) and in 4 celltype conditions;
- at f = 1 rigid keeps less biology than scVI (gain −0.12 to −0.36).
Where rigid does reach the mixing bar, its retention advantage is real (+0.13 to +0.27). gene_rot keeps more biology than both Harmony and scVI in 16 of 24 conditions, but in most of those it does not remove the batch to the pre-registered bar. scVI removes the batch in all 24.

## Verdict
**ACCEPT** — outcome **none**. By the task's own rule, the T04 advantage is not general: it holds for a mild-to-moderate, cell-type-independent batch with a medium-sized response, and the rigid correction under-mixes elsewhere. Per "what would change my mind", T04 is a case study, not a method claim.

## Quarantined (not evidence)
- Executor observations: the 30 condition matrices and 6 oracles were cached before the sweep; the σ 0.5 global matrices equal T04's cached counts; checkpoints are in the gitignored cache.
- deviations_requested: none.
- Planner note (seed only): at f = 8 Harmony also drops below 0.80 mixing (0.72–0.77), so mix_bg on non-responders may be affected by a strong response. Not tested.
