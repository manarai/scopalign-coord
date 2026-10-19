# O1i review

reviewed_commit: f3d5384 (task/O1i; deliverables committed by the planner from the finished run, unchanged)
date: 2026-10-18
verdict: ACCEPT — outcome: E5r **fails** (Test R and Test F)
recomputation: light — gated Test R numbers from `testR.csv` (exact).

## Checks

Checker PASS. `symmetry_max` 2.1e-7. No deviations.

## Results (n 2759, 3 replicates)

| k | err_k (median) | sub_k (median) | eigen gap λk/λk+1 (median, X) |
|---|---|---|---|
| 3 | 0.973 | 0.44 | 1.33 |
| 5 | **0.934** | 0.84 (max 0.98) | 1.11 |
| 8 | 0.951 | 0.94 | 1.10 |

`frac_robust` 1.0 (starts agree — on the wrong answer). `L(Q̂)` is 15–21×
`L(Q*)` in every case: the restricted family `Q(R)` cannot reach the truth.
Test F: `mix_F` 0.265 (no rotation 0.312), `frame_ratio_F` 2.33,
`sil_F` 0.375 → fails.

## Why (planner diagnostic, `planner_subspace_diag.py`, case r 0)

Subspace transfer `sub_k` (0 = the copy's subspace maps exactly onto the
other's under `Q*`) for candidate statistics:

| statistic | k = 3 | k = 5 | k = 8 |
|---|---|---|---|
| drift second moment `Ū` | 0.77 | 0.95 | 0.97 |
| operator second moment `J̃` | 0.46 | 0.90 | 0.88 |
| mean operator | 0.42 | 0.85 | 1.00 |
| `S` (O1i) | 0.37 | 0.98 | 0.99 |
| position covariance | **0.07** | **0.11** | 0.29 |

Operator statistics have nearly flat spectra (gaps 1.03–1.33), so their top
eigenvectors are not reproducible between copies. The planner's choice of
subspace (`S`), not the restriction idea, is what failed. The position
covariance defines a transferable subspace for k ≤ 5.

## Next

A restricted estimator whose subspace comes from each copy's own position
covariance needs a rule 13 amendment (Tom's decision); otherwise O1c per
rule 14.
