# T02c review

reviewed_commit: 0dec531 (task/T02c)
date: 2026-10-15
verdict: ACCEPT — Part 1 outcome: **frame agreement** (rule 1)
recomputation: light — Part 1 recomputed from `centroid_ops.npz` and
`reports/T02b/loadings.npz` (`resid_identity_10`, `resid_load_10`,
`resid_null_10`, `frame_ratio_load`, `subspace_min_sv_10`: all exact).
Part 2 and views checked against the task, not recomputed.

## Checks

Checker PASS (0 fail, 0 warn); branch changes only `reports/T02c/`;
`k2_repro_max_abs` 0.0; no deviations. Figures as specified.

## Part 1 (decides)

3 clusters (B cells, Monocytes, T cells). Top-10 residual: no rotation
0.607 → gene-anchored rotation 0.314; refit null 0.288.
`frame_ratio_load` **1.09 ≤ 1.5 → frame agreement**.

## Part 2 (pre-registered reading)

`n_starts_near_best_R` 1 (< 2) → **operators validate the rotation but do
not estimate it reproducibly at this cell count**. For the record:
`q_op_vs_load_10` 0.229, `r_split` 0.112, `r_seed` 0.640 (fit-seed
variation dominates).

## Views

| view | mix | sil |
|---|---|---|
| joint PCA (uncorrected) | 0.111 | 0.413 |
| separate PCAs, no rotation | 0.643 | 0.418 |
| separate PCAs, gene rotation | **0.899** | 0.421 |

Reading met: mixing rises (0.64 → 0.90) and cell-type separation is kept
(0.421 ≥ 0.9 × 0.418).

## Verdict

**ACCEPT.** Frame agreement under the gene-anchored rotation.

## Planner reading (not findings)

- With the frames lined up by the shared genes, the v2 and v3 wells of the
  three big clusters point the same way within drift-fit noise (1.09×).
  Together with K2 (eigenvalues agree), the axiom holds for this pair at
  matched count and depth, in the top 10 PCs.
- `fig_frames.png`: the stiffest eigen-directions agree well after rotation
  (cosines 0.8–0.99 for ranks 8–10 in all three clusters; monocytes high
  throughout); middle ranks of B and T cells agree less, consistent with a
  crowded part of the spectrum.
- The rotation comes from gene loadings, not operators (Part 2). The joint
  PCA view shows the chemistry effect clearly (mix 0.11); the gene-rotated
  separate frames remove it while keeping cell types apart.
- Limits: one pair, ~970 cells, 3 clusters, top 10 PCs, depth-matched.

## Gated numbers (as reported)

- `frame_ratio_load`: 1.09032
- `k2_repro_max_abs`: 0
- `loss_best_R`: 0.0841254
- `mix_joint`: 0.111146
- `mix_sep_load`: 0.898521
- `mix_sep_noA`: 0.643137
- `n_cells_a`: 969
- `n_cells_b`: 969
- `n_clusters`: 3
- `n_starts_near_best_R`: 1
- `outcome`: frame agreement
- `q_op_vs_load_10`: 0.228887
- `r_seed`: 0.639769
- `r_split`: 0.112076
- `resid_identity_10`: 0.607021
- `resid_load_10`: 0.313649
- `resid_null_10`: 0.287666
- `resid_opfit_10`: 0.290044
- `sil_joint`: 0.413136
- `sil_sep_load`: 0.420815
- `sil_sep_noA`: 0.418157
- `subspace_min_sv_10`: 0.539825
- `symmetry_max`: 2.30366e-07
