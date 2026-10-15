# FINDINGS.md — accepted results only

Appended by the planner on ACCEPT. One entry per task. Gated numbers are the
planner's recomputed values, not the executor's reported ones.

## 2026-10-14 — T00 (ACCEPT)

- scjdo2_sha: 949fb9cc285f99219a95b8cdc9001c6e5d1d4c13
- scjdo_tag: scjdo-v56-submission-20260922-r25
- n_failing_tests: 4
- n_cells_v2_raw: 996
- n_cells_v3_raw: 1222
- Decision rule outcome: ACCEPT (both packages and both scJDO2 imports load; tag matches; both .h5ad have qc_pass and cell_type).
- Commit: task/T00 @ 9ca4e14

## 2026-10-14 — T00b (ACCEPT)

- rank_pca_v2: 21
- rank_pca_v3: 19
- rank_fa_v2: 6
- rank_fa_v3: 6
- rank_ldvae_v2: 19
- rank_ldvae_v3: 13
- frac_fa_at_ceiling: 0.0
- Decision rule outcome: FA is the positional layer.
- Commit: task/T00b @ 0bc9c53

## 2026-10-14 — T01 (ACCEPT; outcome INCONCLUSIVE)

- n_cell_types_tested: 4
- symmetry_v2: 1.156e-07
- symmetry_v3: 1.118e-07
- spectral_w2_pooled: 1233.12
- spectral_w2_null: 61.871
- spectral_ratio: 19.93
- basin_w1_pooled: 3.388
- basin_w1_null: 0.5659
- basin_ratio: 5.988
- class_tv_pooled: 0.00603
- class_tv_null: 0.000581
- class_tv_ratio: 10.38
- spectral_sep_v2: 0.733
- spectral_sep_v3: 1.237
- basin_sep_v2: 10.11
- basin_sep_v3: 1.036
- basin_index_rho (secondary): -0.80
- basin_class_agreement (secondary): 1.0
- Decision rule outcome: INCONCLUSIVE (rule 3 — spectral_sep_v2, spectral_sep_v3, basin_sep_v3 < 2).
- Commit: task/T01 @ 2673fe9
