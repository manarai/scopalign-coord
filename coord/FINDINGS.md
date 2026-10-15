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

## 2026-10-14 — T00c (ACCEPT)

- n_genes_detected_both: 12890
- n_hvg_zero_v2: 0
- n_hvg_zero_v3: 0
- n_fa_tiny_noise_v2: 0
- n_fa_tiny_noise_v3: 0
- rank_pca_v2: 21
- rank_pca_v3: 19
- rank_fa_v2: 11
- rank_fa_v3: 14
- rank_ldvae_v2: 19
- rank_ldvae_v3: 13
- frac_fa_at_ceiling: 0.0
- frac_ldvae_at_ceiling: 0.0
- Decision rule outcome: FA is the positional layer for T01c.
- Commit: task/T00c @ 91a5449

## 2026-10-14 — T01c (ACCEPT; outcome INCONCLUSIVE)

- n_cell_types_tested: 4
- symmetry_v2: 8.69e-08
- symmetry_v3: 1.24e-07
- spectral_w2_pooled: 506.25
- spectral_w2_null: 9.436
- spectral_ratio: 53.65
- basin_w1_pooled: 5.163
- basin_w1_null: 2.243
- basin_ratio: 2.302
- class_tv_pooled: 0.3618
- class_tv_null: 0.05559
- class_tv_ratio: 6.508
- spectral_sep_v2: 2.958
- spectral_sep_v3: 13.77
- basin_sep_v2: 1.221
- basin_sep_v3: 3.530
- basin_index_rho (secondary): 0.20
- basin_class_agreement (secondary): 0.25
- Decision rule outcome: INCONCLUSIVE (rule 3 — basin_sep_v2 < 2).
- Commit: task/T01c @ 13fb1bc

## 2026-10-15 — T01y (ACCEPT; EXPLORATORY — beside T01c, not an axiom test)

- scale_v2: 30.2874
- scale_v3: 4.80332
- scale_ratio: 6.3055
- spectral_w2_pooled_s: 11.778
- spectral_w2_null_s: 0.75438
- spectral_ratio_s: 15.613
- basin_w1_pooled_s: 1.0558
- basin_w1_null_s: 0.19483
- basin_ratio_s: 5.4188
- class_tv_ratio_s: 6.508188
- max_sep_change: 3.61e-06
- Reading-rule outcome: "shape differs" (rule 3).
- Commit: task/T01y @ 9366943

## 2026-10-15 — T01d (ACCEPT; outcome INCONCLUSIVE)

- depth_v2_median: 3366
- depth_v3_median: 7017
- depth_v3thin_median: 3378
- depth_match_max_rel: 0.00726
- n_genes_detected_both: 12512
- n_hvg_zero_v2: 0
- n_hvg_zero_v3: 0
- rank_fa_v2: 12
- rank_fa_v3: 15
- frac_fa_at_ceiling: 0.0
- n_cell_types_tested: 4
- symmetry_v2: 9.46e-08
- symmetry_v3: 8.59e-08
- spectral_w2_pooled: 301.07
- spectral_w2_null: 9.804
- spectral_ratio: 30.71
- basin_w1_pooled: 4.409
- basin_w1_null: 2.498
- basin_ratio: 1.765
- class_tv_pooled: 0.5462
- class_tv_null: 0.05678
- class_tv_ratio: 9.620
- spectral_sep_v2: 3.790
- spectral_sep_v3: 13.06
- basin_sep_v2: 1.719
- basin_sep_v3: 2.841
- basin_index_rho (secondary): 0.80
- basin_class_agreement (secondary): 0.0
- Decision rule outcome: INCONCLUSIVE (rule 3 — basin_sep_v2 < 2).
- Commit: task/T01d @ 721ad25

## 2026-10-15 — T01e (ACCEPT; stopped; outcome INCONCLUSIVE)

- n_cells_b_raw: 11769
- n_cells_b_qc: 11043
- n_cells_a: 1113
- n_cells_bsub: 1113
- count_match_max_abs: 0
- depth_a_median: 6607
- depth_b_median: 6605
- depth_match_max_rel: 0.00158
- n_genes_detected_both: 14475
- n_hvg_zero_a: 0
- n_hvg_zero_b: 0
- rank_fa_a: 8
- rank_fa_b: 9
- frac_fa_at_ceiling: 0.0
- Decision rule outcome: INCONCLUSIVE (rule 0 — FA not admissible; run stopped before operators).
- Commit: task/T01e @ cd0223e

## 2026-10-15 — K0 (ACCEPT; selected_layer none)

- basin_ratio_fa10: 1.5045017257755493
- basin_ratio_fa30: 1.4269434779310777
- basin_ratio_pca30: 1.7902583378565315
- basin_sep_h1_fa10: 7.4132886790852
- basin_sep_h1_fa30: 4.64836922735241
- basin_sep_h1_pca30: 5.742357362492693
- basin_sep_h2_fa10: 7.0873374618521625
- basin_sep_h2_fa30: 6.001543232084296
- basin_sep_h2_pca30: 5.5761102144414
- class_tv_ratio_fa10: 2.136610290440111
- class_tv_ratio_fa30: 2.394813423187858
- class_tv_ratio_pca30: 1.8746014871369239
- count_match_max_abs: 0
- frac_rank_ge10_fa10: 0.3085357013410656
- frac_rank_ge10_fa30: 0.8136100036245016
- frac_rank_ge10_pca30: 0.9757158390721276
- n_cell_types_tested: 6
- n_cells_h1: 5518
- n_cells_h2: 5518
- n_cells_qc: 11043
- n_genes_detected_both: 18082
- n_hvg_zero_h1: 0
- n_hvg_zero_h2: 0
- rank_h1_fa10: 9.0
- rank_h1_fa30: 12.0
- rank_h1_pca30: 20.0
- rank_h2_fa10: 9.0
- rank_h2_fa30: 12.0
- rank_h2_pca30: 20.0
- selected_layer: none
- spectral_ratio_fa10: 7.153349751698642
- spectral_ratio_fa30: 3.0186477939083467
- spectral_ratio_full_fa30: 4.710526750099895
- spectral_ratio_full_pca30: 2.328833799569238
- spectral_ratio_pca30: 1.7164611382311468
- spectral_sep_h1_fa10: 28.316617237383564
- spectral_sep_h1_fa30: 11.986176201735196
- spectral_sep_h1_pca30: 23.37868863118697
- spectral_sep_h2_fa10: 28.667597821699562
- spectral_sep_h2_fa30: 16.025976980544748
- spectral_sep_h2_pca30: 21.75086046954761
- symmetry_h1_fa10: 1.2901113422003618e-07
- symmetry_h1_fa30: 1.1325257176508785e-07
- symmetry_h1_pca30: 1.8742955631978475e-07
- symmetry_h2_fa10: 1.2580034791653816e-07
- symmetry_h2_fa30: 9.608473305888765e-08
- symmetry_h2_pca30: 1.9434391873573986e-07
- Decision rule outcome: pca30 inconclusive (rule 6); fa30 fails (rule 4); fa10 fails (rule 4); selected_layer "none" (no config survives, not all fail).
- Commit: task/K0 @ d360565
