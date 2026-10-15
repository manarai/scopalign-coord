# K0 review

reviewed_commit: d360565 (task/K0)
date: 2026-10-15
verdict: ACCEPT — selection outcome: `selected_layer = "none"` (back to Tom)

## Checker

`check_report.py tasks/K0.md reports/K0/report.json --prev reports/T01e/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/K0/`. Input sha256
f449ff0a… recorded. Task files on main unchanged since issue.

## Recomputation

Planner's own end-to-end script (split, gene filter, HVGs, three layers,
tangent ranks and bases, six drift fits, restricted and full spectra,
classes, draws, all statistics, selection):

- Split, gene list and all six latents: **identical** to the deliverables.
- Pooled draw indices: **identical** to `draw_indices.npz`.
- All 46 numeric gated values agree. Eight spectral ratios/separations
  differ at ≤ 1.0e-5 relative: the executor stores eigenvalues as float32
  (planner: float64; eigenvalue differences 3e-6 to 3e-4). Recomputing the
  statistics with the planner's own draw code on the executor's stored
  float32 spectra reproduces every value to ≤ 1.4e-7. No decision threshold
  lies inside these differences (closest: `spectral_ratio_fa30` 3.01865
  executor vs 3.01866 planner, both > 3.0).

Figures: `fig_spectra.png` (3 configs × 6 tested types, H1 vs H2, null
insets) and `fig_basin.png` (per-type median basin_index H1 vs H2 per
config) are as specified.

## Decision rule outcome

Shared validity: count match 0; no zero-count HVGs; 6 tested types (B cells,
DC, ILC, Megakaryocytes/platelets, Monocytes, T cells).

| config | symmetry | rank ≥ 10 | min separation | spectral | basin | class TV | outcome |
|---|---|---|---|---|---|---|---|
| pca30 | ≤ 1.9e-7 | 0.976 | 5.58 | 1.716 | 1.790 | 1.875 | **inconclusive** (rule 6) |
| fa30 | ≤ 1.1e-7 | 0.814 | 4.65 | **3.019** | 1.427 | 2.395 | **fails** (rule 4) |
| fa10 | ≤ 1.3e-7 | (no gate) | 7.09 | **7.153** | 1.505 | 2.137 | **fails** (rule 4) |

Selection: no config survives; not every eligible config fails (pca30 is
inconclusive) → `selected_layer = "none"`; the planner brings the
per-config outcomes to Tom (as pre-registered).

## Deviations

None listed. None found.

## Verdict

**ACCEPT.** No layer selected; K1 is not written until Tom decides.

## Planner reading (not findings)

- Every config separates cell types well within a half (min separation
  4.6–28.7) and all symmetry gates pass: the operators carry cell-type
  structure. This was never true of the T01-era design.
- **PCA-30 is the closest to reproducible.** Its ratios (1.72 / 1.79 / 1.87)
  sit between the "survives" (≤ 1.5) and "fails" (> 3) lines; its spectra
  overlap closely by eye, its eigenvalue scale is the same in both halves
  (median |λ| 0.336 vs 0.333) and its class counts match (holder 4069 vs
  3945). FA-30 fails narrowly (3.02); FA-10 clearly (7.15).
- **What the ratio measures.** The null compares two samples of one half
  under **one** drift fit; the pooled comparison crosses **two** fits (H1,
  H2). So the ratio contains drift-fit randomness that the null never sees.
  PCA-30's ~1.7 cannot by itself say whether H1 and H2 differ beyond what
  refitting the same cells with another seed would give. That fit-noise
  baseline was not part of K0's design (a planner gap, not an executor one).

## Quarantined (not evidence)

- Observation: input sha256 (confirmed).

## Gated numbers (as reported; recomputed above)

| name | value |
|---|---|
| `basin_ratio_fa10` | 1.5045017257755493 |
| `basin_ratio_fa30` | 1.4269434779310777 |
| `basin_ratio_pca30` | 1.7902583378565315 |
| `basin_sep_h1_fa10` | 7.4132886790852 |
| `basin_sep_h1_fa30` | 4.64836922735241 |
| `basin_sep_h1_pca30` | 5.742357362492693 |
| `basin_sep_h2_fa10` | 7.0873374618521625 |
| `basin_sep_h2_fa30` | 6.001543232084296 |
| `basin_sep_h2_pca30` | 5.5761102144414 |
| `class_tv_ratio_fa10` | 2.136610290440111 |
| `class_tv_ratio_fa30` | 2.394813423187858 |
| `class_tv_ratio_pca30` | 1.8746014871369239 |
| `count_match_max_abs` | 0 |
| `frac_rank_ge10_fa10` | 0.3085357013410656 |
| `frac_rank_ge10_fa30` | 0.8136100036245016 |
| `frac_rank_ge10_pca30` | 0.9757158390721276 |
| `n_cell_types_tested` | 6 |
| `n_cells_h1` | 5518 |
| `n_cells_h2` | 5518 |
| `n_cells_qc` | 11043 |
| `n_genes_detected_both` | 18082 |
| `n_hvg_zero_h1` | 0 |
| `n_hvg_zero_h2` | 0 |
| `rank_h1_fa10` | 9.0 |
| `rank_h1_fa30` | 12.0 |
| `rank_h1_pca30` | 20.0 |
| `rank_h2_fa10` | 9.0 |
| `rank_h2_fa30` | 12.0 |
| `rank_h2_pca30` | 20.0 |
| `selected_layer` | none |
| `spectral_ratio_fa10` | 7.153349751698642 |
| `spectral_ratio_fa30` | 3.0186477939083467 |
| `spectral_ratio_full_fa30` | 4.710526750099895 |
| `spectral_ratio_full_pca30` | 2.328833799569238 |
| `spectral_ratio_pca30` | 1.7164611382311468 |
| `spectral_sep_h1_fa10` | 28.316617237383564 |
| `spectral_sep_h1_fa30` | 11.986176201735196 |
| `spectral_sep_h1_pca30` | 23.37868863118697 |
| `spectral_sep_h2_fa10` | 28.667597821699562 |
| `spectral_sep_h2_fa30` | 16.025976980544748 |
| `spectral_sep_h2_pca30` | 21.75086046954761 |
| `symmetry_h1_fa10` | 1.2901113422003618e-07 |
| `symmetry_h1_fa30` | 1.1325257176508785e-07 |
| `symmetry_h1_pca30` | 1.8742955631978475e-07 |
| `symmetry_h2_fa10` | 1.2580034791653816e-07 |
| `symmetry_h2_fa30` | 9.608473305888765e-08 |
| `symmetry_h2_pca30` | 1.9434391873573986e-07 |
