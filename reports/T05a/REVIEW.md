# T05a review

reviewed_commit: 762fbb9 (task/T05a)
date: 2026-10-15
verdict: ACCEPT — descriptive benchmark (no single outcome)
recomputation: light — `mix` and `sil` for all five embeddings recomputed
from `embeddings.npz`; all five `frame_ratio` values recomputed from
`centroid_ops.npz` (all exact). Embedding construction and drift fits
checked against the task, not rerun.

## Checks

Checker PASS (1 warn: lockfile adds harmonypy, the authorised install);
branch changes only `reports/T05a/`; `ours_repro_max_abs` 0.0;
`symmetry_max` 2.4e-07; 3 clusters. Figure as specified.

## Results

| method | mix | sil | frame ratio (full 30-d) | frame outcome |
|---|---|---|---|---|
| joint (uncorrected) | 0.111 | 0.413 | 1.549 | inconclusive |
| centre | 0.705 | 0.437 | 1.364 | agreement |
| Harmony | **1.001** | **0.447** | **1.349** | agreement |
| scVI | 0.983 | 0.082 | 1.416 | agreement |
| ours (gene rotation) | 0.899 | 0.421 | 1.553 | inconclusive |

## Pre-registered readings

- Mixes (≥ 0.80): Harmony, scVI, ours. Keeps types (≥ 0.9 × 0.413 = 0.372):
  joint, centre, Harmony, ours (scVI 0.082 does not, by this measure).
- **Ours competitive?** No: `mix_ours` 0.899 < 1.001 − 0.05, and
  `sil_ours` 0.421 < 0.447 − 0.02.
- **Ours adds operator information?** No: ours is inconclusive (1.553),
  while Harmony, which mixes and keeps types, reaches agreement.
- **"What would change my mind" is met:** Harmony mixes as well (better),
  keeps types as well (better) and reaches operator agreement. **On this
  pair, the gene-anchored rotation offers nothing extra.**

## Verdict

**ACCEPT.**

## Planner reading (not findings)

- The operator frame test, in full 30-d, barely separates methods: even the
  uncorrected joint PCA sits at 1.55, next to ours (1.55), and Harmony /
  centre / scVI at 1.35–1.42. On this easy, same-donor, depth-matched pair
  it is not a discriminating measure. (T02c's 1.09 for ours used the top 10
  PCs; full 30-d adds the noisier low-variance directions.)
- scVI's silhouette (0.08) contradicts its UMAP, which separates types
  cleanly: silhouette in scVI's latent geometry is not on the same scale as
  in PCA-type spaces, so cross-method silhouette comparisons are weak.
- The honest position: on a pair this easy, the operator idea is not needed
  to align cells, and its test does not yet distinguish good from bad
  alignments. Whether it adds anything must be shown on harder data
  (two labs, real biological differences, a known batch: T04 / T03).

## Gated numbers (as reported)

- `frame_outcome_centre`: agreement
- `frame_outcome_harmony`: agreement
- `frame_outcome_joint`: inconclusive
- `frame_outcome_ours`: inconclusive
- `frame_outcome_scvi`: agreement
- `frame_ratio_centre`: 1.36376
- `frame_ratio_harmony`: 1.34944
- `frame_ratio_joint`: 1.5489
- `frame_ratio_ours`: 1.55313
- `frame_ratio_scvi`: 1.41606
- `mix_centre`: 0.704575
- `mix_harmony`: 1.00131
- `mix_joint`: 0.111146
- `mix_ours`: 0.898521
- `mix_scvi`: 0.982766
- `n_clusters`: 3
- `ours_repro_max_abs`: 0
- `sil_centre`: 0.437443
- `sil_harmony`: 0.447127
- `sil_joint`: 0.413136
- `sil_ours`: 0.420815
- `sil_scvi`: 0.0820816
- `symmetry_max`: 2.36e-07
