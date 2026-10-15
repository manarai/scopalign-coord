# T02d review

reviewed_commit: da59e7f (task/T02d)
date: 2026-10-15
verdict: ACCEPT — outcome: **frame agreement** (rule 1)
recomputation: light — Part 1 recomputed from `floor_ops.npz` and
`reports/T02b/loadings.npz`; `resid_identity_10`, `resid_load_10`,
`resid_null_10`, `frame_ratio_load`, `floor_vs_centroid` exact.

## Checks

Checker PASS; branch changes only `reports/T02d/`; floor = 30 lowest-energy
cells by stable argsort of `drift.energy`, per fit, including the null;
`k2_repro_max_abs` 0.0; no deviations.

## Outcome

Basin-floor operators, top 10 PCs: no rotation 0.598 → gene rotation 0.345;
null 0.287; `frame_ratio_load` **1.20 ≤ 1.5 → frame agreement**.

## Verdict

**ACCEPT.**

## Planner reading (not findings)

- Centroid (T02c, 1.09) and basin floor (T02d, 1.20) agree: both give frame
  agreement; the floor reading is ~10% further from the null
  (`floor_vs_centroid` 1.10). The floor is the density peak (lowest
  energy); the centroid is the mean position.
- The two summaries use largely different cells (overlap 0.38), so the
  agreement is not one set of cells counted twice.

## Gated numbers (as reported)

- `floor_centroid_overlap`: 0.377778
- `floor_vs_centroid`: 1.10223
- `frame_ratio_load`: 1.20179
- `k2_repro_max_abs`: 0
- `n_cells_a`: 969
- `n_cells_b`: 969
- `n_clusters`: 3
- `outcome`: frame agreement
- `resid_identity_10`: 0.597869
- `resid_load_10`: 0.344704
- `resid_null_10`: 0.286825
- `symmetry_max`: 2.30366e-07
