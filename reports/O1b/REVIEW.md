# O1b review

reviewed_commit: 0df3c07 (task/O1b)
date: 2026-10-16
verdict: ACCEPT — outcome: E5 **fails** (as registered)
recomputation: light — Test R medians from `testR.csv` and `err_F` from
`Q_hat_F.npy` vs `Q_load` (all exact).

## Checks

Checker PASS. Four declared deviations, none able to move a gated number:
parallel starts (deterministic, identical results), stdin redirect for
detached workers, a views-only indexing fix, `fit_symmetry.csv` moved to
`artifacts`. `o1_repro_max_abs` 0.0 (operators reproduce O1's exactly).

## Results

| n | err_R_E5 (median of 3) |
|---|---|
| 500 | 0.727 |
| 1000 | 0.327 |
| 2000 | 0.189 |
| 2759 | 0.179 |

`err_F_E5` 0.582; `n_starts_near_best_F` 1. No `err_R ≤ 0.10` and
`err_F > 0.25` → **E5 fails**.

## Planner's spec error — start seeds overlap the true rotation

O1b's starts use `ortho_group(random_state=20261021 + i)`; Test R's true
rotation uses `random_state=20261031 + r`. For start `i = 10 + r` the start
**is** `Q*`. The best start was 10/11/12 in every n ≥ 1000 case (and 11/12
for two n = 500 cases). So O1b's Test R errors are those reached **from the
truth**: the objective's minimum lies ~0.18 from `Q*` at n ≥ 2000 (a bias
from field-fit noise), and a search without that start can only do the
same or worse. O1's E0 had the same overlap (and failed anyway); O1f used
other seeds. Test F (target `Q_load`) is unaffected: 0.58.

## Views (Test F; executor-made)

| view | mix | sil | frame ratio |
|---|---|---|---|
| separate frames, no rotation | 0.312 | 0.381 | 3.62 |
| **E5 estimate** | **0.808** | **0.391** | **1.22** |
| `Q_load` (known answer) | 0.973 | 0.392 | 1.48 |

## Planner reading (not findings)

- **E5 is the first operator-only rotation that aligns cells**: mixing 0.81
  with cell types kept, and the operators agree under it better than under
  the gene rotation (1.22 vs 1.48) — although it is 0.58 from `Q_load`.
  `Q_load` is the gene-loading frame map; operators may constrain the
  low-variance directions differently from the loadings. The registered
  rule still says "fails", and stands.
- Using the flow (odd order) resolves the axis-sign ambiguity that defeated
  curvature-only estimators (O1, O1f).
- Remaining error at n ≥ 2000 (~0.18 from the truth even when started
  there) is attributable to field-fit noise. Targeted follow-up: less noisy
  fields (more drift seeds per copy) and clean start seeds (O1g).

## Gated numbers (as reported)

- `e5_outcome`: fails
- `err_F_E5`: 0.581694
- `err_R_E5_1000`: 0.327357
- `err_R_E5_2000`: 0.188794
- `err_R_E5_2759`: 0.178775
- `err_R_E5_500`: 0.726983
- `err_random`: 0.995846
- `latent_check_max`: 0
- `min_n_e5`: 0
- `n_starts_near_best_F`: 1
- `o1_repro_max_abs`: 0
- `symmetry_max`: 2.64196e-07
