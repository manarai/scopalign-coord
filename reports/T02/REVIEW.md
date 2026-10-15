# T02 review

reviewed_commit: b947a8c (task/T02)
date: 2026-10-15
verdict: ACCEPT (stopped run) — outcome: **inconclusive (stop)**, rule 0
recomputation: **not recomputed (Tom's instruction, 2026-10-15)**

## Checks

- Checker → PASS (0 fail, 0 warn), status `stopped`, `stop_rule` quotes rule 0.
- Branch changes only `reports/T02/`.
- Script read against the task: 20 starts (identity + `ortho_group`
  seeds 20261022…20261040), `Q₀·matrix_exp(S − Sᵀ)`, LBFGS as specified,
  1.01 near-best rule, null generator 20261022, gene subset before
  normalisation. `k2_repro_max_abs` 0.0 (K2's eigenvalues reproduced).
- Stopped before the views stage, so no UMAP figures (allowed for a stop).

## Reported numbers

| gated | value |
|---|---|
| `frame_ratio` | 1.24972 |
| `heldout_ratio` | 2.52649 |
| `heldout_resid` | 0.6752 |
| `k2_repro_max_abs` | 0 |
| `loss_best` | 0.111546 |
| `n_cell_types_tested` | 4 |
| `n_cells_a` | 969 |
| `n_cells_b` | 969 |
| `n_starts_near_best` | 1 |
| `outcome` | stopped |
| `q_op_vs_pos` | 0.953905 |
| `resid_aligned` | 0.333985 |
| `resid_identity` | 0.732471 |
| `resid_null` | 0.267248 |
| `resid_pos` | 0.793384 |
| `symmetry_max` | 2.30366e-07 |

## Decision rule outcome

0. Validity passes (`k2_repro_max_abs` 0.0, symmetry 2.3e-07), but
   `n_starts_near_best` = 1 < 2 → **inconclusive (stop)**: only start 6
   reached the best loss (0.1115); the next best were 0.1192 (identity),
   0.1206, 0.1217 — 6.9–9.1% higher. The optimum was not found
   reproducibly.

## Verdict

**ACCEPT** the stopped run. Outcome inconclusive.

## Planner reading (not findings; numbers as reported)

- Alignment clearly helps: residual 0.73 with no rotation → 0.33 with the
  best one, against a fit-noise null of 0.27 (`frame_ratio` 1.25 would have
  been within 1.5).
- But the fitted rotation generalises poorly to an unseen cell type:
  held-out residual 0.68 (`heldout_ratio` 2.53), close to no alignment.
  With 4 type-mean operators and a 435-parameter rotation, the fit can
  partly match noise.
- The loss surface has several near-equal minima (identity is only 7%
  worse than the best), and the operator-based rotation is far from the
  position-based one (`q_op_vs_pos` 0.95). Both point to the same limit:
  four types do not pin down a 30-d rotation.

## Quarantined (not evidence)

- None.
