# T02b review

reviewed_commit: 4fd6380 (task/T02b)
date: 2026-10-15
verdict: ACCEPT (stopped run) — outcome: **inconclusive (stop)**, rule 0 (`n_groups` 5 < 6)
recomputation: light — group count recomputed from `fine_labels.csv` and
`k2_cells.csv` (5 groups, same names); `latent_maxdiff` 0.0 and
`k2_repro_max_abs` 0.0 as reported. Nothing after the stop was run.

## Checks

- Checker PASS (status `stopped`); branch changes only `reports/T02b/`.
- Deviation: fits rerun once after a file-format fix (same inputs, config,
  seeds) — cannot move a gated number; accepted.

## Groups (≥ 40 in both)

CD16+ NK cells (53/53), Classical monocytes (268/272), MAIT cells (63/70),
Naive B cells (108/100), Tcm/Naive helper T cells (249/221).
Near misses: Tcm/Naive cytotoxic T cells (39/42), Non-classical monocytes
(41/33), Memory B cells (33/35).

## Verdict

**ACCEPT.** Inconclusive: the 1k-cell runs give 5 fine groups at ≥ 40
cells, one short of the floor.

## Planner reading (not findings)

Correspondence count is limited by cell numbers (~970 per dataset), not by
the label model. Lowering the floor now would be chosen knowing the counts.
A design that does not need many groups: take `Q_load` (fixed by the
shared gene loadings, no operator fitting) as the rotation and test whether
the operator frames agree under it.

## Gated numbers (as reported)

- `k2_repro_max_abs`: 0
- `latent_maxdiff`: 0
- `n_cells_a`: 969
- `n_cells_b`: 969
- `n_groups`: 5
- `outcome`: stopped
- `symmetry_max`: 2.30366e-07
