# O1f review

reviewed_commit: f8c5dfb (task/O1f)
date: 2026-10-16
verdict: ACCEPT — outcome: E9 **fails**
recomputation: light — Test R medians from `testR.csv`, `err_F` from
`Q_hat_F.npy` vs `Q_load`, `frac_reached_truth` and `n_alt_minimum` from
`testR.csv`: all exact.

## Checks

Checker PASS; branch changes only `reports/O1f/`. The run was restarted once
after a multiprocessing serialisation crash (no cases had been processed);
the fix is an implementation detail allowed by the task (parallel workers)
and cannot move a gated number. No deviations listed.

## Results

`err_R_E9`: 0.98 (500), 0.70 (1000), 1.36 (2000), 1.36 (2759);
`err_F_E9` 0.53. In all 12 cases the search reached the truth's loss level
(`frac_reached_truth` 1.0) **and** found a lower minimum away from the
truth (`n_alt_minimum` 12) → **E9 fails**; reading: the group-mean
objective has lower minima away from the truth.

Views (Test F): no rotation mix 0.312 / sil 0.380; E9 0.347 / 0.377;
`Q_load` 0.973 / 0.394.

## Planner diagnostic (not findings) — sign ambiguity and residual mixing

From O1f's per-case `Q_hat` and O1's `Q_star`:
- At n ≥ 2000, `Q̂ᵀQ*` is nearly diagonal (mean |diag| ≈ 0.86–0.88 in
  4 of 6 cases): `Q̂ ≈ Q* · diag(±1)` — the estimate has the right axes but
  undetermined **axis signs**. In 4 of 6 cases it is close to `−Q*`
  (err 0.35–0.38).
- Cause: `Q M Qᵀ` is unchanged when `Q` is multiplied by a reflection along
  the operators' own eigenaxes (including `−I`). Hessian (curvature)
  operators alone cannot fix axis directions. Only odd-order information —
  the drift (flow) field, which reverses under a reflection — can.
- With oracle signs the error is still 0.35–0.38 (n ≥ 2000): residual mixing
  between near-degenerate axes. Drift-based sign rules from summary
  statistics agree with the true signs only 40–83% of the time.
- Consequence for the design: an operator-only estimate must use the full
  fitted dynamics (flow + curvature) at the cell level, not curvature
  summaries. O1b (E5, field equivariance) is that estimator.

## Gated numbers (as reported)

- `e9_outcome`: fails
- `err_F_E9`: 0.526842
- `err_R_E9_1000`: 0.70232
- `err_R_E9_2000`: 1.36349
- `err_R_E9_2759`: 1.36286
- `err_R_E9_500`: 0.981497
- `frac_reached_truth`: 1
- `min_n_e9`: 0
- `n_alt_minimum`: 12
