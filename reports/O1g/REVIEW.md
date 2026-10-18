# O1g review

reviewed_commit: ea199bc (task/O1g; deliverables committed by the planner from the finished run, unchanged)
date: 2026-10-18
verdict: ACCEPT — outcome: E5b **fails** (Test R fails; Test F passes functionally)
recomputation: light — Test R medians and `frac_hat_below_star` from
`testR.csv`, Test F rule from `testF.csv` (all exact).

## Checks

Checker PASS. One declared deviation: start workers limited to three
(memory); starts, seeds and optimiser settings unchanged — cannot move a
gated number. `start_overlap` 0, `symmetry_max` 2.2e-7, `latent_check_max` 0.

## Results

| n | err per replicate | err_R_E5b (median) |
|---|---|---|
| 2000 | 0.577, 0.732, 0.663 | 0.663 |
| 2759 | 0.604, 0.668, 0.594 | 0.604 |

`frac_hat_below_star` 0.5. Five cases have `L(Q̂)` within ~6 % of `L(Q*)`;
case 2759/r2 has `L(Q̂)` 0.219 vs `L(Q*)` 0.141 (search missed the truth's
basin there).

Test F (functional, from O1g on): `mix_F` 0.802 ≥ 0.80, `frame_ratio_F`
1.24 ≤ 1.5, `sil_F` 0.387 ≥ 0.9 × 0.381 → **passes** (mixing only just).
`err_F` 0.630 (reported, not in the rule).

No `err_R ≤ 0.10` → **E5b fails**; `min_n_e5b` 0.

## Reading

Quadrupling the field ensemble left Test R where O1b's clean-start cases
were (~0.6). Per "what would change my mind", the leftover error is not
drift-fit noise: the E5 objective is nearly flat over rotations far from
the truth. The estimate still aligns cells (Test F), consistent with the
operators fixing the high-variance directions only — O1h tests this on one
Test R case (top-k axis error). Views: no rotation mix 0.31 / frame 3.62;
E5b 0.80 / 1.24; known answer 0.97 / 1.48.

## Next

O1h (issued). Its reading decides between a restricted-subspace estimator
and O1c (rule 14).
