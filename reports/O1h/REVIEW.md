# O1h review

reviewed_commit: fc5bf89 (task/O1h)
date: 2026-10-18
verdict: ACCEPT — reading: **partial** (as registered)
recomputation: light — all top-k errors recomputed from `Q_hats.npz` with
`Q*` regenerated (`ortho_group`, seed 20261031; match exact); `best` = argmin L.

## Checks

Checker PASS. `repro_err_abs` 0, `repro_L_abs` 0 (O1g's case reproduces
exactly). `axes_sorted` false: a swap in the low-variance tail only; the
variance order of the top 12 axes is 0…11, so the top-k numbers are
unaffected.

## Results (E5b, n 2759, r 0)

| k | variance share | err_k(Q̂) | identity | random median |
|---|---|---|---|---|
| 1 | | 0.013 | 0.98 | 1.01 |
| 3 | 0.72 | **0.053** | 0.93 | 1.00 |
| 5 | 0.80 | **0.125** | 0.92 | 1.00 |
| 10 | 0.89 | 0.213 | 0.97 | 1.00 |
| 30 | 1.00 | 0.604 | 1.01 | 1.00 |

`n_starts_top5_ok` 0; `disp_hat` 0.251 vs random 1.00. Only the best start
(L 0.146; next 0.213) is accurate on the top axes.

## Reading

Rule: `err_top5` 0.125 lies in (0.10, 0.25] → `partial`. The operators
recover the rotation on the directions carrying most of the variance (top 3
to 0.05) and the error grows smoothly with k — information fades with
variance; it is not a search failure. Cells move by 25 % of their spread
from the truth, against 100 % for a random rotation.

## Next (planner judgement, shown to Tom)

A restricted-subspace estimator: rotation estimated on the top axes only
(identity / no change on the rest) is supported by this result. Pre-register
before running.
