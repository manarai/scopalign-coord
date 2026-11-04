# O1m — Review

Reviewer: planner (Claude), 2026-11-04. Branch `task/O1m` @ df0dc41, Mac worktree `~/Documents/scopalign-manus-O1m`.

## Checker
`python checker/check_report.py tasks/O1m.md reports/O1m/report.json --prev reports/O1l/report.json` → **PASS** (0 fail, 0 warn).

## Recomputation
Light recompute: every gated value recomputed exactly from `synth_cases.csv`, `testR.csv`, `testF.csv` and `basin.csv`. Not recomputed: fits, screen, LBFGS.

| gated | value |
|---|---|
| symmetry_max / fits_match / l_star_repro_max | 2.1e-7 / 1 / 0 |
| negctrl_err_shape_d30 | 0.989 |
| errA_d30 γ 1.05 / 1.15 | 0.116 / 0.813 |
| errB_d30 γ 1.05 / 1.15 | 0.591 / 0.273 |
| err_R_2759 (top5 / top10) | 0.654 (0.145 / 0.256) |
| frac_hat_below_star_R | 1/3 |
| basin_radius_R | 0 |
| mix_F / mix_F_noA | 0.751 / 0.312 |
| sil_F / sil_F_noA | 0.390 / 0.381 |
| frame_ratio_F | 1.41 |

Deviations (declared, 4): an ExactEnergy buffer fix before any case was written; a ranking bookkeeping fix (same ordering); `final_L` added to basin rows (needed by the procedure); screen evaluated in batches of 2. None changes a value or a selection.

## Decision rule outcome
Validity: all pass.
Rule 1: Test R fails (0.65 > 0.10).
Rule 2: Test F fails (mix 0.75 < 0.80).
Rule 3: errA_d30_g115 0.81 > 0.10 → **search_d30**.

## Verdict
**ACCEPT** — outcome **search_d30**. At d = 30 the axis-sign search does not recover the rotation even from exact fields at γ 1.15 (L 2.8–3.4 vs L* ≈ 1e-30). PBMC fails Test R and Test F.

Reading (not the outcome):
- **The PBMC minimum is displaced.** Every basin start, even at t = 0.05, converges to one point per case: err 0.149 / 0.161 / 0.176, with L 0.095–0.100 below L(Q*) 0.141–0.146. Even a perfect search would land about 0.16 from the truth. So PBMC at d = 30 has two problems: the search, and fitting noise that moves the minimum.
- **Fitted γ 1.05 at d = 30:** L(Q̂) < L(Q*) in all 3 cases, at err 0.39–0.65, so the fitted loss is flat or shifted. Exact γ 1.05 is almost recovered (0.08–0.14).
- Test F mix 0.75 is below O1g's random-start E5b (0.80).

## Quarantined (not evidence)
- Planner seed: this agrees with the scJDO2 orientation session's diagnosis. Density-curvature (score) operators lose shallow structure above a few dimensions, which motivates a low-d test (2–10 PCs).
