# O1k — Review

Reviewer: planner (Claude), 2026-11-02. Branch `task/O1k` @ 18d4c73 (report data from 7ac85ab), Mac worktree `~/Documents/scopalign-manus-O1k`.

## Checker
`python checker/check_report.py tasks/O1k.md reports/O1k/report.json --prev <O1g report>` → **PASS** (0 fail, 0 warn).

## Recomputation
Light recompute from `cases.csv`: every gated median recomputed exactly, as were `frac_hat_below_star_B` (16/27) and `fe_u_3000` / `fe_J_3000`. Not recomputed: fits, fields, LBFGS.

| gated | value |
|---|---|
| n_cases / start_overlap | 72 / 0 |
| loss_copy_check | 0.0 |
| symmetry_max | 2.18e-7 |
| negctrl_err_shape_A | 0.925 |
| errA_max | 1.002 |
| errB n3000 γ 1.15 / 1.33 / 1.6 | 0.081 / 0.041 / 0.026 |
| errB n10000 γ 1.15 / 1.33 / 1.6 | 0.449 / 0.018 / 0.632 |
| frac_hat_below_star_B | 0.593 |
| fe_u_3000 / fe_J_3000 | 0.266 / 0.300 |

Deviation (undeclared, minor): O1g's cache has no Test R fits, because O1g reused O1b's fits. The copy check therefore read O1b's fits (`scopalign-manus-O1b/reports/O1b/cache/drifts/R_n2759_r0_*`, sha256 in `check.json`) and imported O1g's script by path. The two losses agree exactly (L(I) 2.3117, L(Q*) 0.1451). This does not affect validity.

## Decision rule outcome
Validity: 72 cases, overlap 0, copy check 0, symmetry 2.2e-7, negative control 0.93 ≥ 0.30 → valid.
Rule 1: errA_max 1.00 > 0.10 → **test_fault**.

## Verdict
**ACCEPT** — outcome **test_fault**. With exact fields, L(Q*) ≈ 1e-30, yet the best of 20 random starts stops at L 1–13 (err 0.72–1.00): the search fails, not the information.

Reading (not the outcome): with fitted fields, the same random starts recover the rotation at n = 3000 for every γ > 1 (0.08 / 0.04 / 0.03). At n = 10000 they miss in 2 of 3 replicates at γ 1.15 and 1.6, and in every miss L(Q̂) > L(Q*). So the misses are search misses, not noise moving the minimum. Fitted fields are about 27–30 % off the exact ones and still carry the rotation. The random-start result is therefore unreliable in both directions, which is what O1l tests.

## Quarantined (not evidence)
- Planner seeds: fitted fields seem to give a smoother loss than exact fields (exact Gaussian-mixture Hessians are sharp far from the wells). d = 10, Gaussian wells and a pure rotation, so nothing here transfers to PBMC at d = 30 without a test.
