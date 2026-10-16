# O1 review — pre-finalisation note

date: 2026-10-15

The O1 run completed, but the checker rejected the 12 `cases/case_*.npz`
entries listed under `deliverables` ("deliverable not in task"). This is a
planner spec conflict: O1 asks for the case arrays to be recorded in
`report.json` with sha256, but they are not deliverables.

Binding clarification (no change to any computation or gated number):
record the 12 case arrays under the new optional `artifacts` list
(`[{"path": "cases/case_<n>_<r>.npz", "sha256": ...}]`), added to the
schema and checker on main; keep `deliverables` exactly as the task lists.
Merge `origin/main` first to get the updated checker.

---

# O1 review

reviewed_commit: a6cd92e (task/O1)
date: 2026-10-16
verdict: ACCEPT — selection: **none** (PROJECT.md rule 14 → O1b onward)
recomputation: light — Test R medians recomputed from `testR.csv` (all 20
match); `err_F` for E0–E4 recomputed from `testF.npz` (`Q_hat_E*` vs
`Q_load`, all exact).

## Checks

Checker PASS (0 fail, 0 warn); 12 case arrays recorded under `artifacts`;
two declared deviations (resume of saved cases; detached restart) cannot
move a gated number. `latent_check_max` 0.0; symmetry 2.6e-07.

## Results

All estimators are far from the 0.10 recovery line, at every n: Test R
errors 0.63–1.02 (random rotation: 0.996); Test F errors 0.69–0.92 (best
E3 0.694). **No estimator works → `selected_estimator` "none".**

## Views (planner-made, `planner_views.png`; illustration)

Test F, separate PCAs of the K0 halves:

| view | mix | sil |
|---|---|---|
| no rotation | 0.312 | 0.385 |
| E3 (best O1 estimate) | 0.278 | 0.374 |
| `Q_load` (known answer) | 0.973 | 0.396 |

## Planner diagnostic (not findings) — information present, search failed

From the saved case arrays: the seed-averaged, group-mean operator mismatch
evaluated **at the true rotation `Q*`** is far below its value at the
identity or at random rotations:

| n | L(Q*) | L(I) | L(random), median of 20 |
|---|---|---|---|
| 500 | 0.39–0.56 | 1.01–1.20 | 1.04–1.24 |
| 1000 | 0.14–0.15 | 0.73–0.76 | 0.77–0.80 |
| 2000 | 0.06–0.08 | 0.68–0.75 | 0.69–0.77 |
| 2759 | 0.05–0.06 | 0.68–0.72 | 0.70–0.72 |

So the operators **do** encode the rotation (≈ 13× lower mismatch at the
truth for n ≥ 2000); the estimators did not find it. Twenty random starts
in O(30), single-seed operators (E0), and invariant-spectrum matching
(E1–E4) did not reach the right basin. This motivates a targeted follow-up
(O1f): seed-averaged group operators with a far wider global search.

## Gated numbers (as reported)

- `err_F_E0`: 0.919882
- `err_F_E1`: 0.91727
- `err_F_E2`: 0.868702
- `err_F_E3`: 0.693718
- `err_F_E4`: 0.773814
- `err_R_E0_1000`: 0.980894
- `err_R_E0_2000`: 0.693958
- `err_R_E0_2759`: 0.872722
- `err_R_E0_500`: 0.982208
- `err_R_E1_1000`: 0.900905
- `err_R_E1_2000`: 0.766066
- `err_R_E1_2759`: 0.632662
- `err_R_E1_500`: 0.987065
- `err_R_E2_1000`: 0.936454
- `err_R_E2_2000`: 0.911206
- `err_R_E2_2759`: 0.849372
- `err_R_E2_500`: 1.02102
- `err_R_E3_1000`: 0.860723
- `err_R_E3_2000`: 0.815259
- `err_R_E3_2759`: 0.702923
- `err_R_E3_500`: 0.987437
- `err_R_E4_1000`: 0.962134
- `err_R_E4_2000`: 0.871734
- `err_R_E4_2759`: 0.813474
- `err_R_E4_500`: 0.984677
- `err_random`: 0.995846
- `latent_check_max`: 0
- `min_n_selected`: 0
- `selected_estimator`: none
- `symmetry_max`: 2.64196e-07
