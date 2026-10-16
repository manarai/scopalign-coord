# O1f — Targeted fix: seed-averaged group operators with a global rotation search (E9)

status: pre-registered (2026-10-16, after O1 and its planner diagnostic; before any O1b–O1e result)
blocked_by: none (O1 ACCEPT)
seed: 20261014
planner_commit: dd787af
kind: follow-up (Tom's standing instruction: when a method is close, test a targeted fix)
motivated_by: reports/O1/REVIEW.md — loss at the true rotation is ~13× below random (n ≥ 2000), so the information is present and the search failed

## Purpose

O1's estimators missed the rotation although the seed-averaged, group-mean
operators clearly encode it. E9 keeps that objective and replaces the
search: thousands of short starts, cluster-eigenbasis candidates, full
refinement of the best. It reuses O1's saved operator arrays (no drift
fits). Same known-answer tests, thresholds and selection logic as O1.
Operator-only estimate (PROJECT.md rule 13).

## Scope

### In scope
- **Inputs (record sha256):** O1's case arrays `reports/O1/cases/case_<n>_<r>.npz`
  (12; fields `M_X`, `M_Y` (5, n, 30, 30), `u_X`, `u_Y`, `type_X`, `type_Y`,
  `Q_star`) and `reports/O1/testF.npz` (`Mbar_X`, `Mbar_Y`, `Q_load`), plus
  for Test F group labels `data/pbmc_10k_v3/pbmc_10k_v3.h5ad` and
  `data/pbmc_10k_v3/k0_halves.csv` (K0 order and coarse `cell_type`).
- **Operators:** per cell `M̄ = ½(A + Aᵀ)` with `A = mean over the 5 seeds`
  (cases) or `Mbar_*` as stored (Test F), float64. Groups: coarse types with
  ≥ 20 cells in both copies (sorted). `G_t^X`, `G_t^Y` = group means of `M̄`.
  Weights `w_t ∝ min(n_X,t, n_Y,t)`.
- **Objective:** `L(Q) = Σ_t w_t ‖Q G_t^X Qᵀ − G_t^Y‖_F² / Σ_t w_t ‖G_t^Y‖_F²`.
- **Parametrisation / optimiser:** `Q = Q₀ · matrix_exp(S − Sᵀ)`, torch
  float64 CPU, `S = 0` at start, `LBFGS(lr=1, line_search_fn="strong_wolfe",
  tolerance_grad=1e-10, tolerance_change=1e-12)` with `max_iter` as below.
- **Search:**
  1. **Random stage:** 2000 starts `Q₀ = ortho_group.rvs(30,
     random_state=20262000 + i)`, `i = 0…1999`, each `max_iter = 50`; keep
     the 50 lowest final `L`.
  2. **Eigenbasis candidates:** for each group `t`: `V^X, V^Y` = eigenvectors
     of `G_t^X`, `G_t^Y` (`eigh`, ascending); `Q₀ = V^Y V^Xᵀ`; then greedy
     sign fixing: for `k = 0…29`, flip the sign of column `k` of `V^X` if it
     lowers `L(V^Y V^Xᵀ)`; repeat the sweep until no flip helps (max 10
     sweeps). One candidate per group.
  3. **Refinement:** LBFGS with `max_iter = 500` from each of the 50 + one
     candidate per group; `Q̂` = lowest final `L`.
- **Diagnostics per case:** `L(Q̂)`, `L(Q_target)` (`Q*` or `Q_load`),
  `err`, and whether `L(Q̂) ≤ 1.01 × L(Q_target)` (search reached the
  truth's level).
- Cases may be processed in parallel worker processes (results identical).
- **Views (illustration; not gated):** as O1b (Test F: no rotation, `Q̂`,
  `Q_load`; UMAP, `mix`, `sil`).

### Out of scope
- New drift fits, other operator summaries, other start counts or
  thresholds; labels beyond the groups, gene loadings or positions in the
  estimate (`Q_load` is the Test F target only).

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_E9_500` | median over 3 replicates, n = 500 | `python reports/O1f/scripts/o1f_run.py --stage testR` |
| `err_R_E9_1000` | n = 1000 | same |
| `err_R_E9_2000` | n = 2000 | same |
| `err_R_E9_2759` | n = 2759 | same |
| `err_F_E9` | Test F | `--stage testF` |
| `frac_reached_truth` | fraction of the 12 Test R cases with `L(Q̂) ≤ 1.01 L(Q*)` | `--stage testR` |
| `n_alt_minimum` | Test R cases with `L(Q̂) < 0.99 L(Q*)` **and** `err > 0.10` (a different, lower minimum) | same |
| `e9_outcome` | `"works"` or `"fails"` | `--stage report` |
| `min_n_e9` | smallest n with `err_R_E9_<n> ≤ 0.10`; 0 if none | same |

### Decision rule
1. E9 **works** if `err_F_E9 ≤ 0.25` and some `err_R_E9_<n> ≤ 0.10`;
   `min_n_e9` = smallest such n. Otherwise **fails**.
2. Reading (not the outcome): `n_alt_minimum > 0` means the group-mean
   objective has lower minima away from the truth — more constraints
   (more operator groups) are needed, not a better search.
Selection across the O1 series follows PROJECT.md rule 14 (O1f placed
right after O1 in the fixed order, set before any O1b–O1e result).

### What would change my mind
If the wide search reaches the truth's loss level but the error stays above
0.10, the objective has a flat valley around the truth; if it finds lower
minima elsewhere, group-mean operators do not identify the rotation.

## Procedure

1. `o1f_run.py --stage testR` — `testR.csv` (n, replicate, err, L_hat,
   L_star, best_source ∈ {random, eigen_<group>}).
2. `o1f_run.py --stage testF` — `testF.csv`, `Q_hat_F.npy`; views →
   `views.csv`, `fig_umap.png`.
3. `o1f_run.py --stage report` — `fig_recovery.png` (err_R vs n with O1's
   estimators), `report.json`, checksums, lockfile.

`o1f_run.py` must accept `--stage all`. Run detached (MANUS.md "Long runs").

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `views.csv`
- `fig_umap.png`
- `fig_recovery.png`
- `scripts/o1f_run.py`

## Checker expectations

- `gated` keys: exactly the nine above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
