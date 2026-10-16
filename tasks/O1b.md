# O1b — Operator-estimated rotation, E5: field equivariance (no correspondences)

status: pre-registered (written 2026-10-15 before O1's result; issue only if O1 selects "none", PROJECT.md rule 14)
blocked_by: O1 review (selected_estimator = "none")
seed: 20261014
planner_commit: 2255ea5
kind: planned (operator-central alignment, PROJECT.md rules 13–14)

## Purpose

Estimate the rotation **from the operator fields themselves**, without
pairing cells. The fitted drift of each copy is a smooth function that can
be evaluated anywhere. If copy Y is copy X seen in a rotated frame, then at
every X cell the Y field evaluated at the rotated position equals the
rotated X field:

- flow: `u^Y(Q x) = Q u^X(x)`;
- operator: `J^Y(Q x) = Q J^X(x) Qᵀ`;

and symmetrically for Y cells mapped back with `Qᵀ`. E5 finds the `Q` that
best satisfies these, using every cell. The estimate is then applied to the
latent space (`z ↦ Q z`) in O2, per PROJECT.md rule 13. Same known-answer
tests, cell counts, replicates and pass rule as O1.

## Scope

### In scope
- **Data, Test R and Test F constructions, drift config, seeds:** exactly
  `tasks/O1.md` (subsets from `default_rng(20261030 + r)`, `Q*` from
  `ortho_group.rvs(30, random_state=20261031 + r)`, n ∈ (500, 1000, 2000,
  2759), r ∈ (0, 1, 2); Test F = K0 halves with re-derived PCA, target
  `Q_load`). Five drift fits per copy (seeds 20261014, 20261016–19),
  refit here (deterministic; must reproduce O1's operators if O1 saved
  them — see `o1_repro_max_abs`).
- **Fields as torch functions:** for each fitted adapter `d`, the energy
  network `d._model` (a SiLU MLP, `E: ℝ³⁰ → ℝ`), cast to float64. For a
  point `x`: `u_s(x) = −β ∇E_s(x)`, `J_s(x) = −β ∇²E_s(x)` (β = 1),
  via `torch.func` (`grad`, `hessian`, `vmap`) or `torch.autograd` with
  `create_graph=True`. Ensemble fields: `ū(x) = mean_s u_s(x)`,
  `J̄(x) = mean_s J_s(x)` over the five seeds of that copy.
- **E5 loss** (copies X and Y; cells `x_i` of X, `y_j` of Y):
  ```
  D(Q) = [ Σ_i ‖ū^Y(Q x_i) − Q ū^X(x_i)‖² / Σ_i ‖ū^X(x_i)‖²
         + Σ_j ‖ū^X(Qᵀ y_j) − Qᵀ ū^Y(y_j)‖² / Σ_j ‖ū^Y(y_j)‖² ] / 2
  H(Q) = [ Σ_{i∈A} ‖J̄^Y(Q x_i) − Q J̄^X(x_i) Qᵀ‖_F² / Σ_{i∈A} ‖J̄^X(x_i)‖_F²
         + Σ_{j∈B} ‖J̄^X(Qᵀ y_j) − Qᵀ J̄^Y(y_j) Q‖_F² / Σ_{j∈B} ‖J̄^Y(y_j)‖_F² ] / 2
  L(Q) = D(Q) + H(Q)
  ```
  `D` uses all cells. `H` uses fixed subsets `A`, `B` of
  `min(256, n)` cells each, drawn once per case with
  `numpy.random.default_rng(20261040 + r)` (Test F: `+ 0`),
  `choice` without replacement, X first then Y. The fields of the
  copy being evaluated at a moved point are recomputed at every optimiser
  evaluation (they depend on `Q`); the fields at the copy's own cells are
  computed once.
- **Optimisation:** `Q = Q₀ · matrix_exp(S − Sᵀ)`, torch float64 CPU,
  `LBFGS(lr=1, max_iter=500, tolerance_grad=1e-10, tolerance_change=1e-12,
  line_search_fn="strong_wolfe")`, one `.step` per start. **20 starts:**
  `Q₀ = I`, then `ortho_group.rvs(30, random_state=20261021 + i)`,
  `i = 1…19`. `Q̂` = lowest final `L`. Network weights are frozen.
- **Error measures, reference:** as O1 (`err = ‖Q̂ − Q_target‖_F / sqrt(60)`;
  `err_R` = median over the 3 replicates; `err_random` as O1).

### Out of scope
- Pairing cells, labels, gene loadings or positions in the estimate
  (`Q_load` is the Test F target only).
- Other loss weights, subset sizes, start counts, seeds, cell counts or
  thresholds; retraining or modifying the drift networks.
- Any real-pair application (that is O2).

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_E5_500` | E5, Test R, n = 500 (median of 3) | `python reports/O1b/scripts/o1b_run.py --stage testR` |
| `err_R_E5_1000` | n = 1000 | same |
| `err_R_E5_2000` | n = 2000 | same |
| `err_R_E5_2759` | n = 2759 | same |
| `err_F_E5` | E5, Test F | `--stage testF` |
| `n_starts_near_best_F` | Test F starts with final `L` ≤ 1.01 × best | same |
| `err_random` | as O1 | `--stage testR` |
| `symmetry_max` | max `median_asym` over all fits | same |
| `latent_check_max` | as O1 | `--stage testF` |
| `o1_repro_max_abs` | max \|seed-averaged operator here − O1's saved `M̄`\| over cases O1 saved; `0` if O1 saved none | same |
| `e5_outcome` | `"works"` or `"fails"` by the rule below | `--stage report` |
| `min_n_e5` | smallest n with `err_R_E5_<n> ≤ 0.10`; 0 if none | same |

### Decision rule
0. **Validity:** `symmetry_max > 1e-4`, `latent_check_max > 1e-4` or
   `o1_repro_max_abs > 1e-4` → **REVISE**.
1. E5 **works** if `err_F_E5 ≤ 0.25` **and** some `err_R_E5_<n> ≤ 0.10` at
   n ≤ 2759; `min_n_e5` = smallest such n. Otherwise **fails** → O1c next
   (PROJECT.md rule 14).
2. If E5 works with `min_n_e5 ≤ 1000`, O2 may use the 1k v2/v3 pair;
   otherwise O2 needs ≥ `min_n_e5` cells per dataset.

### What would change my mind
If E5 cannot recover a known rotation to 0.10 even at 2759 cells, the
fitted operator fields are not equivariant enough under refitting for a
correspondence-free estimate; the remaining options rely on summaries
(O1c) or joint fitting (O1d).

## Procedure

1. `o1b_run.py --stage testR` — 120 drift fits (kept locally), E5 on the 12
   cases; `testR.csv` (n, replicate, err, best_start, final_L, D, H).
2. `o1b_run.py --stage testF` — 10 fits, E5; `testF.csv`, `Q_hat_F.npy`.
3. `o1b_run.py --stage report` — `fig_recovery.png` (err_R vs n with O1's
   five estimators if available, lines at 0.10 and `err_random`),
   `report.json`, checksums, lockfile.

`o1b_run.py` must accept `--stage all`.

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `fig_recovery.png`
- `scripts/o1b_run.py`

## Checker expectations

- `gated` keys: exactly the twelve above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
