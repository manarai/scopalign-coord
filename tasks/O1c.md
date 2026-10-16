# O1c — Operator-estimated rotation, E7: joint diagonalisation + commutator refinement

status: pre-registered (written 2026-10-15 before any O1-series result; may run in parallel, PROJECT.md rule 14)
blocked_by: none
seed: 20261014
planner_commit: 85b10e7
kind: planned (operator-central alignment, PROJECT.md rules 13–14)

## Purpose

Estimate the rotation from a **common eigenbasis** of each copy's operators.
If both copies' wells share axes, jointly diagonalising a set of operators
gives each copy's axes; matching axes across copies gives the rotation.
Hessians alone cannot tell an axis from its mirror image, so axis signs are
set from the **skew of the drift field** along each axis (a reflection flips
its sign). A final refinement minimises the **commutator** between rotated
and target operators, which compares directions and ignores well depth.
Same known-answer tests, cell counts, replicates and pass rule as O1.

## Scope

### In scope
- **Data, Test R / Test F constructions, drift config, five seeds,
  seed-averaged fields `M̄`, `ū`, matching groups (coarse types with ≥ 20
  cells in both copies), `Q` parametrisation and LBFGS settings, error
  measures, `err_random`:** exactly `tasks/O1.md`.
- **Operator set per copy:** for each matching group `t` (sorted), rank its
  cells by seed-averaged energy `Ē = mean_s E_s(x)` (ascending, stable) and
  split into 10 equal-size deciles (`numpy.array_split` of the ranked
  list); `N_{t,d}` = mean `M̄` over the cells of decile `d`. Set
  `𝒩 = {N_{t,d}}` (groups × 10), identically indexed in X and Y.
- **Joint diagonalisation:** orthogonal `V` (30×30) minimising
  `off(V) = Σ_{t,d} ‖offdiag(Vᵀ N_{t,d} V)‖_F² / Σ_{t,d} ‖N_{t,d}‖_F²`,
  `V = V₀ · matrix_exp(S − Sᵀ)`, LBFGS as O1, 20 starts (`I`, then
  `ortho_group.rvs(30, random_state=20261021 + i)`, `i = 1…19`), lowest
  `off`. Done separately for X (`V_X`) and Y (`V_Y`).
- **Axis matching:** profile of axis `k` = vector of
  `diag(Vᵀ N_{t,d} V)_k` over all `(t,d)`; cost
  `C_kl = ‖p_k^X − p_l^Y‖²`; `scipy.optimize.linear_sum_assignment` → `π`.
- **Axis signs:** for each axis `k`, skew
  `γ_k = mean_i (v_kᵀ ū_i)³ / (mean_i (v_kᵀ ū_i)²)^{3/2}` over all cells of
  the matching groups in that copy; `s_k = sign(γ_k^X · γ_{π(k)}^Y)`
  (`+1` if either is 0).
- **Initial rotation:** `Q₀ = V_Y[:, π] · diag(s) · V_Xᵀ`.
- **Commutator refinement:** minimise
  `K(Q) = Σ_{t,d} ‖[Q N^X_{t,d} Qᵀ, N^Y_{t,d}]‖_F² / Σ_{t,d} ‖N^X_{t,d}‖_F² ‖N^Y_{t,d}‖_F²`
  (`[A, B] = AB − BA`) from `Q₀` with LBFGS (one start), then add the
  drift term of O1's E2 evaluated on the decile-mean drift vectors
  `ū_{t,d}` (`D(Q) = Σ ‖Q ū^X_{t,d} − ū^Y_{t,d}‖² / Σ ‖ū^Y_{t,d}‖²`) and
  minimise `K + D` from the result (one start). `Q̂` = final.
- **Views (illustration; not gated):** as O1b (Test F: no rotation, `Q̂`,
  `Q_load`; UMAP, `mix`, `sil`, frame ratio).

### Out of scope
- Cell pairing, labels beyond the matching groups, gene loadings or
  positions in the estimate (`Q_load` is the Test F target only).
- Other decile counts, start counts, thresholds or seeds.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_E7_500` | E7, Test R, n = 500 (median of 3) | `python reports/O1c/scripts/o1c_run.py --stage testR` |
| `err_R_E7_1000` | n = 1000 | same |
| `err_R_E7_2000` | n = 2000 | same |
| `err_R_E7_2759` | n = 2759 | same |
| `err_F_E7` | E7, Test F | `--stage testF` |
| `err_R_E7init_1000` | error of `Q₀` (before refinement), n = 1000 | `--stage testR` |
| `offdiag_F` | final `off(V_X)` and `off(V_Y)` mean, Test F | `--stage testF` |
| `err_random` | as O1 | `--stage testR` |
| `symmetry_max` | as O1 | same |
| `latent_check_max` | as O1 | `--stage testF` |
| `e7_outcome` | `"works"` or `"fails"` | `--stage report` |
| `min_n_e7` | smallest n with `err_R_E7_<n> ≤ 0.10`; 0 if none | same |

### Decision rule
0. Validity: `symmetry_max > 1e-4` or `latent_check_max > 1e-4` → REVISE.
1. E7 **works** if `err_F_E7 ≤ 0.25` and some `err_R_E7_<n> ≤ 0.10`;
   `min_n_e7` = smallest such n. Otherwise **fails**.
Selection across the O1 series follows PROJECT.md rule 14.

### What would change my mind
If the operators' common eigenbasis does not transfer across refits
(`err_R_E7_init` near `err_random`), cluster-level operators do not share
stable axes, and summary-based estimators cannot work.

## Procedure

1. `o1c_run.py --stage testR` — 120 fits (kept locally), E7 on the 12
   cases; `testR.csv` (n, replicate, err_init, err, off_X, off_Y).
2. `o1c_run.py --stage testF` — 10 fits, E7; `testF.csv`, `Q_hat_F.npy`;
   views → `views.csv`, `fig_umap.png`.
3. `o1c_run.py --stage report` — `fig_recovery.png`, `report.json`,
   checksums, lockfile.

`o1c_run.py` must accept `--stage all`.

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `views.csv`
- `fig_umap.png`
- `fig_recovery.png`
- `scripts/o1c_run.py`

## Checker expectations

- `gated` keys: exactly the twelve above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
