# O1m — Transfer the fixed search to d = 30: synthetic bridge, then PBMC

status: pre-registered (approved by Tom 2026-11-02)
blocked_by: none (O1k, O1l ACCEPT)
seed: 20261014
planner_commit: 49a6784
kind: follow-up (O1l: from axis starts, fitted operators recover a known rotation at d = 10)
motivated_by: Tom 2026-11-02 — "Yes, draft it" (PBMC transfer test). O1l passed on fitted fields at d = 10. O1j (unreviewed) used PBMC position axes but only set the complement's signs and never refined the full 30-d rotation. O1g's PBMC estimates sit at L(Q̂) ≈ L(Q*) with err 0.60, which points to a flat loss rather than a missed minimum.

## Purpose

Three things separate O1l's success from the PBMC failures: **dimension**
(10 vs 30), **search** (O1g used random starts) and **data** (real copies may
differ by more than a rotation). O1m changes one at a time:

- **Part B (bridge):** O1k's landscape built at d = 30, exact and fitted
  fields, with the E5s search adapted to d = 30. This shows whether
  dimension alone breaks recovery.
- **Part P (PBMC):** the same search on O1g's Test R (known Q*) and Test F
  (functional), with O1g's fields unchanged.
- **Basin (PBMC Test R):** how far from Q* the search still returns, and
  whether the loss is flat there.

## Scope

### In scope
- **Search E5s-30 (both parts, per case):**
  1. Axes: `U_X`, `U_Y` from each copy's own `numpy.cov` (rows = cells),
     `eigh`, descending (rule 13a).
  2. Complement signs: with axes 1–10 at `s = +1`, set `s_j` for
     `j = 11…30` by greedy flow sweeps exactly as O1j step 2 (full `D`, all
     cells; max 5 sweeps).
  3. Screen: all 1024 sign patterns of axes 1–10 (binary order as O1l),
     complement fixed from step 2, `Q₀ = U_Y diag(s) U_Xᵀ`, ranked by full
     `L = D + H`.
  4. Refine: full 30-d LBFGS (O1g's E5b `Problem`, `Closure` and settings)
     from the 8 lowest (ties broken by index); `Q̂` = lowest final `L`.
     Start workers ≤ 3.
- **Part B — synthetic d = 30:** O1k's rule with `d = 30`:
  - Means: `μ_k` in coordinates 1–3, from `default_rng(20264200)` as O1k.
  - Shapes: `R_k = ortho_group.rvs(30, random_state=20264210 + k)`;
    `σ_j² = γ^(15.5 − j)`, j = 1…30.
  - Rotations: `Q*_r = ortho_group.rvs(30, random_state=20264100 + r)`.
  - Cells: `default_rng(20264300 + r)`, O1k's draw order.
  - Cases: **γ ∈ (1.0, 1.05, 1.15)** (condition numbers 1, 4.1, 58);
    n = 3000; r ∈ (0, 1, 2).
  - Rung A (exact fields): all 9 cases; γ = 1.0 is the negative control.
  - Rung B (fitted fields): γ ∈ (1.05, 1.15) only. Fits use O1g's `cfg`
    with input dimension 30 and 5 seeds per copy, giving 6 cases × 2 copies
    × 5 seeds = **60 fits** (cached, resumable).
  - Errors as O1k with `√(2d)`; `err_shape` on coordinates 4–30.
- **Part P — PBMC:**
  - Test R: n = 2759, r ∈ (0, 1, 2), constructed as O1g. Test F as O1g.
  - Fields: O1g's 20-seed ensembles. The fits are read-only, and each path
    and sha256 must match O1j's `fits_used.csv` (O1b and O1g caches).
  - Errors: `err` (= `err_30`) and O1h's `err_5` and `err_10`, plus `L(Q̂)`
    and `L(Q*)`.
  - Test F: `mix`, `sil` and `frame` as O1g. UMAP views in 4 rows: no
    rotation, O1g's E5b, E5s-30, and the gene-loading `Q_load`.
- **Basin (Part P Test R, 3 cases):** O1l's procedure at d = 30.
  - Directions: `A_i` from `default_rng(20264400 + i).standard_normal((30, 30))`,
    i = 1…5.
  - Distances: t ∈ (0.05, 0.10, 0.20, 0.40).
  - Record each start's final `err` and `L`. A start returns if its final
    err ≤ 0.10.

### Out of scope
- New PBMC fits, other γ, n, d, seeds or starts.
- Changes to the loss or LBFGS.
- Labels, pairing or gene loadings in the estimate.
- Any start using `Q*` outside the basin diagnostic.
- O1j's review (separate).

## Pre-registration

### Gated numbers
Medians are over r.

| name | definition | computed by |
|------|-----------|-------------|
| `symmetry_max` | max `median_asym` over all fits loaded or made | `python reports/O1m/scripts/o1m_run.py --stage check` |
| `fits_match` | 1 if all PBMC fit sha256 match O1j's `fits_used.csv`, else 0 | same |
| `l_star_repro_max` | max \|L(Q*) − O1g's `testR.csv` L_star\| over the 3 n = 2759 cases | `--stage pbmc` |
| `negctrl_err_shape_d30` | B, rung A, γ = 1.0: median `err_shape` | `--stage synth` |
| `errA_d30_g105` | B, rung A, γ = 1.05: median `err` | same |
| `errA_d30_g115` | B, rung A, γ = 1.15 | same |
| `errB_d30_g105` | B, rung B, γ = 1.05 | same |
| `errB_d30_g115` | B, rung B, γ = 1.15 | same |
| `err_R_2759` | P, Test R: median `err` | `--stage pbmc` |
| `err_R_top5_2759` | median `err_5` | same |
| `err_R_top10_2759` | median `err_10` | same |
| `frac_hat_below_star_R` | fraction of the 3 cases with `L(Q̂) ≤ L(Q*)` | same |
| `basin_radius_R` | median basin radius over the 3 Test R cases | `--stage basin` |
| `mix_F` | Test F mixing, E5s-30 | `--stage pbmc` |
| `mix_F_noA` | Test F mixing, no rotation | same |
| `sil_F` | Test F silhouette, E5s-30 | same |
| `sil_F_noA` | Test F silhouette, no rotation | same |
| `frame_ratio_F` | Test F frame ratio, E5s-30 | same |
| `o1m_outcome` | `"works"`, `"functional_only"`, `"search_d30"`, `"dimension_limited"` or `"data_limited"` | `--stage report` |

### Decision rule
0. **Validity → REVISE:** `symmetry_max > 1e-4`, `fits_match = 0`,
   `l_star_repro_max > 1e-8`, or `negctrl_err_shape_d30 < 0.30`.
1. **works:** Test R passes (`err_R_2759 ≤ 0.10`) and Test F passes
   (`mix_F ≥ 0.80`, `frame_ratio_F ≤ 1.5`, `sil_F ≥ 0.9 × sil_F_noA`).
   E5s-30 then becomes the O1-series estimator for O2 (rule 14, inserted
   after O1j).
2. **functional_only:** Test F passes, Test R does not.
3. Otherwise PBMC fails, and Part B says why:
   - **search_d30** if `errA_d30_g115 > 0.10` (the search fails at d = 30
     even with exact fields);
   - **dimension_limited** if `errB_d30_g115 > 0.10` (fitted fields lose
     the rotation at d = 30);
   - **data_limited** otherwise (synthetic d = 30 works, PBMC does not:
     the copies differ by more than a rotation, or PBMC's loss is flat).

Reading (not the outcome):
- `frac_hat_below_star_R` near 1 together with a large `err_R` means a
  flat or shifted PBMC minimum, not a search miss.
- `basin_radius_R` shows how close a start must be on real fields.
- `err_R_top5` and `err_R_top10` vs O1h and O1j show how far orientation
  reaches beyond the top axes.
- γ = 1.05 vs 1.15 brackets the anisotropy needed at d = 30.

### What would change my mind
If synthetic d = 30 works but PBMC fails, with `L(Q̂) ≤ L(Q*)`, then
operator-estimated rotation is limited by the data and not the method. The
operator-only plan would then rest on anchoring plus per-group directions,
not one global rotation.

## Procedure

1. `--stage check` — fit sha256 and symmetry; `inputs.json`.
2. `--stage fits` — 60 synthetic fits; `fit_symmetry.csv`.
3. `--stage synth` — 15 cases; `synth_cases.csv`.
4. `--stage pbmc` — Test R and Test F; `testR.csv`, `testF.csv`,
   `screen.csv`, `Q_hats.npz`, `Q_hat_F.npy`.
5. `--stage basin` — 3 × 5 × 4 = 60 refinements; `basin.csv`.
6. `--stage report` — `fig_ladder.png`, `fig_basin.png`, `fig_umap.png`,
   `report.json`, checksums, lockfile.

`--stage all`; detached, stdin `/dev/null`, resumable, checkpoint per case.
Mac worktree `~/Documents/scopalign-manus-O1m`, branch `task/O1m` from
`main`; commit and push after each stage. Caches of O1b, O1g, O1j and O1k
are read-only.

## Deliverables

- `report.json`
- `requirements.lock`
- `synth_cases.csv`
- `testR.csv`
- `testF.csv`
- `basin.csv`
- `fit_symmetry.csv`
- `fig_ladder.png`
- `fig_basin.png`
- `fig_umap.png`
- `scripts/o1m_run.py`

`screen.csv`, `Q_hats.npz`, `Q_hat_F.npy` and `inputs.json` go in `artifacts`.

## Checker expectations

- `gated` keys: exactly the nineteen above (or a subset with `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
