# O1n — Low-dimension sweep: does the operator rotation work in 2–10 PCs?

status: pre-registered (approved by Tom 2026-11-04)
blocked_by: none (O1m ACCEPT)
seed: 20261014
planner_commit: 9b616b0
kind: follow-up (O1m: fails at d = 30, even with exact fields; O1l: works at d = 10 on synthetic data)
motivated_by: Tom 2026-11-04 — "it seem to work with 2D PC, can you instruct manus to test?" The scJDO2 orientation session found the fitted energy tracks density (ρ vs −log KDE 0.84–0.92) in 2-D but not in 3-D or above on aging tissue.

## Purpose

Repeat O1's PBMC Test R (known rotation) and Test F (real halves) on the
first `d` PCs only, **d ∈ (2, 3, 5, 10)**. This asks whether
operator-estimated rotation works where score operators are still well
estimated.

One control is essential. In low d, each copy's own principal axes already
fix the rotation up to axis signs. So a **positions-only baseline** that
picks each sign from the skew of the cells along that axis must be beaten,
or the operators add nothing.

## Scope

### In scope
- **Data:** K0's PCA-30 latents `z_H1`, `z_H2` (as O1), truncated to
  their first `d` columns. sha256 recorded; read-only.
- **Test R** (n = 2759, r ∈ (0, 1, 2)): O1's construction (subsets from
  `default_rng(20261030 + r)`), with
  `Q*_{d,r} = ortho_group.rvs(dim=d, random_state=20261031 + r)` and
  `Y = z_H1[S_Y, :d] · Q*ᵀ`.
- **Test F:** `X = z_H1[:, :d]`, `Y = z_H2[:, :d]`. The halves' PCA is
  re-derived as O1 (`latent_check_rel`). The target is
  `Q_load,d = U Vᵀ` of `svd(W_2[:d] W_1[:d]ᵀ)`.
- **Fits:** O1's drift config with input dimension `d`, 5 seeds per copy
  (20261014, 20261016–19), ensemble mean. 4 d × (3 Test R cases + Test F)
  × 2 copies × 5 seeds = **160 fits**, cached and resumable.
  `symmetry_check` per fit.
- **Estimator E5s** (O1l exactly): O1b's E5 loss (`D + H`; `H` subsets
  `min(256, n)` from `default_rng(20261040 + r)`, Test F `+ 0`).
  - Starts `Q₀(s) = U_Y diag(s) U_Xᵀ` (each copy's own `numpy.cov`,
    `eigh` descending) over **all** `2^d` sign patterns.
  - Screen by `L`; LBFGS (O1b settings) from the `min(8, 2^d)` lowest;
    `Q̂` = lowest final `L`. Start workers ≤ 3.
- **Baseline (positions only, no operators):** `s_j = sign(m3_X,j) ·
  sign(m3_Y,j)`, where `m3_·,j` is the third central moment of that copy's
  cells along its axis `j`. `Q̂_skew = U_Y diag(s) U_Xᵀ`, with no
  refinement.
- **Errors:** `err = ‖Q̂ − Q_target‖_F / √(2d)`, for `Q̂` and `Q̂_skew`.
  Also record `L(Q̂)`, `L(Q*)` (Test R) and `screen_rank_best`.
- **Test F functional:** O1g's `mix` (30 nearest other-copy cells, ÷ 0.5),
  computed in the d-dim space for no rotation, `Q̂_skew`, `Q̂` and
  `Q_load,d`. `sil` (coarse types) is reported, not gated.
- **Fit quality (reading):** for d ∈ (2, 3), Test R r = 0, copy X, each
  seed: Spearman ρ between the fitted energy at the cells and
  `−log scipy.stats.gaussian_kde(cells.T)(cells.T)` (default bandwidth).
  The gated value is the median over seeds.
- **Views:** Test F, d = 2, as scatter plots (no UMAP is needed in 2-D) for
  no rotation / skew / E5s / `Q_load`, coloured by copy and by coarse
  type, with `mix` and `sil` in the titles.

### Out of scope
- Other d, n, seeds, epochs or starts; changes to the loss or LBFGS.
- Labels, pairing or gene loadings in any estimate.
- The aging data, and correcting the dimensions beyond `d`.

## Pre-registration

### Gated numbers
For each `d` ∈ (2, 3, 5, 10), written `_d2`, `_d3`, `_d5`, `_d10`:
- `err_R_d*`: Test R, E5s, median over r.
- `err_R_skew_d*`: Test R, skew baseline, median over r.
- `mix_F_d*`: Test F mix, E5s.
- `mix_F_skew_d*`: Test F mix, skew baseline.
- `mix_F_noA_d*`: Test F mix, no rotation.

| name | definition | computed by |
|------|-----------|-------------|
| `symmetry_max` | max `median_asym` over the 160 fits | `python reports/O1n/scripts/o1n_run.py --stage fits` |
| `latent_check_rel` | max relative difference, re-derived vs K0 latents | `--stage check` |
| `err_R_d2` | see above | `--stage testR` |
| `err_R_d3` | | same |
| `err_R_d5` | | same |
| `err_R_d10` | | same |
| `err_R_skew_d2` | | same |
| `err_R_skew_d3` | | same |
| `err_R_skew_d5` | | same |
| `err_R_skew_d10` | | same |
| `mix_F_d2` | | `--stage testF` |
| `mix_F_d3` | | same |
| `mix_F_d5` | | same |
| `mix_F_d10` | | same |
| `mix_F_skew_d2` | | same |
| `mix_F_skew_d3` | | same |
| `mix_F_skew_d5` | | same |
| `mix_F_skew_d10` | | same |
| `mix_F_noA_d2` | | same |
| `mix_F_noA_d3` | | same |
| `mix_F_noA_d5` | | same |
| `mix_F_noA_d10` | | same |
| `rho_kde_d2` | fit-quality ρ, d = 2 | `--stage testR` |
| `rho_kde_d3` | fit-quality ρ, d = 3 | same |
| `d_max_op` | largest d meeting rule 1 (0 if none) | `--stage report` |
| `o1n_outcome` | `"op_works_lowd"`, `"positions_suffice"` or `"fails"` | same |

### Decision rule
0. **Validity → REVISE:** `symmetry_max > 1e-4` or `latent_check_rel > 1e-4`.
1. For each d, the operator estimate **passes** if
   - `err_R_d ≤ 0.10` (Test R),
   - `mix_F_d ≥ 0.80` (Test F), and
   - it beats positions: `err_R_d ≤ err_R_skew_d − 0.05` **or**
     `mix_F_d ≥ mix_F_skew_d + 0.10`.

   `d_max_op` = the largest d that passes.
2. **op_works_lowd** if some d passes.
3. **positions_suffice** if no d passes, but at some d E5s meets both
   thresholds (Test R and Test F) without beating the baseline.
4. Otherwise **fails**.

Reading (not the outcome):
- `rho_kde_d2` and `rho_kde_d3` show whether the fits track density, as in
  the orientation session.
- `L(Q̂) < L(Q*)` with a large err means a displaced minimum, as in O1m.
- With op_works_lowd, a low-d operator rotation can only align the top d
  PCs. Whether that helps a full correction is a separate test.

### What would change my mind
If operators cannot beat a positions-only sign rule even in 2–3 PCs,
where the fits track density, then the operators carry no orientation
information beyond what cell positions give. Operator-estimated rotation
should be dropped, and scOPAlign repositioned around validation and
anchoring.

## Procedure

1. `--stage check` — inputs sha256, `latent_check_rel`, `Q*` and `Q_load`
   per d; `inputs.json`.
2. `--stage fits` — 160 fits; `fit_symmetry.csv`.
3. `--stage testR` — 12 cases; `testR.csv` (d, r, err, err_skew, L_hat,
   L_star, screen_rank_best); `rho_kde.csv`; `Q_hats.npz`.
4. `--stage testF` — 4 cases; `testF.csv` (d, mix and sil for no rotation /
   skew / E5s / Q_load, err, err_skew).
5. `--stage report` — `fig_sweep.png` (err_R and mix_F vs d, E5s vs skew,
   lines at 0.10 and 0.80), `fig_views_d2.png`, `report.json`, checksums,
   lockfile.

`--stage all`; detached, stdin `/dev/null`, resumable, checkpoint per case.
Mac (**Manus Desktop, not the cloud sandbox**), worktree
`~/Documents/scopalign-manus-O1n`, branch `task/O1n` from `main`; commit
and push after each stage.

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `rho_kde.csv`
- `fit_symmetry.csv`
- `fig_sweep.png`
- `fig_views_d2.png`
- `scripts/o1n_run.py`

`Q_hats.npz` and `inputs.json` go in `artifacts`.

## Checker expectations

- `gated` keys: exactly the twenty-six above (or a subset with `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
