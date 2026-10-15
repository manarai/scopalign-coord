# K0c — Fit-noise control and seed-averaged operators for FA-30 (K0b design)

status: pre-registered (approved by Tom 2026-10-15)
blocked_by: none (K0, K0b ACCEPT)
seed: 20261014
planner_commit: 8c769af
kind: follow-up
motivated_by: reports/K0b/REVIEW.md; Tom: give FA-30 the same fit-noise check as PCA-30

## Purpose

FA-30 failed K0 narrowly (spectral 3.019 > 3.0; basin 1.427, class TV
2.395). K0b showed for PCA-30 that K0's within-fit null leaves out drift-fit
noise. The same reason applies to FA-30, so this task is **K0b's design
unchanged, with FA-30 in place of PCA-30**: K0's FA-30 latents and draw
indices, the same five seeds per half, the same refit draws and the same
rule.

- **A. Measures fit-to-fit noise** for FA-30.
- **B. Repeats K0's split-half test on seed-averaged FA-30 operators**,
  with K0's rule unchanged.

Part B is FA-30's outcome; Part A explains it. K0's FA-30 outcome
("fails") stays as recorded. The layer in use stays PCA-30 (Tom's design
choice, PROJECT.md rule 1) unless Tom revisits it on this result.

## Scope

### In scope
- **Inputs from K0, as committed (record sha256 of each):**
  `data/pbmc_10k_v3/k0_halves.csv` (the split), `reports/K0/genes.txt`,
  `reports/K0/latent_fa30_h1.npy`, `reports/K0/latent_fa30_h2.npy`,
  `reports/K0/draw_indices.npz`, `reports/K0/report.json`. No refit of FA.
  Cell order and labels: as in K0 (each half in file order of
  `pbmc_10k_v3.h5ad` QC cells; `cell_type` from that file).
- **Tangent bases:** recomputed exactly as in K0 (k = 40 other cells,
  `np.cov`, `np.linalg.eigh`, last 10 columns).
- **Drift fits:** for each half and each seed `s` in
  `SEEDS = (20261014, 20261016, 20261017, 20261018, 20261019)`:
  `EnergyDriftAdapter(dim=30, config=EnergyDriftConfig(hidden=96, depth=3,
  n_epochs=1500, seed=s, batch_size=256, loss="esm", lr=5e-4, beta=1.0,
  esm_laplacian="hutchinson", device="cpu")).fit(z)`;
  `J_s = jacobian_autodiff(z)` (float64 from here on);
  `symmetry_check(z, n_cells=100)`. Ten fits in total. Seed 20261014
  repeats K0's fit.
- **Per-seed spectra:** `λ_s = eigvalsh(Vᵀ·½(J_s + J_sᵀ)·V)` (10 per cell);
  `basin_index_s = −max λ_s`; class by `lead_eig` with
  `tol = 0.1 × median |λ_s|` over that half — as in K0.
- **Averaged operator:** `J̄ = mean_s ½(J_s + J_sᵀ)` per cell;
  `λ̄ = eigvalsh(Vᵀ J̄ V)`; `basin_index`, class and `tol` from `λ̄` as above.
- **Part A, refit draws:** `rng = numpy.random.default_rng(20261020)`.
  For half H1, then H2, 100 draws each, draw `d = 0…99`: seed pair
  `(a, b) = (20261014, SEEDS[1 + d % 4])`; for each tested type in sorted
  order, `rng.choice` of `2·m_t` of that type's cells (without replacement,
  half's file order); the first `m_t` cells are read under fit `a`, the rest
  under fit `b`; per type spectral W2 (on `λ_a`, `λ_b`) and basin W1; mean
  over types. Tested types and `m_t` exactly as in K0 (stored in
  `draw_indices.npz`).
- **K0 statistics on seed-20261014 fits and on averaged operators:** using
  K0's stored draw indices (`draw_indices.npz`; no new draws), recompute
  K0's pooled, null and separation statistics (a) on the seed-20261014
  spectra (to reproduce K0) and (b) on `λ̄`. Definitions, means, medians and
  ratio convention exactly as K0.

### Out of scope
- PCA-30, FA-10, any other layer, `d`, `k`, `r`, epochs, seeds or
  averaging scheme (no weighting, no outlier removal).
- New split, new gene list, new FA fit, new draws for the K0 statistics.
- Any alignment; any other dataset.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `symmetry_max` | max `median_asym` over the 10 fits | `python reports/K0c/scripts/k0c_run.py --stage fits` |
| `k0_repro_max_rel` | max relative difference between K0's reported `spectral_ratio_fa30`, `basin_ratio_fa30`, `class_tv_ratio_fa30` and the same recomputed from the seed-20261014 fits | `--stage compare` |
| `refit_w2_h1` | median spectral W2 over the 100 refit draws, H1 | same |
| `refit_w2_h2` | same, H2 | same |
| `fitnoise_spectral_h1` | `refit_w2_h1` ÷ median spectral W2 over K0's null draws on H1 (seed-20261014 spectra) | same |
| `fitnoise_spectral_h2` | same, H2 | same |
| `fitnoise_basin_h1` | median basin W1 over refit draws on H1 ÷ median basin W1 over K0's H1 null draws | same |
| `fitnoise_basin_h2` | same, H2 | same |
| `k0_vs_refit_spectral` | median spectral W2 over K0's pooled draws (seed 20261014) ÷ median spectral W2 over all 200 refit draws | same |
| `k0_vs_refit_basin` | same for basin W1 | same |
| `spectral_ratio_avg` | K0's spectral ratio on `λ̄` | same |
| `basin_ratio_avg` | K0's basin ratio on `λ̄` | same |
| `class_tv_ratio_avg` | K0's class TV ratio on `λ̄` | same |
| `spectral_sep_h1_avg` | K0's H1 spectral separation on `λ̄` | same |
| `spectral_sep_h2_avg` | same, H2 | same |
| `basin_sep_h1_avg` | K0's H1 basin separation on `λ̄` | same |
| `basin_sep_h2_avg` | same, H2 | same |
| `outcome_avg` | K0's per-config rule applied to the `_avg` numbers: `"survives"`, `"fails"` or `"inconclusive"` | same |

### Decision rule
Top to bottom; the first matching line decides.
0. **Validity:** `symmetry_max > 1e-4` → **REVISE**.
   `k0_repro_max_rel > 1e-4` → **REVISE** (seed-20261014 fits do not
   reproduce K0).
1. **Part B (decides).** Apply K0's per-config rule to the `_avg` numbers:
   any `*_sep_*_avg` < 2 → inconclusive; `spectral_ratio_avg > 3.0` or
   `basin_ratio_avg > 3.0` → fails; all three `_avg` ratios ≤ 1.5 →
   survives; otherwise inconclusive. Record as `outcome_avg`.
   - Whatever the outcome, the planner reports it to Tom beside K0b's
     PCA-30 result; Tom decides whether the layer in use (PCA-30) stays.
     No automatic switch.
2. **Part A (explains; not decisive).** Reported beside the outcome:
   `k0_vs_refit_spectral ≤ 1.5` and `k0_vs_refit_basin ≤ 1.5` reads as
   "H1 and H2 differ no more than two fits of one half"; > 1.5 on either
   reads as "the halves differ beyond fit noise".

### What would change my mind
If FA-30's two halves differ clearly more than two fits of one half
(`k0_vs_refit_spectral` or `k0_vs_refit_basin` > 1.5), FA-30's K0 failure
is not drift-fit noise, and FA-30 is less reproducible than PCA-30.

## Procedure

1. `k0c_run.py --stage fits` — tangent bases; ten drift fits (kept
   locally); per-seed and averaged spectra → `ops_avg_h1.npz`,
   `ops_avg_h2.npz` (arrays: `lam_avg` (n,10), `lam_seeds` (5,n,10),
   `basin_index_avg`, `cls_avg`, `cell_type`, `barcode`), all float64.
2. `k0c_run.py --stage compare` — refit draws (`refit_draws.csv`: half,
   draw, seed_a, seed_b, spectral_w2, basin_w1), K0 statistics on seed
   20261014 and on `λ̄` (`k0_stats.csv`: spectra ∈ {seed20261014, avg},
   kind, draw, spectral_w2, basin_w1, class_tv), gated numbers,
   `fig_fitnoise.png` (per half: histograms of null, refit and K0 pooled
   spectral W2 on one axis, medians marked) and `fig_spectra_avg.png`
   (per tested type: averaged restricted-eigenvalue density, H1 vs H2).
3. `k0c_run.py --stage report` — `report.json`, checksums, lockfile.

`k0c_run.py` must accept `--stage all`. A validity stop is reported as
`"status": "stopped"` (MANUS.md 8b).

## Deliverables

- `report.json`
- `requirements.lock`
- `ops_avg_h1.npz`
- `ops_avg_h2.npz`
- `refit_draws.csv`
- `k0_stats.csv`
- `fig_fitnoise.png`
- `fig_spectra_avg.png`
- `scripts/k0c_run.py`

## Checker expectations

- `gated` keys: exactly the eighteen above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
