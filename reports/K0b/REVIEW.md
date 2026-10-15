# K0b review

reviewed_commit: e3ba7e9 (task/K0b)
date: 2026-10-15
verdict: ACCEPT — decision-rule outcome: `outcome_avg = "inconclusive"` (back to Tom)

## Checker

`check_report.py tasks/K0b.md reports/K0b/report.json --prev reports/K0/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/K0b/`. Input sha256
values for the K0 split, genes, latents, draw indices and report are
recorded.

## Recomputation

Planner's own script: tangent bases from K0's latents, ten drift fits (five
seeds × two halves), per-seed and averaged restricted spectra, K0's
statistics on K0's stored draw indices (seed 20261014 and averaged), the
refit draws on `default_rng(20261020)`. Averaged eigenvalues agree with the
deliverables to ≤ 9.4e-8; every gated number agrees to ≤ 1e-6 relative
(values below are the recomputed ones; report identical to shown precision).

| gated | report | recomputed | match |
|---|---|---|---|
| `basin_ratio_avg` | 1.48809 | 1.48809 | yes |
| `basin_sep_h1_avg` | 6.61314 | 6.61314 | yes |
| `basin_sep_h2_avg` | 5.94479 | 5.94479 | yes |
| `class_tv_ratio_avg` | 1.47497 | 1.47497 | yes |
| `fitnoise_basin_h1` | 2.25529 | 2.25529 | yes |
| `fitnoise_basin_h2` | 1.50136 | 1.50136 | yes |
| `fitnoise_spectral_h1` | 2.15 | 2.15 | yes |
| `fitnoise_spectral_h2` | 1.86577 | 1.86577 | yes |
| `k0_repro_max_rel` | 9.7017e-07 | 9.7017e-07 | yes |
| `k0_vs_refit_basin` | 0.979622 | 0.979622 | yes |
| `k0_vs_refit_spectral` | 0.847724 | 0.847724 | yes |
| `outcome_avg` | inconclusive | inconclusive | yes |
| `refit_w2_h1` | 0.0359505 | 0.0359505 | yes |
| `refit_w2_h2` | 0.0351126 | 0.0351126 | yes |
| `spectral_ratio_avg` | 1.8075 | 1.8075 | yes |
| `spectral_sep_h1_avg` | 27.1127 | 27.1127 | yes |
| `spectral_sep_h2_avg` | 24.6951 | 24.6951 | yes |
| `symmetry_max` | 2.00214e-07 | 2.00214e-07 | yes |

Figures: `fig_fitnoise.png` (per half: null, refit and K0-pooled spectral
W2 with medians) and `fig_spectra_avg.png` (averaged restricted spectra per
tested type, H1 vs H2) are as specified.

## Decision rule outcome

0. `symmetry_max` 2.0e-07 ≤ 1e-4; `k0_repro_max_rel` 9.7e-07 ≤ 1e-4 → valid.
1. Part B on averaged operators: separations 5.9–27.1 (≥ 2); spectral 1.807
   (> 1.5, ≤ 3), basin 1.488 (≤ 1.5), class TV 1.475 (≤ 1.5) →
   **inconclusive** → back to Tom.
2. Part A (explains): `k0_vs_refit_spectral` 0.848 and
   `k0_vs_refit_basin` 0.980 (both ≤ 1.5) → **"H1 and H2 differ no more
   than two fits of one half."**

## Deviations

None listed. None found.

## Verdict

**ACCEPT.** Outcome inconclusive; K1 not written until Tom decides.

## Planner reading (not findings)

- **The yardstick, not PCA-30, is the problem.** Two drift fits of the
  *same* half, differing only in seed, are ~1.5–2.3× further apart than
  K0's within-fit null (`fitnoise_*`). The H1-vs-H2 distance is *smaller*
  than that (0.85 spectral, 0.98 basin). So any two independently fitted
  operator sets, even on identical cells, would miss K0's "≤ 1.5 × null"
  line. K0's rule asks for more agreement than the estimator can give
  itself.
- **Averaging five seeds barely moves anything** (pooled spectral W2
  0.0304 → 0.0291, null 0.0177 → 0.0161). If fit noise were independent
  across seeds, it would have dropped by ~√5. So the fit-to-fit spread is
  largely shared across seeds, and seed averaging is not the remedy.
- **Basin and class now pass** on averaged operators (1.49, 1.47); only the
  spectral ratio (1.81) stays above 1.5. Eigenvalue scale and class
  counts agree between halves.
- A refit-calibrated criterion (compare to two fits of the same data) is
  the mechanistically right yardstick, but it was not pre-registered before
  K0b ran, so it cannot change K0 or K0b. It can govern K1/K2, decided
  before they run.

## Quarantined (not evidence)

- Observations: input sha256 values (confirmed).
