# K0c review

reviewed_commit: e01d6f3 (task/K0c)
date: 2026-10-15
verdict: ACCEPT — decision-rule outcome: `outcome_avg = "fails"` (reported to Tom beside K0b)

## Checker

`check_report.py tasks/K0c.md reports/K0c/report.json --prev reports/K0b/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/K0c/`. The script is
K0b's with only the layer (FA-30) and paths changed, as the task requires
(the copied `started` timestamp is cosmetic).

## Recomputation

Planner's own script (K0b's recompute with FA-30 inputs): ten drift fits,
per-seed and averaged restricted spectra, K0 statistics on K0's stored
draw indices, refit draws. Averaged eigenvalues agree to ≤ 8.9e-6
(FA-30 eigenvalues are ~30× larger than PCA-30's); every gated number
agrees to ≤ 1e-6 relative.

| gated | report | recomputed | match |
|---|---|---|---|
| `basin_ratio_avg` | 1.95532 | 1.95532 | yes |
| `basin_sep_h1_avg` | 6.14062 | 6.14062 | yes |
| `basin_sep_h2_avg` | 5.924 | 5.924 | yes |
| `class_tv_ratio_avg` | 2.63395 | 2.63395 | yes |
| `fitnoise_basin_h1` | 1.79109 | 1.79109 | yes |
| `fitnoise_basin_h2` | 1.79224 | 1.79224 | yes |
| `fitnoise_spectral_h1` | 3.15788 | 3.15788 | yes |
| `fitnoise_spectral_h2` | 2.56568 | 2.56568 | yes |
| `k0_repro_max_rel` | 5.27681e-06 | 5.27681e-06 | yes |
| `k0_vs_refit_basin` | 0.778362 | 0.778362 | yes |
| `k0_vs_refit_spectral` | 1.06515 | 1.06515 | yes |
| `outcome_avg` | fails | fails | yes |
| `refit_w2_h1` | 3.98901 | 3.98901 | yes |
| `refit_w2_h2` | 2.82436 | 2.82436 | yes |
| `spectral_ratio_avg` | 4.39076 | 4.39076 | yes |
| `spectral_sep_h1_avg` | 14.3415 | 14.3415 | yes |
| `spectral_sep_h2_avg` | 16.2123 | 16.2123 | yes |
| `symmetry_max` | 1.24278e-07 | 1.24278e-07 | yes |

Figures `fig_fitnoise.png` and `fig_spectra_avg.png`: as specified.

## Decision rule outcome

0. `symmetry_max` 1.2e-07, `k0_repro_max_rel` 5.3e-06 → valid.
1. Part B: separations ≥ 5.9; `spectral_ratio_avg` 4.39 > 3.0 → **fails**.
   Reported to Tom beside K0b; no automatic layer switch.
2. Part A: `k0_vs_refit_spectral` 1.07, `k0_vs_refit_basin` 0.78 (both
   ≤ 1.5) → "H1 and H2 differ no more than two fits of one half".

## Deviations

None listed. None found.

## Verdict

**ACCEPT.**

## Planner reading (not findings)

| | PCA-30 (K0b) | FA-30 (K0c) |
|---|---|---|
| halves vs refit noise, spectral / basin | 0.85 / 0.98 | 1.07 / 0.78 |
| fit noise ÷ within-fit null, spectral (H1, H2) | 2.15, 1.87 | 3.16, 2.57 |
| averaged-operator ratios spectral / basin / class | 1.81 / 1.49 / 1.47 | 4.39 / 1.96 / 2.63 |
| effect of 5-seed averaging on pooled W2 | 0.0304 → 0.0291 | 3.53 → **5.10** |

- Both layers' halves sit within their own refit noise, so neither shows a
  real split-half difference. They differ in how noisy and stable the fit
  is: FA-30's fit-to-fit noise is larger, its refit distances are broad and
  lumpy (some seed pairs disagree strongly), and averaging seeds makes its
  halves *more* different, which suggests seed-dependent, distinct fit
  solutions. PCA-30's fits are tighter and averaging is close to neutral.
- This supports keeping PCA-30 (Tom's design choice, PROJECT.md rule 1).

## Quarantined (not evidence)

- Observations: input sha256 values (6; confirmed by recomputation inputs).
