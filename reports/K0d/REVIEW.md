# K0d review

reviewed_commit: bc2a6ab (task/K0d)
date: 2026-10-15
verdict: ACCEPT — decision-rule outcome: `outcome_avg = "fails"` (reported to Tom beside K0b, K0c)

## Checker

`check_report.py tasks/K0d.md reports/K0d/report.json --prev reports/K0c/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/K0d/`. The script is
K0c's with the layer set to FA-10, `d = 10` and the tangent basis removed,
as the task requires. Cosmetic: `fig_spectra_avg.png`'s axis says
"restricted" although FA-10 uses full 10-d spectra (as specified).

## Recomputation

Planner's own script (K0b/K0c recompute with `d = 10`, full spectra, FA-10
inputs): ten drift fits, averaged spectra, K0 statistics on K0's stored
draw indices, refit draws. Averaged eigenvalues and every gated number are
**identical** to the report.

| gated | report | recomputed | match |
|---|---|---|---|
| `basin_ratio_avg` | 1.33017 | 1.33017 | yes |
| `basin_sep_h1_avg` | 8.11543 | 8.11543 | yes |
| `basin_sep_h2_avg` | 6.53974 | 6.53974 | yes |
| `class_tv_ratio_avg` | 1.6645 | 1.6645 | yes |
| `fitnoise_basin_h1` | 2.60867 | 2.60867 | yes |
| `fitnoise_basin_h2` | 3.07086 | 3.07086 | yes |
| `fitnoise_spectral_h1` | 7.59103 | 7.59103 | yes |
| `fitnoise_spectral_h2` | 7.91489 | 7.91489 | yes |
| `k0_repro_max_rel` | 4.2538e-06 | 4.2538e-06 | yes |
| `k0_vs_refit_basin` | 0.530178 | 0.530178 | yes |
| `k0_vs_refit_spectral` | 0.9274 | 0.9274 | yes |
| `outcome_avg` | fails | fails | yes |
| `refit_w2_h1` | 15.8233 | 15.8233 | yes |
| `refit_w2_h2` | 18.1074 | 18.1074 | yes |
| `spectral_ratio_avg` | 4.08945 | 4.08945 | yes |
| `spectral_sep_h1_avg` | 30.9924 | 30.9924 | yes |
| `spectral_sep_h2_avg` | 36.0283 | 36.0283 | yes |
| `symmetry_max` | 1.50185e-07 | 1.50185e-07 | yes |

## Decision rule outcome

0. `symmetry_max` 1.5e-07, `k0_repro_max_rel` 4.3e-06 → valid.
1. Part B: separations ≥ 6.5; `spectral_ratio_avg` 4.09 > 3.0 → **fails**.
2. Part A: `k0_vs_refit_spectral` 0.93, `k0_vs_refit_basin` 0.53 → halves
   within refit noise.

## Deviations

None listed. None found.

## Verdict

**ACCEPT.**

## Planner reading (not findings) — the three layers side by side

| | PCA-30 (K0b) | FA-30 (K0c) | FA-10 (K0d) |
|---|---|---|---|
| K0 outcome (within-fit null) | inconclusive | fails (3.02) | fails (7.15) |
| averaged ratios spectral / basin / class | 1.81 / 1.49 / 1.47 | 4.39 / 1.96 / 2.63 | 4.09 / 1.33 / 1.66 |
| averaged outcome | inconclusive | fails | fails |
| halves vs refit, spectral / basin | 0.85 / 0.98 | 1.07 / 0.78 | 0.93 / 0.53 |
| fit noise ÷ within-fit null, spectral | 1.9–2.2 | 2.6–3.2 | **7.6–7.9** |
| refit distances | one tight cluster | broad, lumpy | **distinct clusters per seed** |
| pooled W2 after 5-seed averaging | 0.030 → 0.029 | 3.53 → 5.10 | 15.8 → 7.4 |

- All three layers' halves sit within their own refit noise: none shows a
  real split-half difference.
- They differ in how stable the drift fit is. PCA-30's fits agree with each
  other; FA-30's scatter; FA-10's fall into separate solutions per seed. A
  layer whose operators depend on the seed this much cannot support a
  between-dataset comparison, whatever the threshold.
- Only PCA-30 is not a "fails" on either yardstick. This supports Tom's
  design choice of PCA-30 (PROJECT.md rule 1).

## Quarantined (not evidence)

- Observations: input sha256 values (confirmed).
