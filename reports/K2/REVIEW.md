# K2 review

reviewed_commit: 0af20a1 (task/K2)
date: 2026-10-15
verdict: ACCEPT — decision-rule outcome: **survives** (rule 5)

## Checker

`check_report.py tasks/K2.md reports/K2/report.json --prev reports/K1/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/K2/`. The script
follows the task: count matching on `default_rng(20261014)` over sorted
shared types, 2000-gene subset before `normalize_total` and `log1p`, K1's
draw order on `default_rng(20261015)` / `default_rng(20261020)`.

## Recomputation

Planner's own end-to-end script (count matching, latents, tangent bases,
ten drift fits, spectra, classes, all draws and statistics). PCA latents
**identical** to the deliverables; per-seed eigenvalues agree to ≤ 6e-8;
all 29 gated numbers agree to ≤ 1e-6 relative.

| gated | report | recomputed | match |
|---|---|---|---|
| `basin_ratio_refit` | 0.907431 | 0.907431 | yes |
| `basin_ratio_withinfit` | 1.44644 | 1.44644 | yes |
| `basin_sep_a` | 3.52039 | 3.52039 | yes |
| `basin_sep_b` | 3.37286 | 3.37286 | yes |
| `basin_vs_k1` | 1.04087 | 1.04087 | yes |
| `basin_w1_pooled` | 0.0231399 | 0.0231399 | yes |
| `basin_w1_refit` | 0.0255005 | 0.0255005 | yes |
| `class_tv_pooled` | 0.0755023 | 0.0755023 | yes |
| `class_tv_ratio_refit` | 0.713011 | 0.713011 | yes |
| `class_tv_ratio_withinfit` | 1.10437 | 1.10437 | yes |
| `class_tv_refit` | 0.105892 | 0.105892 | yes |
| `class_tv_vs_k1` | 0.651633 | 0.651633 | yes |
| `count_match_max_abs` | 0 | 0 | yes |
| `depth_match_max_rel` | 0.00996867 | 0.00996867 | yes |
| `frac_rank_ge10` | 0.9613 | 0.9613 | yes |
| `n_cell_types_tested` | 4 | 4 | yes |
| `n_cells_a` | 969 | 969 | yes |
| `n_cells_b` | 969 | 969 | yes |
| `outcome` | survives | survives | yes |
| `rank_pca_a` | 21 | 21 | yes |
| `rank_pca_b` | 20 | 20 | yes |
| `spectral_ratio_refit` | 1.20933 | 1.20933 | yes |
| `spectral_ratio_withinfit` | 2.37648 | 2.37648 | yes |
| `spectral_sep_a` | 8.73915 | 8.73915 | yes |
| `spectral_sep_b` | 9.98935 | 9.98935 | yes |
| `spectral_vs_k1` | 0.790688 | 0.790688 | yes |
| `spectral_w2_pooled` | 0.0374546 | 0.0374546 | yes |
| `spectral_w2_refit` | 0.0309714 | 0.0309714 | yes |
| `symmetry_max` | 2.30366e-07 | 2.30366e-07 | yes |

Figures: `fig_spectra.png` (per tested type, restricted spectra v2 vs v3;
refit and pooled W2 inset) and `fig_basin.png` (per-type median
basin_index) are as specified.

## Decision rule outcome

0. `count_match_max_abs` 0; `depth_match_max_rel` 0.010 ≤ 0.05;
   `frac_rank_ge10` 0.961 ≥ 0.80 → valid.
1. `symmetry_max` 2.3e-07. 2. 4 tested types (B cells, ILC, Monocytes,
   T cells; 969 cells per dataset).
3. Separations 3.37–9.99 (≥ 2).
4. Not met (`spectral_ratio_refit` 1.209, `basin_ratio_refit` 0.907).
5. `spectral_ratio_refit` 1.209 ≤ 1.5, `basin_ratio_refit` 0.907 ≤ 1.5,
   `class_tv_ratio_refit` 0.713 ≤ 1.5 → **survives**.

Secondary reading: `spectral_vs_k1` 0.79, `basin_vs_k1` 1.04,
`class_tv_vs_k1` 0.65 — all ≤ 1.25: the v2/v3 difference is no larger than
K1's run-to-run (same-chemistry) difference at this sample size.

## Deviations

None listed. None found.

## Verdict

**ACCEPT.** Outcome: **survives.** Per K2's pre-registration, the planner
brings a T02 proposal to Tom; T02 stays Tom's decision because K1 did not
pass.

## Planner reading (not findings) — what this does and does not show

- **Shows:** at matched cell count and depth, PCA-30 basin operators of the
  same donor's cells differ between v2 and v3 chemistry by no more than two
  drift fits of one dataset differ from each other (spectral 1.21, basin
  0.91, class 0.71), and no more than two v3 runs differ (K1). Spectra
  overlap closely in all four types; eigenvalue scales (0.363 vs 0.326)
  and class counts (holder 545 vs 569) agree; per-type basin medians lie
  near the identity line in the same order.
- **Does not show:** that operators are identical. "Survives" means no
  difference above fit noise was detectable with ~970 cells and 4 cell
  types. Small chemistry effects below that noise are not excluded.
- **Context the result depends on:**
  - the refit yardstick (rule 12), adopted after K0b but before K1/K2 ran;
    under the T01c within-fit yardstick the spectral ratio is 2.38
    (inconclusive, not failing);
  - depth matching by thinning (T01d) — the result is about chemistry at
    equal depth, not raw v2 vs v3 data;
  - PCA-30 as Tom's design choice (K0/K0b inconclusive for it);
  - K1 inconclusive by 0.03, so the two-run baseline is not itself a pass.
- The T01–T01e "shape differs" history came from the retired FA-30 /
  full-spectrum design; it does not contradict this result.

## Quarantined (not evidence)

- None (report `observations` list input sha256 values only; confirmed).
