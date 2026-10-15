# K1 review — attempt 1

reviewed_commit: 5cbbd67 (task/K1)
date: 2026-10-15
verdict: REVISE — normalisation order (planner's spec ambiguity); rerun required

## Checker

`check_report.py tasks/K1.md reports/K1/report.json --prev reports/K0d/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/K1/`.

## Recomputation

The planner's independent recompute does **not** reproduce the report: 19
of 24 gated numbers differ (e.g. `spectral_ratio_refit` 1.725 reported vs
1.529 recomputed; `spectral_sep_a` 12.4 vs 22.0), and the PCA latents differ
by up to 58.

Cause, confirmed exactly: the script (`norm_pca`, line 33) runs
`normalize_total(1e4)` and `log1p` on **all ~33,500 genes** and then subsets
to the 2000 HVGs. Doing that reproduces the executor's latent bit for bit
(max |Δ| = 0.0). The planner normalised the 2000-gene subset, the
convention used and verified in T00c, T01c, T01d, T01e and K0.

K1's wording — "`normalize_total(target_sum=1e4)` then `log1p` on the 2000
genes" — can be read either way. **This is the planner's spec ambiguity,
not an executor error.** It changes every downstream number, so the run
cannot stand.

No numbers from attempt 1 enter FINDINGS.

## Required changes (binding)

1. In `k1_run.py`, the positional layer is computed as: **subset each
   matched file to the 2000 genes of `reports/T01e/genes.txt`, in that
   order, first; then** `sc.pp.normalize_total(target_sum=1e4)` **on that
   2000-gene matrix; then** `sc.pp.log1p`; then PCA as specified. Nothing
   else in the task changes.
2. Rerun all stages on the same branch `task/K1` (in a fresh session),
   overwriting `reports/K1/`; record this change under `deviations` as a
   `{"what", "why"}` object citing this review.
3. Checker must PASS; commit `reports/K1`; `git push origin task/K1`;
   reply with the DONE line.

## Quarantined (not evidence)

- Observations: input sha256 values (A_matched 3106ea2e…, B_matched
  29b81d31…, genes 4a30c174…) — will be checked on the revised run.

---

# K1 review — revision 1

reviewed_commit: af89531 (task/K1)
date: 2026-10-15
verdict: ACCEPT — decision-rule outcome: **inconclusive** (rule 6)

## Checker

PASS (0 fail, 0 warn). Branch changes only `reports/K1/` (plus this review
via the merge of main). The script now subsets to the 2000 genes before
`normalize_total` and `log1p`; the change is declared as a deviation citing
this review. Latents and operators changed, so all stages were rerun.

## Recomputation

Planner's independent script (unchanged from attempt 1): PCA latents
**identical** to the deliverables; per-seed eigenvalues agree to ≤ 1.9e-7;
all 24 gated numbers agree to ≤ 1e-6 relative.

| gated | report | recomputed | match |
|---|---|---|---|
| `basin_ratio_refit` | 0.871801 | 0.871801 | yes |
| `basin_ratio_withinfit` | 1.33373 | 1.33373 | yes |
| `basin_sep_a` | 4.30385 | 4.30385 | yes |
| `basin_sep_b` | 4.29259 | 4.29259 | yes |
| `basin_w1_pooled` | 0.0194761 | 0.0194761 | yes |
| `basin_w1_refit` | 0.0223401 | 0.0223401 | yes |
| `class_tv_pooled` | 0.10054 | 0.10054 | yes |
| `class_tv_ratio_refit` | 1.09419 | 1.09419 | yes |
| `class_tv_ratio_withinfit` | 1.54848 | 1.54848 | yes |
| `class_tv_refit` | 0.0918849 | 0.0918849 | yes |
| `frac_rank_ge10` | 0.921833 | 0.921833 | yes |
| `n_cell_types_tested` | 4 | 4 | yes |
| `n_cells_a` | 1113 | 1113 | yes |
| `n_cells_b` | 1113 | 1113 | yes |
| `outcome` | inconclusive | inconclusive | yes |
| `rank_pca_a` | 19 | 19 | yes |
| `rank_pca_b` | 19 | 19 | yes |
| `spectral_ratio_refit` | 1.52947 | 1.52947 | yes |
| `spectral_ratio_withinfit` | 2.48885 | 2.48885 | yes |
| `spectral_sep_a` | 22.0494 | 22.0494 | yes |
| `spectral_sep_b` | 25.5166 | 25.5166 | yes |
| `spectral_w2_pooled` | 0.0465244 | 0.0465244 | yes |
| `spectral_w2_refit` | 0.0304187 | 0.0304187 | yes |
| `symmetry_max` | 2.37802e-07 | 2.37802e-07 | yes |

Figures: `fig_spectra.png` (per tested type, A vs B; refit and pooled W2
inset) and `fig_basin.png` (per-type median basin_index A vs B) are as
specified.

## Decision rule outcome

0. 1113 = 1113; `frac_rank_ge10` 0.922 ≥ 0.80 → valid.
1. `symmetry_max` 2.4e-07. 2. 4 tested types (B, ILC, Monocytes, T cells).
3. Separations 4.29–25.5 (≥ 2).
4. `spectral_ratio_refit` 1.529, `basin_ratio_refit` 0.872 (≤ 3).
5. Not met: `spectral_ratio_refit` 1.529 > 1.5 (basin 0.872, class TV
   1.094 ≤ 1.5).
6. → **inconclusive**; planner brings the numbers to Tom.

## Deviations

1. Normalisation order per REVIEW item 1 — authorised; the specified
   behaviour.

## Verdict

**ACCEPT.** Outcome inconclusive; K2 not written until Tom decides.

## Planner reading (not findings)

- K1 misses "survives" by 0.03 on one of three ratios; the other two pass
  with room (basin 0.87, class TV 1.09). Within-fit ratios (secondary) are
  2.49 / 1.33 / 1.55.
- By eye the two runs' restricted spectra overlap closely in all four types
  and per-type basin medians lie on the identity line in the same order;
  eigenvalue scale agrees (0.287 vs 0.267); classes are similar (holder 606
  vs 680).
- The comparison is at the power limit of the 1k run: 1113 cells, 4 types.
  No larger same-donor v3 run exists to repeat it with more cells.

## Quarantined (not evidence)

- Observations: input sha256 values (confirmed).
