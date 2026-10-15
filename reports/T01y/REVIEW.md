# T01y review (exploratory)

reviewed_commit: 9366943 (task/T01y)
date: 2026-10-15
verdict: ACCEPT — reading-rule outcome: "shape differs" (exploratory)

**Exploratory. Reported beside T01c; it does not change T01c's outcome
(INCONCLUSIVE, rule 3) and is not a test of the axiom.**

## Checker

`check_report.py tasks/T01y.md reports/T01y/report.json --prev reports/T01c/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/T01y/`. Input ops
sha256 (ab5238dc… / 7b1fb01d…) match the committed T01c files. The script
compares separations with the written-out T01c reference values, as the
task requires.

## Recomputation

Planner's own script, reading `reports/T01c/ops_*.npz` as committed:

| gated | report | recomputed | match |
|---|---|---|---|
| `scale_v2` | 30.2874 | 30.2874 | yes |
| `scale_v3` | 4.80332 | 4.80332 | yes |
| `scale_ratio` | 6.3055 | 6.3055 | yes |
| `spectral_w2_pooled_s` | 11.778 | 11.778 | yes |
| `spectral_w2_null_s` | 0.75438 | 0.75438 | yes |
| `spectral_ratio_s` | 15.613 | 15.613 | yes |
| `basin_w1_pooled_s` | 1.0558 | 1.0558 | yes |
| `basin_w1_null_s` | 0.19483 | 0.19483 | yes |
| `basin_ratio_s` | 5.4188 | 5.4188 | yes |
| `class_tv_ratio_s` | 6.508188 | 6.508188 | yes |
| `max_sep_change` | 3.61e-06 | 3.61e-06 | yes |

Figure `fig_spectra_scaled.png`: per-type scaled eigenvalue densities v2 vs
v3 with null inset; as specified.

## Reading-rule outcome

1. `class_tv_ratio_s` equals T01c's; `max_sep_change` 3.6e-06 ≤ 1e-4 →
   valid.
2. Not met (ratios > 1.5).
3. `spectral_ratio_s` 15.6 > 3.0 (and `basin_ratio_s` 5.4 > 3.0) →
   **"shape differs"**: one scale per dataset does not explain the v2/v3
   difference.

## Verdict

**ACCEPT.**

## Planner reading (not findings)

- Dividing out the 6.3× scale cuts `spectral_ratio` from 53.7 to 15.6, but
  `basin_ratio` rises from 2.3 to 5.4: on a common scale, v2's large
  within-type spread in basin_index shrinks, so the null shrinks more than
  the between-chemistry distance. The remaining
  difference is several times the within-chemistry null.
- In the figure, v2's scaled spectra have a long negative tail (−57…0) with
  the bulk near 0; v3's are concentrated near 0 with a sparse tail to
  −115. Leading eigenvalues (scaled) cluster near +0.03…+0.06 for every v2
  type but range −0.10…+0.65 across v3 types. v2's wells look uniformly
  shallow-topped; v3's differ by type.
- One untested explanation: sequencing depth. v2 cells have far fewer UMIs
  than v3; a noisier, sparser v2 latent could yield a flatter, more uniform
  learned energy near cells and steeper walls. Depth is technical but is not
  a rotation, so if it drives this, the axiom as written still fails for
  this pair, while a depth-matched comparison would tell whether the wells
  agree once depth is equal.

## Quarantined (not evidence)

- Observations: input ops sha256 values (confirmed).
