# T01d review — attempt 1

reviewed_commit: none (nothing committed; `task/T01d` in the executor folder
is still at 4c08220 with `reports/T01d/` untracked)
date: 2026-10-15
verdict: REVISE

## Checker

Run read-only by the planner on the uncommitted files:
`check_report.py tasks/T01d.md reports/T01d/report.json --prev reports/T01y/report.json`
→ FAIL (1): `schema: '<string>' is not of type 'object'` — the single
`deviations` entry is a plain string; the schema requires
`{"what": ..., "why": ...}`.

## Recomputation

Not done: the protocol reviews committed work only. Gated numbers will be
recomputed on the revised, committed run.

## Required changes (binding)

1. In `reports/T01d/report.json`, write the deviation as an object:
   `{"what": "Plasma cells (v3-thin, 8 cells, not a tested type) have a tied
   majority class (fork=4, drift=4); written as an empty class in
   basin_table.csv.", "why": "The task does not say how to show a tie for an
   untested type; tested-type statistics are unaffected."}`
   This deviation is accepted: it cannot move a gated number.
2. Do not rerun anything else and do not change any number or file other
   than `report.json` (and its own checksum entries if the checker needs them).
3. Run the checker; it must PASS.
4. `git add reports/T01d && git commit` on `task/T01d`, then
   `git push origin task/T01d`. Reply with the DONE line.

## Quarantined (not evidence)

- Observation: thinned input sha256 f0385115… (will be checked on review).

---

# T01d review — revision 1

reviewed_commit: 721ad25 (task/T01d)
date: 2026-10-15
verdict: ACCEPT (run) — decision-rule outcome: INCONCLUSIVE (rule 3)

## Checker

`check_report.py tasks/T01d.md reports/T01d/report.json --prev reports/T01y/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/T01d/`. The deviation
is now a `{what, why}` object, as required. Thinned input sha256 f0385115…
matches the reported value.

## Recomputation

Planner's own end-to-end script: thinning, gene filter, HVGs, FA, tangent
ranks, drifts, Jacobians, classification, draws.

- Thinned counts: identical to `data/pbmc_v2_v3/pbmc_v3_thin.h5ad`.
- **Spec gap (planner's):** three genes tie exactly at the 2000th HVG
  normalised variance, so which ones are selected depends on the column
  order passed to `highly_variable_genes`. T00c/T01d say "genes sorted
  alphabetically" without saying whether before or after selection. The
  executor sorted before selection; the planner first did not, and got a
  list differing by 2 genes, with every downstream number ~1% off (outcome
  unchanged: `spectral_ratio` 30.45, `basin_sep_v2` 1.70, still rule 3).
  With sorting before selection, the planner's gene list, both FA latents
  and both Jacobian tensors are **identical** to the deliverables, and:

| gated | report | recomputed | match |
|---|---|---|---|
| `depth_v2_median` | 3366 | 3366 | yes |
| `depth_v3_median` | 7017 | 7017 | yes |
| `depth_v3thin_median` | 3378 | 3378 | yes |
| `depth_match_max_rel` | 0.00726 | 0.00726 | yes |
| `n_genes_detected_both` | 12512 | 12512 | yes |
| `n_hvg_zero_v2` | 0 | 0 | yes |
| `n_hvg_zero_v3` | 0 | 0 | yes |
| `rank_fa_v2` | 12 | 12 | yes |
| `rank_fa_v3` | 15 | 15 | yes |
| `frac_fa_at_ceiling` | 0.0 | 0.0 | yes |
| `n_cell_types_tested` | 4 | 4 | yes |
| `symmetry_v2` | 9.46e-08 | 9.46e-08 | yes |
| `symmetry_v3` | 8.59e-08 | 8.59e-08 | yes |
| `spectral_w2_pooled` | 301.07 | 301.07 | yes |
| `spectral_w2_null` | 9.804 | 9.804 | yes |
| `spectral_ratio` | 30.71 | 30.71 | yes |
| `basin_w1_pooled` | 4.409 | 4.409 | yes |
| `basin_w1_null` | 2.498 | 2.498 | yes |
| `basin_ratio` | 1.765 | 1.765 | yes |
| `class_tv_pooled` | 0.5462 | 0.5462 | yes |
| `class_tv_null` | 0.05678 | 0.05678 | yes |
| `class_tv_ratio` | 9.620 | 9.620 | yes |
| `spectral_sep_v2` | 3.790 | 3.790 | yes |
| `spectral_sep_v3` | 13.06 | 13.06 | yes |
| `basin_sep_v2` | 1.719 | 1.719 | yes |
| `basin_sep_v3` | 2.841 | 2.841 | yes |
| `basin_index_rho` (secondary) | 0.80 | 0.80 | yes |
| `basin_class_agreement` (secondary) | 0.0 | 0.0 | yes |

Figures `fig_spectra.png` and `fig_basin.png`: as specified.

## Decision rule outcome

0. Validity: depth matched within 0.7% (≤ 5%); no zero-count HVGs; FA
   admissible (12 / 15, gap 3, ceiling 0) → valid.
1. Symmetry ≤ 1e-4 → valid run. 2. 4 types ≥ 3.
3. `basin_sep_v2` = 1.72 < 2 → **INCONCLUSIVE.** Rules 4–6 not reached.

## Deviations

Plasma cells tie (untested type) — accepted in attempt 1; cannot move a
gated number.

## Verdict

**ACCEPT** the run. Outcome: **inconclusive**. T02 stays blocked.

## Planner reading (not findings)

- Matching depth halves the scale gap (median |λ| 21.9 vs 6.8, ratio 3.2;
  T01c: 6.3) and lowers `spectral_ratio` from 53.7 to 30.7, but it stays
  ten times the 3.0 line, with spectra that do separate cell types in both
  datasets (3.79, 13.06).
- v2 still has no holders (0 / 992) while v3-thin has 390 / 1113; v2's
  basin_index still barely separates types. So v2's wells stay
  uniformly shallow-topped at equal depth: not a depth artefact.
- Across T01, T01c and T01d the same guard (v2 basin separation < 2) has
  blocked a verdict while every between-chemistry ratio sat far above 3.
- Follow-up recorded: future tasks must say "sort before HVG selection".

## Quarantined (not evidence)

- Observation: thinned input sha256 (confirmed).
