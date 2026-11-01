# K3 — Review

Reviewer: planner (Claude), 2026-11-01. Branch `task/K3` @ 59ea202 (report data from e0b1876), Manus cloud sandbox (T04/T04b clone).

## Checker
`python checker/check_report.py tasks/K3.md reports/K3/report.json --prev reports/T04b/report.json` → **PASS** (0 fail, 0 warn).

## Recomputation
Light recompute (working mode):
- Every ratio in `ratios.csv` (20 rows) equals observed ÷ null. Each gated ratio equals the median over types (tech, genebatch) or the single row (type, bio). All 8 ratios match the report exactly.
- `latent_check_rel` 4.65e-5 = 0.00146 / 31.42 (`latent_check.json`); component |r| min 0.99999998. `symmetry_max` 2.0e-7 matches `fit_symmetry.csv`. `n_types` 4 / 4 from `type_counts.csv` (B cells, ILC, Monocytes, T cells; ≥ 100 cells).
- End-to-end spot-check on the Mac from the raw data: correlation-network distance, monocytes H1 vs H2 (np.corrcoef on log-normalised 2000 genes, zero-variance genes set to 0) = **0.680504**, report 0.680504.
- Not recomputed: the operator networks (fits and per-cell operators live in the cloud cache) and the nulls.

| gated | value |
|---|---|
| n_types_P1 / P2 | 4 / 4 |
| latent_check_rel | 4.65e-5 |
| symmetry_max | 2.02e-7 |
| net_ratio_tech | 1.469 |
| net_ratio_genebatch | **2.212** |
| net_ratio_type | 9.730 |
| net_ratio_bio | 0.455 |
| corr_ratio_tech | 0.829 |
| corr_ratio_genebatch | 0.889 |
| corr_ratio_type | 1.368 |
| corr_ratio_bio | 0.863 |

Figure: `fig_network_ratios.png` shows the four ratios for both network kinds with lines at 1.5 and 3, as asked (the legend overlaps the type bar; cosmetic).
Deviations: none.

## Decision rule outcome
Validity: 4 types ≥ 3, latent 4.65e-5 ≤ 1e-4, symmetry 2.0e-7 ≤ 1e-4 → valid.
Rule 1: `net_ratio_genebatch` 2.21 > 2.0 → **batch_sensitive**. Rules 2–4 not reached. For reference only, rule 2 would also fail: `net_ratio_bio` 0.46 < 1.5.

## Verdict
**ACCEPT** — outcome **batch_sensitive**. The operator-derived gene network separates cell types strongly (9.7× noise) but changes under the gene-level batch (2.2×) and does not see the interferon response (0.46×, i.e. below noise). Network space as built here is not a usable alignment space. The plain correlation network is batch-robust (0.83, 0.89) but barely separates types (1.37) and is also blind to the response (0.86).

## Quarantined (not evidence)
- Executor observation: networks were evaluated on the fly, not saved.
- Planner note (seed only): the high per-type ratios come from monocytes (tech 4.06, genebatch 5.87), the largest group, which has the smallest null (0.17 / 0.14). The other types sit at 1.4–2.4. The ratio may partly scale with group size, because larger groups have tighter nulls. Not tested.
