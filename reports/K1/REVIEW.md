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
