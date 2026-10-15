# T01e review — attempt 1

reviewed_commit: none (nothing committed; `task/T01e` in the executor folder
is at 0feb6c1 with `reports/T01e/` untracked)
date: 2026-10-15
verdict: REVISE

## What happened

The run completed the stage, match and latent stages and then stopped. The
planner, computing T00c's tangent rank read-only on the untracked FA
latents, gets median ranks 8 (A) and 9 (B): below the pre-registered floor
of 10, so decision rule 0 ("FA not admissible → inconclusive (stop)")
applies. Stopping was correct. No `report.json` was written, so the checker
could not run.

That is partly the planner's design gap: the checker required all gated
numbers and had no way to report a pre-registered early stop. Fixed on main
(this commit): report status `"stopped"` with a `stop_rule`, carrying only
the gated numbers reached; `coord/MANUS.md` item 8b says how.

Nothing here is evidence yet: the numbers are recomputed on the committed
report.

## Required changes (binding)

1. `git fetch origin && git merge --no-edit origin/main` on `task/T01e` (to
   get the new checker, schema and MANUS.md item 8b). Do not rerun any
   analysis stage; the existing files under `reports/T01e/` and
   `data/pbmc_10k_v3/` stand.
2. Write `reports/T01e/report.json` with `"status": "stopped"` and
   `"stop_rule": "0. FA not admissible by T00c's rule (median rank ≥ 10 in both, gap ≤ 5, ceiling ≤ 0.20) → inconclusive (stop)"`,
   and in `gated` exactly the numbers the completed stages produce:
   `n_cells_b_raw`, `n_cells_b_qc`, `n_cells_a`, `n_cells_bsub`,
   `count_match_max_abs`, `depth_a_median`, `depth_b_median`,
   `depth_match_max_rel`, `n_genes_detected_both`, `n_hvg_zero_a`,
   `n_hvg_zero_b`, `rank_fa_a`, `rank_fa_b`, `frac_fa_at_ceiling` — each
   produced by its regeneration command, not typed.
3. List as deliverables only files that exist: `requirements.lock`,
   `data_checksums.txt`, `match_table.csv`, `genes.txt`,
   `latent_fa_a.npy`, `latent_fa_b.npy`, `scripts/t01e_run.py`.
4. Deviations as `{"what", "why"}` objects (if none, `[]`).
5. Checker must PASS; commit `reports/T01e`; `git push origin task/T01e`.
   Reply with the DONE line.
