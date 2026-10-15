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

---

# T01e review — revision 1

reviewed_commit: cd0223e (task/T01e)
date: 2026-10-15
verdict: ACCEPT (stopped run) — decision-rule outcome: INCONCLUSIVE (rule 0: FA not admissible)

## Checker

`check_report.py tasks/T01e.md reports/T01e/report.json --prev reports/T01d/report.json`
→ PASS (0 fail, 0 warn), status `stopped`, `stop_rule` quotes rule 0.
Branch changes only `reports/T01e/`. Raw 10k file sha256 ebc5dedc… matches.

## Recomputation

Planner's own end-to-end script (download checksum, read, QC, CellTypist on
all 11,769 cells, B-sub draw, depth thinning, gene filter, sorted HVGs, FA,
tangent rank). QC flags and cell-type labels equal the staged file
cell-for-cell.

**Spec gap (planner's):** T01e says thinning goes "cells in file order"
and also "thin A's cells first (types in sorted order), then B-sub's". The
executor read this as: types in sorted order, cells in file order within
each type. Under that reading, the planner's matched files, gene list and
FA latents (≤ 9.5e-7) are identical to the executor's, and:

| gated | report | recomputed | match |
|---|---|---|---|
| `n_cells_b_raw` | 11769 | 11769 | yes |
| `n_cells_b_qc` | 11043 | 11043 | yes |
| `n_cells_a` | 1113 | 1113 | yes |
| `n_cells_bsub` | 1113 | 1113 | yes |
| `count_match_max_abs` | 0 | 0 | yes |
| `depth_a_median` | 6607 | 6607 | yes |
| `depth_b_median` | 6605 | 6605 | yes |
| `depth_match_max_rel` | 0.00158 | 0.00158 | yes |
| `n_genes_detected_both` | 14475 | 14475 | yes |
| `n_hvg_zero_a` | 0 | 0 | yes |
| `n_hvg_zero_b` | 0 | 0 | yes |
| `rank_fa_a` | 8 | 8 | yes |
| `rank_fa_b` | 9 | 9 | yes |
| `frac_fa_at_ceiling` | 0.0 | 0.0 | yes |

Under the other reading (all cells of a dataset in file order) the planner
got `rank_fa_a` 9 and `rank_fa_b` 10 (A still < 10): **rule 0 fires under
both readings**, so the outcome does not depend on the gap.

## Decision rule outcome

0. `count_match_max_abs` 0, `depth_match_max_rel` 0.0016 ≤ 0.05, no
   zero-count HVGs → valid. FA median ranks 8 and 9 < 10 → **not
   admissible → INCONCLUSIVE (stop)**. Rules 1–6 not reached.

## Verdict

**ACCEPT** the stopped run. Outcome: **inconclusive** (no admissible
positional layer for the same-chemistry pair). T02 stays blocked.

## Planner reading (not findings)

- The same 1113 v3 cells reached FA rank 14 in T00c (genes shared with v2)
  and 8 here (genes shared with 10k v3, after depth matching). The
  positional layer's local dimensionality moves by ~6 with the gene list
  and depth, enough to decide by itself whether the kill test can run.
- The test design cannot now separate "chemistry reshapes wells" from
  "operators do not reproduce across runs": the control that would decide
  it did not get past the positional-layer gate.
- Recommendation: pause the ladder and redesign before further patches.

## Quarantined (not evidence)

- "The pre-registered stop rule was triggered after FA tangent-rank
  evaluation." — restates the outcome.
