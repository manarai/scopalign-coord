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
