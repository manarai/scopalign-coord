# T00 review — attempt 1

reviewed_commit: 247adc1 (task/T00)
date: 2026-10-14
verdict: REVISE

## Checker

`python checker/check_report.py tasks/T00.md reports/T00/report.json` →
FAIL (15 fail): every deliverable missing, no lockfile, no regeneration
commands. All of these follow from `status: blocked`; none is a separate defect.

## Recomputation

None. `gated` is empty. Confirmed independently: `python3.11` is not on
PATH (`python3` is 3.13.2, from miniconda); Homebrew is present
(`/opt/homebrew/bin/brew`) and offers `python@3.11` 3.11.17 (bottled).

## Decision rule outcome

Not reached. The block was correct: T00 names Python 3.11 and does not
authorise installing it.

## Verdict

**REVISE.** Binding answer to the blocker, treated as part of T00:

1. Install Python 3.11 with `brew install python@3.11`. Create the venv with
   `$(brew --prefix python@3.11)/bin/python3.11 -m venv .venv`. Do not use
   miniconda's Python or any conda environment. Record
   `.venv/bin/python --version` in the report under `deviations` (as an
   authorised install, citing this review).
2. Delete the blocked `reports/T00/report.json` contents and run T00 in full
   from Procedure step 1 on the same branch `task/T00`, in a fresh session.
3. Everything else in `tasks/T00.md` (c765636) is unchanged.

## Quarantined (not evidence)

- Observation: "python3.11 is not present on PATH; python3 reports Python
  3.13.2." Confirmed; addressed by item 1.

## Revision 2 (2026-10-14) — supersedes items 1–3 above

Tom changed the velocity source to scJDO2's drift field before any data
was staged. `tasks/T00.md` has been re-issued as revision 2 under the same
id (see its amendment log). Binding:

1. Run `tasks/T00.md` revision 2 in full, from Procedure step 1, on branch
   `task/T00`, in a fresh session. First bring the amended task onto the
   branch: `git fetch origin && git checkout task/T00 && git merge --no-edit
   origin/main`. Overwrite the blocked `reports/T00/report.json`.
2. Python 3.11 from `brew install python@3.11` remains authorised; record it
   under `deviations` citing this review.
3. Do not download BAMs, the GTF, or install velocyto/scVelo/samtools.

## Revision 3 (2026-10-14) — supersedes Revision 2

Tom rewrote the project around scJDO2's basin operator (956c901, ported
here). `tasks/T00.md` is re-issued as revision 3 under the same id; nothing
was staged in earlier attempts. Binding:

1. On branch `task/T00`, in a fresh session:
   `git fetch origin && git checkout task/T00 && git merge --no-edit origin/main`,
   then run `tasks/T00.md` revision 3 in full from Procedure step 1.
   Overwrite the blocked `reports/T00/report.json`.
2. Python 3.11 from `brew install python@3.11` is authorised; record it
   under `deviations` citing this review.
3. scJDO v1 is installed only through the worktree route in the task
   file. Never check out or tag `~/Documents/scJDO`.

---

# T00 review — attempt 3 (revision 3)

reviewed_commit: 9ca4e14 (task/T00)
date: 2026-10-14
verdict: ACCEPT

## Checker

`python checker/check_report.py tasks/T00.md reports/T00/report.json` →
FAIL (1): `known_failures changed` (4 newly failing). This is the
pre-registered T00 exception (no baseline yet); overridden. No other failure.

## Recomputation

| gated | report | recomputed | how | match |
|---|---|---|---|---|
| `scjdo2_sha` | 949fb9c… | 949fb9c… | `git -C ~/Documents/scJDO2 rev-parse HEAD` (unchanged since install) | yes |
| `scjdo_tag` | scjdo-v56-submission-20260922-r25 | same | `git describe --tags` in the worktree; worktree HEAD cbfd591 = `r25^{commit}` | yes |
| `n_failing_tests` | 4 | 4 | both pytest commands rerun by the planner (scJDO2 4 failed / 38 passed; v1 155 passed) | yes |
| `n_cells_v2_raw` | 996 | 996 | `t00_stage.py --count`, and by independent read of the raw `.h5` | yes |
| `n_cells_v3_raw` | 1222 | 1222 | same | yes |

Independent checks (planner's own code, not the executor's script):
both raw `.h5` sha256 match `data_checksums.txt`; staged `.X` equals the raw
`.h5` counts exactly (int32), same barcodes; `qc_pass` recomputed from the
raw files matches cell-for-cell (v2 992 / v3 1113 pass); `cell_type`
present; `~/Documents/scJDO` HEAD (a0b8ef4) and working tree unchanged.
Staging script follows the spec (URLs, QC thresholds, CellTypist model,
no majority voting, header-only `label_map.csv`).

The four scJDO2 failures are `XPASS(strict)` documented negative controls
in `test_phase_b_pitchfork.py`, not broken assertions. Pinned by ID in
`KNOWN_FAILURES` (PROJECT.md, this merge).

Correction: T00 and PROJECT.md gave the tag's hash as `a47729b`; that is the
annotated tag object. The tagged commit is `cbfd591`. PROJECT.md fixed;
T00's text left as issued.

## Decision rule outcome

Both packages import, `EnergyDriftAdapter` and `classify_basin_latent`
import, `scjdo_tag` = `scjdo-v56-submission-20260922-r25`, both `.h5ad`
have `qc_pass` and `cell_type` → **ACCEPT**.

## Deviations

1. Homebrew Python 3.11.17 — authorised (Revision 3 item 2). Cannot move a gated number.
2. `PYTHONDONTWRITEBYTECODE=1` for pytest — keeps the checkouts unwritten as
   the task requires. Cannot move a gated number. Accepted.

## Verdict

**ACCEPT.** `KNOWN_FAILURES` set from this report.

## Quarantined (not evidence)

- "scJDO2 pytest reported 4 failing test IDs." — restates a gated number.
- "scJDO v1 pytest reported no failing test IDs." — consistent with the planner's rerun.
