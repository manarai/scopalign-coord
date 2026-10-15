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
