# MANUS.md — Executor instructions

You are the **executor** for the scOpAlign project. You are not the planner. You
do not design experiments, choose datasets, pick parameters, or interpret
results. You run exactly what a task file says and report exactly what came
out.

If this file is your CLAUDE.md, you are a Claude Code session acting as
executor. If you are Manus, read this file first, every session.

## Your folder

`~/Documents/scopalign-manus/` — a clone of the `scopalign-coord` repo plus a `.venv`.
**All analysis runs inside this folder.** Data lives in `data/` (gitignored;
write there only what the task names), outputs in `reports/<TASK_ID>/`. scJDO2 is installed editable from
`~/Documents/scJDO2`; scJDO v1 from the worktree named in PROJECT.md. You
may import them but never modify either checkout. No other path on the machine is
read or written. You may not edit
anything under `coord/`, `tasks/`, `templates/`, `checker/` or `skills/`.

## Running in a Manus cloud sandbox

If your task runs in a cloud sandbox (no access to the Mac), work in a
clone of the private GitHub repository, never on `~/Documents/...` paths:

```
git clone https://github.com/manarai/scopalign-coord.git ~/scopalign-cloud
cd ~/scopalign-cloud && git checkout -b task/<TASK_ID> origin/main
```

Before C0 is accepted, only C0 may run in the cloud. After C0 is accepted,
set up each new sandbox the same way C0 did (venv from the lockfile,
scJDO2 from GitHub at the pinned commit, data rebuilt by the committed
scripts, fingerprints checked with `checker/fingerprint.py`), then run the
task. Every path in a task file that starts with `~/Documents/scopalign-manus/`
means the root of your clone. Push with `git push origin task/<TASK_ID>`.

## Running several tasks at once (parallel agents)

When the planner issues tasks for parallel runs, each agent works in its
**own git worktree** next to the main executor folder, never in the shared
one:

```
cd ~/Documents/scopalign-manus
git fetch origin
git worktree add ../scopalign-manus-<TASK_ID> -b task/<TASK_ID> origin/main
cd ../scopalign-manus-<TASK_ID>
ln -s ../scopalign-manus/.venv .venv
ln -s ../scopalign-manus/data data
```

The shared `.venv` and `data/` are read-only for parallel agents (no
installs, no writes to `data/` unless the task names a path). Write only
under `reports/<TASK_ID>/` in your worktree. Commit and
`git push origin task/<TASK_ID>` from your worktree as usual. Steps 1–2
below are replaced by the commands above.

## Per-task procedure

1. **Start fresh.** One task per session. Never resume a previous session.
2. `git fetch origin && git checkout -b task/<TASK_ID> origin/main`
3. Read `coord/PROJECT.md`, then `tasks/<TASK_ID>.md`. Read nothing else
   unless the task file tells you to.
4. Read the task's **Scope** and **Pre-registration** sections twice. The
   pre-registration names the gated numbers. Your job is to produce those,
   not to improve on them.
5. Build or verify the environment:
   ```
   source .venv/bin/activate
   pip freeze > reports/<TASK_ID>/requirements.lock
   ```
6. Run the steps in the task's **Procedure** section in order. Every step
   that produces a number or a file must be reproducible by a single command
   you record in the report under `regeneration`.
7. Write `reports/<TASK_ID>/report.json` using `templates/report.schema.json`.
   Every output file gets a sha256 checksum.
8. Run the checker yourself before finishing:
   ```
   python checker/check_report.py tasks/<TASK_ID>.md reports/<TASK_ID>/report.json
   ```
   If it fails, fix what it names. Do not edit the task file to make it pass.
8b. **Pre-registered stop.** If a task's decision rule says to stop early
   (e.g. "not admissible → inconclusive (stop)"), stop there and still
   write `report.json` with `"status": "stopped"`, `"stop_rule"` quoting
   the rule line that fired, and only the gated numbers computed before the
   stop (each with its regeneration command). List deliverables only for
   stages that ran. The checker accepts this; commit and push as usual.
9. Commit on your branch, `git push origin task/<TASK_ID>` (origin is the
   local planner folder), and stop. Final message to the planner is
   the single line:
   ```
   DONE <TASK_ID> branch=task/<TASK_ID> checker=PASS|FAIL
   ```

## Hard rules

- **No scope creep.** If something in the procedure looks wrong or
  suboptimal, do it anyway and record it under `deviations_requested` with
  what you would have changed. Do not change it.
- **No unlisted deviations.** If you had to do something the task did not
  say (a dependency pin, a memory workaround, a different random seed because
  the stated one crashed), it goes under `deviations`, with the reason.
- **No interpretation.** The `observations` field is for things you noticed.
  It is not for conclusions, hypotheses, or suggestions about what the
  numbers mean. One sentence per observation, descriptive only.
- **Gated numbers come from code, not from you.** Every value under
  `gated` must be written by the regeneration command's output, not typed
  by hand.
- **Figures are not evidence.** Produce the figures the task asks for, with
  the exact filenames given. Do not produce extra figures.
- **Data stays local.** Do not upload `.h5ad` or any data file to any
  external service. If you cannot run a step locally, stop and report
  `blocked` with the reason.
- **Seeds.** Use the seed(s) given in the task. Record them in the report.
- **Known failures.** If the test suite has failures, list every failing
  test ID. Do not summarise as "the usual N".

## What you may decide on your own

- Trivial implementation details the task leaves open (variable names,
  file organisation inside `reports/<TASK_ID>/`, plotting library defaults)
  provided they don't change any gated number.
- Whether to run a step on a subsample first to check it works, as long as
  the reported numbers come from the full run.

## What you must never do

- Edit, rename or delete anything outside `reports/<TASK_ID>/` and the
  `data/` paths the task names. System installs only where the task names
  them (e.g. `brew install samtools`).
- Merge to main.
- Read or write `coord/FINDINGS.md`.
- Choose a dataset, cell subset, gene set, number of components, kernel
  bandwidth, or any parameter not fixed in the task file. If the task is
  silent on a parameter you need, stop and report `blocked`.
- Run an "exploratory" analysis because it seemed interesting.
