# scopalign-coord

Planner/executor protocol for scOpAlign (operator-level single-cell dataset alignment).

- Planner: Claude, using `skills/scopalign-planner/`. Writes tasks, reviews reports.
- Executor: role "Manus" — Manus AI, a fresh Claude Code session, or a human — bound by `coord/MANUS.md`.

## Setup

Local-first: no GitHub remote is needed. Git is used only on this machine
(task branches, lockfile diffs, merge on ACCEPT). A remote can be added later
for backup.

Planner side:
```
cd ~/Documents/scopalign-coord
ln -s ~/Documents/scopalign-coord/skills/scopalign-planner ~/.claude/skills/scopalign-planner
pip install jsonschema
```

scJDO2 is expected at `~/Documents/scJDO2`; scJDO v1 at `~/Documents/scJDO` (installed from a worktree at its preprint tag, see PROJECT.md).

Executor side (a local clone; its `origin` is the planner folder):
```
git clone ~/Documents/scopalign-coord ~/Documents/scopalign-manus
cd ~/Documents/scopalign-manus
cp coord/MANUS.md CLAUDE.md        # only if executor is Claude Code
python3.11 -m venv .venv
```

## Loop

1. `/scopalign-planner new T00` → commits the task, prints the executor paste. Order: T00 → T00b → T01 → T02 → T04 → T03 → T05 (→ T06 if needed).
2. Paste into a **fresh** executor session. Executor pushes `task/T01` back to the planner folder, replies with the DONE line.
3. `/scopalign-planner review T01` → checker, recompute, verdict, FINDINGS, merge, backup copy.
4. Repeat.

## Why the pieces exist

| piece | protects against |
|---|---|
| pre-registration in task file | moving the goalposts after seeing data |
| checker | scope creep, missing deliverables, silent env drift |
| planner recomputes gated numbers | trusting the executor's typing |
| observations quarantined | executor narrative leaking into findings |
| fresh session per task | executor inheriting planner reasoning |
| lockfile per report | "it changed but I don't know why" |
| known-failure set pinned both directions | allowlist rot |
| backup copy on every ACCEPT | one-laptop risk |
| skill inside the repo | unversioned protocol |
