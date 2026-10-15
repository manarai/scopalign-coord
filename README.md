# scopalign-coord

Planner/executor protocol for scOpAlign (operator-level single-cell dataset alignment).

- Planner: Claude, using `skills/scopalign-planner/`. Writes tasks, reviews reports.
- Executor: role "Manus" — Manus AI, a fresh Claude Code session, or a human — bound by `coord/MANUS.md`.

## Setup

Planner side:
```
# unpack to ~/Documents/scopalign-coord, then
cd ~/Documents/scopalign-coord
git remote add origin <github url> && git push -u origin main
ln -s ~/Documents/scopalign-coord/skills/scopalign-planner ~/.claude/skills/scopalign-planner
pip install jsonschema
```

scJDO2 and scJDO v1 are expected at `~/Documents/scJDO2` and `~/Documents/scJDO` (the latter on the preprint tag).

Executor side:
```
git clone <origin> ~/Documents/scopalign-manus
cd ~/Documents/scopalign-manus
cp coord/MANUS.md CLAUDE.md        # only if executor is Claude Code
python3.11 -m venv .venv
```

## Loop

1. `/scopalign-planner new T01` → commits `tasks/T01.md`, prints the executor paste.
2. Paste into a **fresh** executor session. Executor pushes `task/T01`, replies with the DONE line.
3. `/scopalign-planner review T01` → checker, recompute, verdict, FINDINGS, merge, push.
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
| push on every ACCEPT | one-laptop risk |
| skill inside the repo | unversioned protocol |
