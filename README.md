# tahoe-coord

Planner/executor protocol for Tahoe (operator-level single-cell dataset alignment).

- Planner: Claude, using `skills/tahoe-planner/`. Writes tasks, reviews reports.
- Executor: role "Manus" — Manus AI, a fresh Claude Code session, or a human — bound by `coord/MANUS.md`.

## Setup

Planner side (this repo):
```
git init && git add . && git commit -m "protocol: initial"
# push to origin
cp -r skills/tahoe-planner ~/.claude/skills/tahoe-planner   # or symlink
pip install jsonschema
```

Executor side:
```
git clone <origin> ~/Documents/tahoe-manus
cd ~/Documents/tahoe-manus
cp coord/MANUS.md CLAUDE.md        # only if executor is Claude Code
python3.11 -m venv .venv
```

## Loop

1. `/tahoe-planner new T01` → commits `tasks/T01.md`, prints the executor paste.
2. Paste into a **fresh** executor session. Executor pushes `task/T01`, replies with the DONE line.
3. `/tahoe-planner review T01` → checker, recompute, verdict, FINDINGS, merge, push.
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
