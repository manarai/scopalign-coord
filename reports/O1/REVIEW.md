# O1 review — pre-finalisation note

date: 2026-10-15

The O1 run completed, but the checker rejected the 12 `cases/case_*.npz`
entries listed under `deliverables` ("deliverable not in task"). This is a
planner spec conflict: O1 asks for the case arrays to be recorded in
`report.json` with sha256, but they are not deliverables.

Binding clarification (no change to any computation or gated number):
record the 12 case arrays under the new optional `artifacts` list
(`[{"path": "cases/case_<n>_<r>.npz", "sha256": ...}]`), added to the
schema and checker on main; keep `deliverables` exactly as the task lists.
Merge `origin/main` first to get the updated checker.
