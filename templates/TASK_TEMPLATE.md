# <TASK_ID> — <short title>

status: pre-registered | running | ACCEPT | REVISE | REJECT
blocked_by: <task ids or none>
seed: <int>
planner_commit: <sha of main when this task was written>

## Purpose

One paragraph. What question does this task answer, and what happens to the
project under each answer.

## Scope

### In scope
- Dataset(s), by exact file path or accession.
- Cell subset rule (or "all cells").
- Gene set rule (e.g. "top 2000 HVGs, seurat_v3 flavor, batch_key=None").
- Number of latent components.
- kNN size for operator estimation.
- Every parameter the executor needs.

### Out of scope
- Everything not listed above. Name the tempting ones explicitly.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `name_1` | exact definition | `python scripts/<x>.py --flag` |

### Decision rule
Written before any data is touched. Example:
- If `spectral_w2 < 0.10` and `stability_agreement > 0.90` → axiom survives.
- Otherwise → axiom fails for this replicate pair; T02–T05 are cancelled.

### What would change my mind
The one result that would make the planner abandon or revise the axiom.

## Procedure

Numbered steps. Each step that produces output names the output file and
the single command that regenerates it.

1. ...
2. ...

## Deliverables

Exact filenames under `reports/<TASK_ID>/`:
- `report.json`
- `requirements.lock`
- `<figure_1>.png`
- ...

## Checker expectations

- `gated` keys: the exact list from the table above.
- `deliverables`: the exact list above.
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
