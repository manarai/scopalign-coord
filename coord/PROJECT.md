# PROJECT.md — scOpAlign (operator-level dataset alignment)

Owner: Tom Terooatea (planner / reviewer)
Executor role: "Manus" (see MANUS.md) — a role name, not a vendor. Any agent or
person bound by MANUS.md and the checker can fill it.

## What this project is

scOpAlign aligns single-cell datasets on their **local dynamics** (per-cell operators
from scJDO2's neighbourhood operator field) rather than on cell positions.

Axiom under test: *batch effects change where cells are, not how they move.*

If the axiom fails on a pure technical replicate, the project stops. That is
the point of task T01.

## Fixed design decisions (do not revisit inside a task)

1. Positional layer is **linear** (PCA or FactorAnalysis). No scVI, no UMAP,
   no nonlinear embedding anywhere upstream of operator estimation.
   Reason: a linear map keeps latent Jacobians a similarity transform of
   gene-space Jacobians, and lets the alignment map be carried back to the
   embedding exactly.
2. Number of latent components is fixed across datasets within a task.
3. Operators are per-cell Jacobians from kNN balls (scJDO2 estimator).
4. Alignment map A is a **change of basis only**: orthogonal Procrustes
   first; locally linear A's with a graph-smoothness penalty only if a task
   explicitly authorises it. Eigenvalues must never be altered by alignment.
5. Spectral / stability-class agreement is checked **before** alignment and
   is the primary outcome. Eigenvector-frame agreement after alignment is
   secondary.
6. Reference libraries: scanpy, scikit-learn, scipy, POT (Python Optimal
   Transport), scJDO2 (local checkout `~/Documents/scJDO2`, installed
   editable, commit sha recorded per task). scJDO v1 is not used.
7. Any dataset or parameter not named in a task file is out of scope.

## Task ladder

| ID  | Purpose | Status |
|-----|---------|--------|
| T01 | Kill test: technical replicate, spectra agree pre-alignment? | pre-registered |
| T02 | Global Procrustes alignment on T01 data; frame agreement post-alignment | blocked by T01 |
| T03 | Back-translation to embedding; residual = biology check on a two-lab dataset | blocked by T02 |
| T04 | Simulator ground truth: inject nonlinear batch distortion, recover operators | blocked by T02 |
| T05 | Benchmark vs Harmony / scVI on scib metrics + preservation of a known biological difference | blocked by T03 |
| T06 | Locally linear A's (only if T02 residual is uneven) | conditional |

## Repository layout

```
scopalign-coord/
  coord/PROJECT.md        this file
  coord/MANUS.md          executor instructions (copy to executor folder as CLAUDE.md if executor is Claude)
  coord/FINDINGS.md       accepted results only; appended by planner on ACCEPT
  tasks/Txx.md            one pre-registered task per file
  templates/              task and report templates
  checker/check_report.py mechanical gate
  reports/Txx/            executor output: report.json, lockfile, figures
  skills/scopalign-planner/   Claude planner skill
```

## Folders (all under ~/Documents)

```
~/Documents/scopalign-coord/   planner side — this repo
~/Documents/scopalign-manus/   executor side — clone of this repo + .venv + data/ + reports/
~/Documents/scJDO2/            scJDO2 source (editable install into the executor .venv)
```

Every analysis runs inside `~/Documents/scopalign-manus/`. Nothing is read
from or written to any other location except the scJDO2 checkout, which
is read-only to the executor. Planner never edits the executor folder.
Executor never edits `coord/` or `tasks/`.

## Known-failure allowlist

Pinned by test ID. Checker fails if the set changes in either direction.

```
KNOWN_FAILURES = []
```

(Populate after the first environment build in T00; thereafter edit only via
an accepted task.)

## Environment

- Executor runs on Tom's Mac, in the project `.venv`. Large `.h5ad` files do
  not leave the machine.
- Every report commits a lockfile (`pip freeze > reports/Txx/requirements.lock`).
  Checker flags diffs against the previous accepted task.
- Local-first: no remote required. The executor folder is a local clone of
  this repo and pushes task branches back here. Add a backup remote later.

## Vocabulary

- **Gated number**: a value named in a task's pre-registration. Only gated
  numbers can enter FINDINGS.md.
- **Deviation**: anything the executor did that the task file did not say.
  Must be listed. An unlisted deviation is grounds for REJECT.
- **Observation**: anything the executor noticed that was not asked for.
  Quarantined. Never enters FINDINGS.md until a new task tests it.
- **ACCEPT / REVISE / REJECT**: the planner's only three verdicts.
