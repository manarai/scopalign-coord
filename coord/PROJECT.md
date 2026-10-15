# PROJECT.md — scOpAlign (operator-level dataset alignment)

Owner: Tom Terooatea (planner / reviewer)
Executor role: "Manus" (see MANUS.md) — a role name, not a vendor.

## What this project is

scOpAlign aligns single-cell datasets on their **basin operators** rather
than on cell positions. The basin operator is scJDO2's symmetric latent
Hessian: for each cell, `J_c = −β · Hess E_θ(z_c)`, where `E_θ` is a scalar
energy fit by exact score matching on the dataset's own latent coordinates.
Its eigenvalues are real and basis-free; its leading eigenvalue defines the
cell's `basin_index`.

Axiom under test: *batch effects move cells, they do not reshape the local
wells they sit in.* If that holds, eigenvalue spectra and `basin_index`
agree across batches before any alignment, and alignment is a single
orthogonal change of basis that never touches eigenvalues.

If the axiom fails on a pure technical replicate, the project stops (T01).

## What scJDO2 actually provides (verified from source, Oct 14 2026)

- `scjdo2.energy_drift.EnergyDriftAdapter` — energy-parameterised score,
  Hyvärinen exact score matching by default (`loss="esm"`), no pseudotime.
- `.jacobian_autodiff(z)` — per-cell `(d, d)` symmetric Hessian in latent.
- `scjdo2.opspace.classify_basin_latent(rule="lead_eig")` — writes
  `op_basin_index = −lead_re` and `op_behaviour_basin ∈ {holder, fork, drift}`.
- Gene-space Jacobians (Phase C) exist but are rank-limited to `d` and are
  **not used** by scOpAlign.
- Known limitation recorded in scJDO2 SUMMARY §5.7: PCA latents saturate
  the tangent rank at `d`. This is why the positional layer here is not PCA.

## Fixed design decisions

1. **Positional layer is a linear-decoder model.** Default: Factor Analysis
   (`sklearn.decomposition.FactorAnalysis`, `d = 30`, fit **separately**
   per dataset). LDVAE (scvi-tools, linear decoder) is the authorised
   upgrade if T00b shows FA still saturates tangent rank. PCA is kept only
   as the T00b comparator. scVI and any nonlinear decoder are excluded:
   the alignment map must carry back to cell positions exactly, which
   requires a linear decoder.
2. `d` is fixed across datasets within a task.
3. **Operators** are scJDO2 latent symmetric Hessians, autodiff, from an
   energy drift fit with ESM at fixed seed and epochs.
4. **Alignment** is an orthogonal map `A` with `J_v3 ≈ A J_v2 Aᵀ` on
   matched neighbourhoods (Procrustes on eigenvector frames, or direct
   minimisation of the Frobenius mismatch; T02 pre-registers which).
   Eigenvalues are untouched by construction. Locally varying `A`
   (orthogonal per neighbourhood, graph-smoothness penalised) only under
   T06.
5. **Correction** lands in the latent embedding: `z ↦ A z`. Operators are
   never corrected; they are the evidence used to find `A`. Whatever
   positional difference remains after `A` is the candidate biology.
6. Gene-space lift of the correction (through the FA loadings) is out of
   scope until after T03.
7. **Primary outcome** at every stage is agreement of eigenvalue spectra
   and `basin_index` **before** alignment. Eigenvector-frame agreement
   after alignment is secondary.
8. No RNA velocity, no velocyto, no scVelo anywhere.
9. Any dataset or parameter not named in a task file is out of scope.

## Known risks (pinned so task writers see them)

- **Identifiability.** Count-level batch distortion is nonlinear and bends
  the learned energy. Operators alone cannot separate "chemistry bent the
  well" from "biology changed the well." T01 bounds the chemistry effect
  on a pair with no biology; T04 is the only task with both a known batch
  and a known biological difference, so it runs before T03.
- **Correspondence.** Procrustes needs matched pairs. Cell-type labels
  are the weak link; spectral-similarity matching via optimal transport
  is the reserve.
- **Tangent rank.** A linear decoder may still saturate rank at `d`.
  T00b measures this before any operator is compared.
- **Local freedom.** Locally varying `A` can bend anything into anything.
  T06 must report residual removed beyond the global `A`; a large value is
  a warning, not a success.

## Task ladder

| ID   | Purpose | Status |
|------|---------|--------|
| T00  | Environment, scJDO2/scJDO installs, data staging | ACCEPT (2026-10-14) |
| T00b | Positional layer: PCA vs FA vs LDVAE tangent-rank on PBMC pair; pick FA unless it saturates | ACCEPT (2026-10-14) — FA by rule; held, see reports/T00b/REVIEW.md |
| T01  | Kill test: eigenvalue spectra + basin_index agree across v2/v3 chemistry before alignment? | ACCEPT (2026-10-14) — outcome INCONCLUSIVE (rule 3: operators do not separate cell types); T02 stays blocked |
| T00c | Positional layer redone on genes detected in both chemistries; adds a rank floor | ACCEPT (2026-10-14) — FA (ranks 11/14) |
| T01c | Kill test rerun on the T00c layer (T01 design unchanged) | ACCEPT (2026-10-14) — outcome INCONCLUSIVE (rule 3: basin_sep_v2 1.22); T02 stays blocked |
| T02  | Global orthogonal `A`; eigenvector-frame agreement after alignment; pre-register estimator | blocked by T01c |
| T04  | Simulator: known nonlinear batch + known biological difference; does `A` remove one and keep the other? | blocked by T02 |
| T03  | Two-lab dataset: carry `A` to `z`; residual = biology check | blocked by T04 |
| T05  | Benchmark vs Harmony / scVI on scib metrics + preservation of a known biological difference | blocked by T03 |
| T06  | Locally varying orthogonal `A` (only if T02 residual is uneven) | conditional |

## Folders (all under ~/Documents)

```
~/Documents/scopalign-coord/   planner side — this repo
~/Documents/scopalign-manus/   executor side — clone + .venv + data/ + reports/
~/Documents/scJDO2/            scJDO2 source, editable install, read-only to executor
~/Documents/scJDO/             scJDO v1 repo, read-only; never checked out
~/Documents/scopalign-manus/vendor/scJDO-v1-r25/
                               git worktree of scJDO at tag scjdo-v56-submission-20260922-r25
                               (commit cbfd591), editable install — scJDO2's own install route
                               (scripts/install_v1_compare.sh), so the live scJDO HEAD never moves
```

Every analysis runs inside `~/Documents/scopalign-manus/`. Planner never
edits that folder. Executor never edits `coord/`, `tasks/`, `templates/`,
`checker/`, `skills/`.

## Known-failure allowlist

Pinned by test ID. Checker fails if the set changes in either direction.

```
KNOWN_FAILURES = [
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.1]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.3]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_saddle_localisation[0.5]",
    "tests/test_phase_b_pitchfork.py::test_pitchfork_half_split_reproducibility",
]
```

(Set from T00's report on ACCEPT; thereafter edited only via an accepted task.)

## Environment

- Executor runs on Tom's Mac in the project `.venv`. Data never leaves the
  machine.
- Every report commits `requirements.lock`; checker flags diffs against the
  previous accepted task.
- Git is local only: the executor folder is a clone of this repo and pushes
  task branches back here. No GitHub remote.
- Backup on every ACCEPT: planner's merge step copies the repo to the
  backup location below. Not optional.

  BACKUP_LOCATION = ~/Library/Mobile Documents/com~apple~CloudDocs/scopalign-backup

  Backup is a single git bundle per ACCEPT (iCloud syncs one file reliably,
  a live `.git` folder it does not):
  ```
  mkdir -p "$HOME/Library/Mobile Documents/com~apple~CloudDocs/scopalign-backup"
  git bundle create "$HOME/Library/Mobile Documents/com~apple~CloudDocs/scopalign-backup/scopalign-coord-$(date +%Y%m%d)-<TASK_ID>.bundle" --all
  ```
  Restore with `git clone <bundle> scopalign-coord`.

## Vocabulary

- **Gated number**: a value named in a task's pre-registration. Only gated
  numbers can enter FINDINGS.md.
- **Deviation**: anything the executor did that the task file did not say.
  Must be listed; unlisted → REJECT.
- **Observation**: anything noticed but not asked for. Quarantined; never
  enters FINDINGS until a task tests it.
- **ACCEPT / REVISE / REJECT**: the planner's only three verdicts.
