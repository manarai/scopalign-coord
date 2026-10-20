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

If operators do not reproduce across two halves of one run (K0) or two runs
of one chemistry (K1), the method stops before any chemistry claim. If they
reproduce but fail across chemistries (K2), the axiom fails for chemistry.

## What scJDO2 actually provides (verified from source, Oct 14 2026)

- `scjdo2.energy_drift.EnergyDriftAdapter` — energy-parameterised score,
  Hyvärinen exact score matching by default (`loss="esm"`), no pseudotime.
- `.jacobian_autodiff(z)` — per-cell `(d, d)` symmetric Hessian in latent.
- `scjdo2.opspace.classify_basin_latent(rule="lead_eig")` — writes
  `op_basin_index = −lead_re` and `op_behaviour_basin ∈ {holder, fork, drift}`.
- Gene-space Jacobians (Phase C) exist but are rank-limited to `d` and are
  **not used** by scOpAlign.
- scJDO2 SUMMARY §5.7 records that PCA latents saturate the tangent rank
  at `d`. **Not observed on this project's data:** PCA's median tangent rank
  was 21 / 19 in both T00b and T00c with no cell at the ceiling, while FA's
  moved between 6 and 15 with the gene list and depth (T00b, T00c, T01d,
  T01e). Stable rank does not make PCA's operators better; no PCA operator
  has been fitted yet. Rule 1 therefore lets K0 choose (2026-10-15).

## Fixed design decisions

1. **Positional layer is chosen by K0** (amended 2026-10-15; was FA by
   T00b/T00c). Three linear candidates, each fit **separately** per dataset
   on `normalize_total(1e4)` + `log1p` of the task's HVGs, scores unscaled:
   - **PCA-30:** `PCA(n_components=30, svd_solver="full", random_state=<seed>)`,
     tangent-restricted spectra (rule 3a);
   - **FA-30:** `FactorAnalysis(30, random_state=<seed>, max_iter=2000)`,
     tangent-restricted spectra (rule 3a);
   - **FA-10:** `FactorAnalysis(10, random_state=<seed>, max_iter=2000)`,
     full 10-d spectra (d is already near FA's observed local rank, so no
     restriction).
   K0 selected none (pca30 inconclusive; fa30, fa10 failed). **Layer in use
   from 2026-10-15: PCA-30, by Tom's design choice** — supported by K0b
   Part A (PCA-30's two halves differ no more than two drift fits of the
   same half: 0.85 spectral, 0.98 basin), **not** a K0/K0b pass. K0c
   applies the same check to FA-30; Tom may revisit the choice on its
   result before K1 is issued. scVI and any nonlinear decoder stay
   excluded: the alignment map must carry back to cell positions exactly,
   which requires a linear map.
2. `d` is fixed across datasets within a task (30 for PCA-30).
3. **Operators** are scJDO2 latent symmetric Hessians, autodiff, from an
   energy drift fit with ESM at fixed seed and epochs.
3a. **Tangent-restricted spectra are the primary operator statistic for
   d = 30 layers** (added 2026-10-15). For each cell: `V` = the eigenvectors of the `r = 10`
   largest eigenvalues of the sample covariance of its `k = 40` nearest
   other cells in latent; `J_r = Vᵀ J V` (10 × 10); `λ_r = eigvalsh(J_r)`.
   `basin_index = −max λ_r`; class `fork` / `holder` / `drift` by the
   `lead_eig` rule on `max λ_r` with `tol = 0.1 × median |λ_r|` over that
   dataset. Full 30-d spectra are reported as secondary only. Reason: with
   local rank ≈ 10–20 of 30, the off-tangent directions carry very steep,
   uninformative curvature that dominated every full-spectrum W2 in
   T01–T01d.
4. **Alignment** is an orthogonal map `A` with `J_v3 ≈ A J_v2 Aᵀ` on
   matched neighbourhoods (Procrustes on eigenvector frames, or direct
   minimisation of the Frobenius mismatch; T02 pre-registers which).
   Eigenvalues are untouched by construction. Locally varying `A`
   (orthogonal per neighbourhood, graph-smoothness penalised) only under
   T06.
5. **Correction** lands in the latent embedding: `z ↦ A z`. Operators are
   never corrected; they are the evidence used to find `A`. Whatever
   positional difference remains after `A` is the candidate biology.
6. Gene-space lift of the correction (through the selected layer's loadings) is out of
   scope until after T03.
7. **Primary outcome** at every stage is agreement of tangent-restricted
   eigenvalue spectra and `basin_index` **before** alignment.
   Eigenvector-frame agreement after alignment is secondary.
8. No RNA velocity, no velocyto, no scVelo anywhere.
9. Any dataset or parameter not named in a task file is out of scope.
10. **Kill test is a three-step ladder** (added 2026-10-15), each step issued
   only if the previous one survives: **K0** two random halves of one run
   (estimator reproducibility); **K1** two runs, same chemistry, same donor,
   count- and depth-matched (run reproducibility); **K2** v2 vs v3,
   count- and depth-matched (chemistry). K0 used T01c's statistics and
   rule for all three candidate layers; K1 and K2 use the layer in rule 1
   and the **refit yardstick** (rule 12).
12. **Refit yardstick for K1 and K2** (added 2026-10-15, before either is
   run). K0b showed that two drift fits of the *same* cells, differing only
   in seed, are 1.5–2.3× further apart than T01c's within-fit null, so a
   within-fit null asks for more agreement than the estimator gives itself.
   From K1 on, the decision ratios compare the between-dataset distance with
   a **refit null**: two disjoint cell samples of one dataset, the first read
   under the seed-20261014 drift fit and the second under another seed's fit
   (five seeds per dataset). `ratio_refit = median pooled ÷ median refit
   null`, for spectral W2, basin W1 and class TV. Thresholds unchanged:
   survives if all three ≤ 1.5; fails if spectral or basin > 3.0. The
   within-fit separation guard (≥ 2) stays. Within-fit ratios are reported
   as secondary. **Not retroactive:** K0 and K0b keep their recorded
   outcomes; for reference only, K0b's Part A values (0.85, 0.98) are what
   this yardstick would have shown for PCA-30's halves.
11. **Spec hygiene** (added 2026-10-15) — every task states: genes sorted
   alphabetically **before** `highly_variable_genes`, and the selected list
   kept in that order; for every random draw, the generator, its seed and
   the exact loop order (outer loop, inner loop, which cells consume draws);
   a pre-registered early stop reported as `"status": "stopped"`.

13. **Operator-central alignment** (Tom, 2026-10-15; restates the original
   plan). The rotation `Q` is **estimated from operators and their drift
   field only**, then applied to the latent space (`z_A ↦ Q z_A`). Gene
   loadings and cell positions are never used to estimate `Q`; `Q_load`
   appears only as the known answer in tests. Rationale: a correction
   derived from the local dynamics should be more faithful to the biology
   than one derived from positions or expression.
13a. **Axes from positions, orientation from operators** (Tom, 2026-10-18;
   after O1i showed operator statistics have flat spectra and do not define
   a transferable subspace, while each copy's position covariance does).
   An estimator may use **each copy's own cell covariance** (its principal
   axes and their variance order, computed within that copy only) to
   choose *which* directions are aligned. Everything those axes leave open
   — the rotation within a block of axes of similar variance, and every
   axis sign — must come from operators and drift only. Still excluded from
   the estimate: cross-copy pairing, labels, gene loadings and any other
   use of positions. Claims are stated accordingly: "positions locate the
   axes; operators orient them". Applies from O1j on; rule 13 otherwise
   unchanged.
14. **Estimator search order and stopping rule** (Tom, 2026-10-15): test
   estimators in order of predicted success until one works, on the same
   known-answer tests and pass rule as O1 (works at n if `err_R ≤ 0.10`;
   `err_F ≤ 0.25`): **O1** (matching family E0–E4) → **O1b** (E5 field
   equivariance) → **O1c** (E7 joint diagonalisation + commutator) →
   **O1d** (E6 shared energy with learnable Q) → **O1e** (E8 eigen-subspace
   matching with fused Gromov–Wasserstein pairing). **Amended 2026-10-15
   (Tom), before any result:** the steps may run **in parallel** (separate
   executor worktrees, at most three at once). Selection is unchanged in
   effect: the first estimator **in the fixed order** O1 (E0…E4), O1f (E9; inserted 2026-10-16 after O1's
   diagnostic, before any O1b–O1e result), O1b (E5), O1g (E5b; inserted 2026-10-16 after O1b, before any
   O1c–O1e result), O1i (E5r; inserted 2026-10-18 after O1g/O1h, before any
   O1c–O1e result), O1j (E5p; inserted 2026-10-18 after O1i, before any
   O1c–O1e result),
   O1c (E7), O1d (E6), O1e (E8) that works is the one used, whatever order
   the results arrive in. It goes to **O2** (real pairs, or larger data if
   it needs more cells than a pair has).
   **Test F criterion from O1g onward** (Tom, 2026-10-16): Test F is
   judged functionally — `mix ≥ 0.80`, operator frame ratio ≤ 1.5 and
   silhouette ≥ 0.9 × no-rotation — instead of distance to `Q_load`
   (reported only). Test R (distance ≤ 0.10 to the exact `Q*`) is
   unchanged. Not retroactive: O1–O1f keep their recorded outcomes.

## Known risks (pinned so task writers see them)

- **Identifiability.** Count-level batch distortion is nonlinear and bends
  the learned energy. Operators alone cannot separate "chemistry bent the
  well" from "biology changed the well." T01 bounds the chemistry effect
  on a pair with no biology; T04 is the only task with both a known batch
  and a known biological difference, so it runs before T03.
- **Correspondence.** Procrustes needs matched pairs. Cell-type labels
  are the weak link; spectral-similarity matching via optimal transport
  is the reserve.
- **Tangent rank.** The restricted spectra need local rank ≥ `r = 10`;
  every K task gates d = 30 layers on the fraction of cells meeting it.
- **Estimator reproducibility.** Untested before K0. A difference between
  datasets means nothing until two halves of one run agree.
- **Local freedom.** Locally varying `A` can bend anything into anything.
  T06 must report residual removed beyond the global `A`; a large value is
  a warning, not a success.

## Task ladder

**Redesign adopted 2026-10-15 (Tom):** positional layer chosen by K0 among
PCA-30, FA-30 and FA-10 (rule 1); tangent-restricted spectra for d = 30
(rule 3a); K0 → K1 → K2 kill-test ladder (rule 10). T02 is blocked by K2.

**2026-10-15, after K0b (Tom):** K0 and K0b were inconclusive for PCA-30.
PCA-30 is the layer in use as Tom's design choice, not a test pass (rule 1).
K1 and K2 are judged with the refit yardstick (rule 12). K0c gives FA-30
the same fit-noise check before K1 is issued.

**Status of T01–T01e:** none of them is a verdict on the axiom. All were
inconclusive or stopped under a design whose positional layer (FA-30, full
spectra) proved unstable; their between-chemistry ratios are recorded
history, not evidence for or against the axiom. Only K1 and K2, on a layer
that K0 shows to be reproducible, can decide it.

| ID   | Purpose | Status |
|------|---------|--------|
| T00  | Environment, scJDO2/scJDO installs, data staging | ACCEPT (2026-10-14) |
| T00b | Positional layer: PCA vs FA vs LDVAE tangent-rank on PBMC pair; pick FA unless it saturates | ACCEPT (2026-10-14) — FA by rule; held, see reports/T00b/REVIEW.md |
| T01  | Kill test: eigenvalue spectra + basin_index agree across v2/v3 chemistry before alignment? | ACCEPT (2026-10-14) — outcome INCONCLUSIVE (rule 3: operators do not separate cell types); T02 stays blocked |
| T00c | Positional layer redone on genes detected in both chemistries; adds a rank floor | ACCEPT (2026-10-14) — FA (ranks 11/14) |
| T01c | Kill test rerun on the T00c layer (T01 design unchanged) | ACCEPT (2026-10-14) — outcome INCONCLUSIVE (rule 3: basin_sep_v2 1.22); T02 stays blocked |
| T01x | *Exploratory:* is the T01c v2/v3 spectral difference one global scale? (reuses T01c operators; cannot change T01c's outcome) | REJECT (2026-10-15) — max_sep_change computed off-definition; planner's validity threshold too tight; see reports/T01x/REVIEW.md |
| T01y | *Exploratory:* T01x reissued with relative validity checks against T01c's reported values | ACCEPT (2026-10-15) — "shape differs" (spectral_ratio_s 15.6) |
| T01d | Kill test at matched depth: v3 thinned to v2's per-type depth; T00c genes/FA and T01c test unchanged | ACCEPT (2026-10-15) — outcome INCONCLUSIVE (rule 3: basin_sep_v2 1.72); T02 stays blocked |
| T01e | Same-chemistry control: pbmc_1k_v3 vs pbmc_10k_v3 (count- and depth-matched), T01d pipeline and rule | ACCEPT (2026-10-15) — stopped, INCONCLUSIVE (rule 0: FA ranks 8/9 < 10) |
| K0   | Split-half control on pbmc_10k_v3 for PCA-30, FA-30 and FA-10: which layer gives reproducible operators? Selects the layer | ACCEPT (2026-10-15) — pca30 inconclusive (1.72/1.79), fa30 fails (3.02), fa10 fails (7.15); selected_layer none → Tom |
| K0b  | PCA-30 only: fit-to-fit noise (5 seeds per half) and K0's split-half test on seed-averaged operators | ACCEPT (2026-10-15) — inconclusive (avg spectral 1.81; basin 1.49, class 1.47); H1 vs H2 ≤ refit noise (0.85/0.98) → Tom |
| K0c  | FA-30 only: the K0b design (fit-to-fit noise, seed-averaged split-half test) | ACCEPT (2026-10-15) — avg fails (spectral 4.39); halves within refit noise (1.07/0.78); larger, less stable fit noise than PCA-30 |
| K0d  | FA-10 only: the K0b design (fit-to-fit noise, seed-averaged split-half test), for completeness | ACCEPT (2026-10-15) — avg fails (spectral 4.09); halves within refit noise (0.93/0.53); fit noise 7.6–7.9×, seed-specific solutions |
| K1   | Same chemistry, two runs: pbmc_1k_v3 vs pbmc_10k_v3 subsample (count/depth matched), PCA-30, refit yardstick | ACCEPT (2026-10-15, rev 1) — inconclusive: spectral 1.529 (> 1.5), basin 0.87, class 1.09 → Tom |
| K2   | Chemistry: pbmc_1k_v2 vs pbmc_1k_v3 (count/depth matched), PCA-30, refit yardstick; ratios also vs K1 | ACCEPT (2026-10-15) — **survives**: spectral 1.21, basin 0.91, class 0.71 (refit); ≤ K1's ratios → T02 proposal to Tom |
| T02  | Global orthogonal `A` on the K2 pair (direct Frobenius estimator on type-mean operators); frame agreement vs refit null; before/after UMAPs | ACCEPT (2026-10-15) — stopped, inconclusive (rule 0: optimum reached by 1 of 20 starts); frame 1.25, held-out 2.53 |
| T02b | T02 with fine-label groups (≥ 6) and a gene-loading rotation as start/comparator | ACCEPT (2026-10-15) — stopped, inconclusive (5 groups < 6) |
| T02c | Part 1: gene-anchored rotation tested by big-cluster centroid operators (top 10 PCs); Part 2: operator-fitted 10-d rotation vs gene one | ACCEPT (2026-10-15) — **frame agreement** (1.09); Part 2: operators validate, don't estimate; mix 0.11 → 0.90 |
| T02d | T02c Part 1 with basin-floor operators (30 lowest-energy cells per cluster); written before T02c's result | ACCEPT (2026-10-15) — frame agreement (1.20; centroid 1.09) |
| T04  | Hypothesis test: PBMC halves + known gene-level batch + known interferon response in half of B's monocytes; 6 methods scored on batch removal, biology retention, operator flag | ACCEPT (2026-10-20) — not_supported: E5 rotation fails under a gene-level batch (mix 0.07, retention 0.68); operator flag works for no method; rigid gene rotation alone keeps biology (0.91 vs Harmony/scVI 0.63) |
| T04c | Diagnostic: within-batch operator contrast (B responders vs non-responders; S0 control) from T04's cached fits | issued (2026-10-20; Manus cloud, same sandbox) |
| T04b | Sweep: biology effect size × batch strength; rigid correction vs Harmony vs scVI | planned (after T04c) |
| T03  | Two-lab dataset: carry `A` to `z`; residual = biology check | blocked by T04 |
| T05  | Benchmark vs Harmony / scVI on scib metrics + preservation of a known biological difference | blocked by T03 |
| T05a | Pilot benchmark on the K2 pair: joint PCA, centring, Harmony, scVI, gene-anchored rotation; mixing, silhouette, operator frame agreement | ACCEPT (2026-10-15) — Harmony best (mix 1.00, sil 0.447, frame 1.35); ours not competitive on this pair |
| C0   | Cloud bootstrap: rebuild env + data in a Manus sandbox; content fingerprints must match the Mac | ACCEPT (2026-10-18) — 8/8 fingerprints, drift repro 2e-6; latent tolerance corrected to relative (planner error; PCA float32 round-off 5e-5) |
| O1   | Operator-estimated rotation on known answers (known Q*, known gene frame); five estimators (baseline, operator ICP, + drift, + seed averaging, + gap weighting); selects estimator and minimum cell count | ACCEPT (2026-10-16) — none works (errors 0.63–1.0); diagnostic: L(Q*) ≪ L(random) → search failure, not missing information |
| O1f  | Targeted fix: seed-averaged group operators + global rotation search (reuses O1 arrays) | ACCEPT (2026-10-16) — fails; estimates ≈ Q*·diag(±1): curvature alone cannot fix axis signs → flow needed (O1b) |
| O1b  | E5 field equivariance (no correspondences); O1's known-answer tests and pass rule | ACCEPT (2026-10-16) — fails (err_R 0.18 at 2759, from a start at the truth; err_F 0.58); Test F view mix 0.81, frame 1.22 |
| O1g  | E5 with 20-seed fields and clean starts; Test F functional | ACCEPT (2026-10-18) — fails (err_R 0.60 at 2759, 0.66 at 2000); Test F passes functionally (mix 0.80, frame 1.24); more fits did not help → flat objective |
| O1h  | Diagnostic: E5b top-k axis error on one cached Test R case (n 2759, r 0); saves all start rotations | ACCEPT (2026-10-18) — partial: top-3 axes err 0.05, top-5 0.12, top-10 0.21, all 30 0.60 → information fades with variance |
| O1i  | E5r: rotation restricted to the operator-defined top-k subspace (k = 5 gated; 3, 8 reported), complement left unchanged; screened sign starts | ACCEPT (2026-10-18) — fails (top-5 err 0.93; Test F mix 0.26): operator-defined subspace does not transfer (flat spectra); position-covariance subspace would |
| O1j  | E5p: axes from each copy's position covariance (rule 13a); top-5 block rotation + remaining axis signs from operators | issued (2026-10-18) |
| O1c  | E7 joint diagonalisation (energy-decile operators) + drift-skew axis signs + commutator refinement | pre-registered (parallel allowed) |
| O1d  | E6 one shared energy with a learnable rotation (ESM) | pre-registered (parallel allowed) |
| O1e  | E8 fused Gromov–Wasserstein initial matching + E3 ICP | pre-registered (parallel allowed) |
| O2   | Selected operator estimator on real technical batches: run-to-run (K1 pair), chemistry (K2 pair); Q from operators → latent; UMAP + mix/sil/frame | pre-registered (issue after O1-series selection) |
| T05b | Benchmark rerun with the operator-estimated rotation as "ours" | planned |
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

- Executor runs on Tom's Mac in the project `.venv`, **or** (amended
  2026-10-15, Tom) in a **Manus cloud sandbox** once C0 has shown the cloud
  rebuild matches the Mac. The datasets are public 10x PBMC data and may be
  processed in the cloud; nothing private leaves the Mac. A task runs
  entirely on one platform; cloud and Mac numbers are not mixed in one
  comparison unless C0 shows the drift fit reproduces across platforms.
- Every report commits `requirements.lock`; checker flags diffs against the
  previous accepted task.
- Git: the local executor folder is a clone of this repo and pushes task
  branches back here. Cloud executors clone the private GitHub copy
  `manarai/scopalign-coord` and push task branches there; the planner pushes
  `main` there when issuing and fetches cloud branches for review.
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
