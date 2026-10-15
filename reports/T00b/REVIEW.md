# T00b review

reviewed_commit: 0bc9c53 (task/T00b; 1796a4e excluded, see Checker)
date: 2026-10-14
verdict: ACCEPT

## Checker

- Executor's DONE line reported `checker=FAIL`. Cause: two planner bugs —
  the report schema rejected the task id `T00b`, and `parse_known_failures`
  stopped at the first `]` inside a parametrised test id. Fixed by the
  planner on main (262ed7b).
- Commit 1796a4e on `task/T00b` edits `checker/check_report.py` and
  `templates/report.schema.json`, which MANUS.md forbids the executor to
  touch. It was made after the DONE line, at a follow-up prompt, and does
  not touch any output. **Not merged**: main's checker and schema were kept.
  Recorded as a protocol breach; future checker failures come to the planner.
- Main's checker on the report with `--prev reports/T00/report.json` →
  PASS (0 fail, 1 warn: lockfile anndata 0.10.9 → 0.12.19, declared below).

## Recomputation

Planner's own script (independent HVG, normalisation, refits, tangent rank),
not the executor's:

| gated | report | recomputed | match |
|---|---|---|---|
| `rank_pca_v2` | 21 | 21 | yes |
| `rank_pca_v3` | 19 | 19 | yes |
| `rank_fa_v2` | 6 | 6 | yes |
| `rank_fa_v3` | 6 | 6 | yes |
| `rank_ldvae_v2` | 19 | 19 (seeds 20/19/19) | yes |
| `rank_ldvae_v3` | 13 | 13 (seeds 13/13/13; executor 13/13/12) | yes |
| `frac_fa_at_ceiling` | 0.0 | 0.0 | yes |

Rank is unchanged whether or not the cell itself is counted among its 40
neighbours. The executor's regeneration command reproduces both FA latent
deliverables bit-for-bit.

## Decision rule outcome

`frac_fa_at_ceiling = 0.0 ≤ 0.20` → **FA is the positional layer**.

## Deviations

- anndata upgraded 0.10.9 → 0.12.19 so `scvi.model.LinearSCVI` imports.
  Declared. FA/PCA ranks recompute identically under it; cannot move a
  gated number. Accepted. (The executor also adds `1e-6·I` to each local
  covariance before eigendecomposition; numerically immaterial at a 0.05
  relative floor; recomputation without it matches.)

## Verdict

**ACCEPT.** The executor ran the task as written and every gated number
recomputes.

## Planner notes — hold before T01 (not findings)

These are planner checks on the deliverables, not evidence. They concern the
pre-registered procedure itself (planner's design), not the executor's run.

1. **Chemistry-specific HVGs.** HVGs are selected on the concatenated
   counts. 135 of the 2000 are all-zero in v2 QC cells (9 in v3), e.g.
   ANKRD65, AJAP1, TRNP1. v3's higher sensitivity puts v3-only genes into
   the shared gene list, so the "same genes" are not measured the same way
   in both chemistries.
2. **Degenerate FA noise model.** FA assigns noise variance ≈ 1e-12 to 239
   genes in v2 and 49 in v3. Genes with near-zero variance then get very
   large weight in the factor scores; v2's FA log-likelihood is +1.75e6
   versus −2.6e5 for v3. FA's median tangent rank of 6 (PCA 21, LDVAE
   13–19) may be this artefact rather than a property of the cells.
3. **FA v2 is not uniquely determined.** Refitting with the genes in a
   different column order (both fits converge, ~43 iterations) gives a v2
   latent that differs from the deliverable by ~6% beyond any linear map
   (one sign flip, weakest factor r = 0.956); v3 is identical. T01 reads
   the deliverable files, so T01 itself stays reproducible, but its
   operators would rest on one of several FA optima.

Consequence: the decision rule picked FA from ranks that may reflect points
1–2. T01 should not run on these latents until a planner task re-checks the
choice with a gene filter that requires detection in both chemistries.
Proposed to Tom as T00c.

## Quarantined (not evidence)

(The report's `observations` and `deviations_requested` are empty.)

## Decision (Tom, 2026-10-14)

T00c not written. Proceed to T01 on the T00b FA latents as pre-registered.
The planner notes above stand as recorded caveats for interpreting T01.
