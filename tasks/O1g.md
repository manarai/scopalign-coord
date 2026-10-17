# O1g — E5 with less noisy fields and clean starts; Test F judged functionally

status: pre-registered (approved by Tom 2026-10-16)
blocked_by: none (O1b ACCEPT)
seed: 20261014
planner_commit: 20d275f
kind: follow-up (targeted fix of E5, Tom's standing instruction)
motivated_by: reports/O1b/REVIEW.md — E5 aligns cells on Test F (mix 0.81, frame 1.22) but its Test R minimum sits ~0.18 from the truth (field-fit noise), and O1b's starts overlapped the true rotation

## Purpose

Keep E5 (field equivariance, flow + curvature, every cell, no pairing) and
fix the two problems O1b exposed:

1. **Clean starts:** start seeds that cannot coincide with `Q*`.
2. **Less noisy fields:** ensemble of **20** drift fits per copy instead of
   5 (fit noise in the ensemble mean drops by about half).

Test R is judged as before (distance to the exact `Q*`). **Test F is judged
functionally** (Tom, 2026-10-16): the estimate must align the two halves —
mixing, operator agreement and kept cell types — rather than match the
gene-loading rotation `Q_load`, which need not be the operators' own frame
map in low-variance directions. The estimate is applied to the latent space
as in PROJECT.md rule 13.

## Scope

### In scope
- **Everything as `tasks/O1b.md`** (data, Test R / Test F constructions,
  `Q*` = `ortho_group.rvs(30, random_state=20261031 + r)`, drift config,
  fields as torch functions of the frozen energy networks, E5 loss `D + H`
  with `H` on 256-cell subsets from `default_rng(20261040 + r)`, LBFGS
  settings), with these changes only:
  - **Cell counts:** `n ∈ (2000, 2759)`, replicates `r ∈ (0, 1, 2)`.
  - **Field ensemble:** 20 seeds per copy, `SEEDS20 = (20261014,
    20261016, 20261017, …, 20261034)` (20261014 plus 20261016–20261034).
    O1b's cached fits for the first five seeds may be reused (they are
    identical); the other 15 are fitted here. `ū`, `J̄` = means over all 20.
  - **Starts:** 20 — `Q₀ = I` and `ortho_group.rvs(30, random_state=20263000 + i)`
    for `i = 1…19` (none can equal `Q*`; check and record
    `start_overlap = max_i,r 1[‖Q₀,i − Q*_r‖_F < 1e-8]`, must be 0).
    Starts may run in parallel worker processes (deterministic, identical
    results).
- **Test R diagnostics per case:** `err`, `L(Q̂)`, `L(Q*)`, `best_start`.
- **Test F functional measures** (H1 rotated by `Q̂` vs H2, as O1b's views):
  `mix_F` (30 nearest other cells from the other half, ÷ 0.5), `sil_F`
  (K0 tested types), `frame_ratio_F` (full-30-d centroid-operator residual
  under `Q̂` ÷ T02c's refit null, `default_rng(20261022)`, seed-20261014 vs
  20261016–19 fits), and the same three for **no rotation** (`*_noA`).
  `err_F` (distance to `Q_load`) is reported, not gated in the rule.
- **Views:** O1b's three Test F panels (no rotation, `Q̂`, `Q_load`) with
  mix, sil, frame ratio in the titles.

### Out of scope
- Any change to E5's objective, loss weights, subset sizes, start count or
  thresholds beyond the above; pairing, labels beyond K0 tested types for
  scoring, gene loadings or positions in the estimate.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_E5b_2000` | median over 3 replicates, n = 2000 | `python reports/O1g/scripts/o1g_run.py --stage testR` |
| `err_R_E5b_2759` | n = 2759 | same |
| `frac_hat_below_star` | fraction of the 6 Test R cases with `L(Q̂) ≤ L(Q*)` | same |
| `start_overlap` | must be 0 | same |
| `mix_F` | Test F mixing under `Q̂` | `--stage testF` |
| `mix_F_noA` | Test F mixing, no rotation | same |
| `sil_F` | Test F silhouette under `Q̂` | same |
| `sil_F_noA` | Test F silhouette, no rotation | same |
| `frame_ratio_F` | Test F frame ratio under `Q̂` | same |
| `err_F` | distance to `Q_load` (reported, not in the rule) | same |
| `symmetry_max` | max `median_asym` over all fits | same |
| `latent_check_max` | as O1 | same |
| `e5b_outcome` | `"works"` or `"fails"` | `--stage report` |
| `min_n_e5b` | smallest n with `err_R_E5b_<n> ≤ 0.10`; 0 if none | same |

### Decision rule
0. Validity: `start_overlap ≠ 0`, `symmetry_max > 1e-4` or
   `latent_check_max > 1e-4` → **REVISE**.
1. **Test R passes** if some `err_R_E5b_<n> ≤ 0.10`; `min_n_e5b` = smallest
   such n.
2. **Test F passes** if `mix_F ≥ 0.80` **and** `frame_ratio_F ≤ 1.5`
   **and** `sil_F ≥ 0.9 × sil_F_noA`.
3. E5b **works** if Test R and Test F both pass; otherwise **fails**.

Reading (not the outcome): `frac_hat_below_star` near 1 with large `err`
means the objective's minimum is away from the truth (noise-limited);
near 0 means the search missed the truth's basin.

### What would change my mind
If quadrupling the field ensemble does not bring Test R below 0.10 at 2759
cells, the remaining error is not drift-fit noise but a systematic
difference between operator fields fitted on different cells; then a
rotation can be estimated to ~0.15–0.2 at best at these cell counts.

## Procedure

1. `o1g_run.py --stage fits` — 15 new seeds per copy for 6 cases + Test F
   (fits in parallel workers; cached); reuse O1b's five cached seeds if
   present (record sha256 of reused files).
2. `o1g_run.py --stage testR` — E5 with 20-seed fields; `testR.csv` (n,
   replicate, err, L_hat, L_star, best_start).
3. `o1g_run.py --stage testF` — E5; functional measures; `testF.csv`,
   `Q_hat_F.npy`; views → `views.csv`, `fig_umap.png`.
4. `o1g_run.py --stage report` — `fig_recovery.png` (err_R vs n with O1b's
   E5), `report.json`, checksums, lockfile.

`o1g_run.py` must accept `--stage all`. Run detached with stdin from
`/dev/null` (MANUS.md "Long runs").

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `views.csv`
- `fig_umap.png`
- `fig_recovery.png`
- `scripts/o1g_run.py`

## Checker expectations

- `gated` keys: exactly the fourteen above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
