# O1p — Anchored per-group alignment: operator-matched clusters, each with its own shift and operator rotation

status: pre-registered (approved by Tom 2026-11-04)
blocked_by: none (O1o, K4 ACCEPT)
seed: 20261014
planner_commit: 06b4dae
kind: follow-up (O1o: one global rotation corrects the chemistry pair but not donor pairs; K4: operator fingerprints match clusters across batches)
motivated_by: Tom 2026-11-04 — "one rotation per matched cell group, combining the K4 anchoring with this method". O1o: donor effects differ by cell type, so no single rotation (operator or gene) works, while Harmony (per-cluster correction) does.

## Purpose

Correct each cell group separately, using only operators and each copy's
own cells:

1. **Anchor:** cluster each copy on its own, and match clusters across
   copies by operator fingerprints (K4). No labels or genes.
2. **Shift:** move each matched cluster's centroid onto its partner's.
3. **Rotate:** within each matched pair, estimate a rotation with E5s from
   the operators of those cells.

A **shift-only** version (steps 1–2) is the key control. If it does as
well, the per-group operator rotation adds nothing. An **oracle** version
(clusters matched by labels) separates matching errors from rotation
errors.

## Scope

### In scope
- **Pairs, cells, latents, fits:** O1o's, read-only from
  `~/Documents/scopalign-manus-O1o` (cached inputs and the 50 fits,
  sha256 must match O1o's `fit_symmetry.csv`; **no new fits**). Pairs:
  K2 (chemistry) and the four aging donor pairs. d = 10. Copy 1 is mapped
  into copy 2's frame.
- **Clusters (per copy):** K4's Leiden settings on the copy's 10-d latent
  (`n_neighbors=15`, `resolution=1.0`, igraph, `n_iterations=2`,
  `random_state=20261014`).
  - Clusters with ≥ 30 cells are kept. Cells of smaller clusters are
    *unassigned*.
- **Anchoring (operators):** K4's fingerprint with d = 10.
  - Per cluster and seed: `eigvalsh(centroid_operator(z, M_t,
    n_cells_centroid=30))`, sorted descending.
  - Fingerprint = mean over the 5 seeds, standardised per coordinate over
    both copies' clusters pooled; distance = Euclidean.
  - Matching: Hungarian (`linear_sum_assignment`), rectangular allowed.
    Unmatched clusters are *unassigned*.
- **Oracle anchoring (reading control):** the same clusters, matched by
  Hungarian on `1 − (fraction of shared majority coarse label)`, i.e.
  clusters with the same majority label are matched first.
- **Per matched pair g** (clusters `C1_g`, `C2_g`; centroids `c1_g`,
  `c2_g`):
  - **shift:** `x' = x − c1_g + c2_g`.
  - **anchor_op:** `x' = Q_g (x − c1_g) + c2_g`. `Q_g` comes from E5s
    (O1n/O1o exactly: axis-sign starts from the two clusters' own
    covariances, all 1024 sign patterns, refine the 8 lowest). The loss
    uses only the pair's cells, with fields from the copies' ensembles
    evaluated at absolute coordinates: the models are wrapped as
    `E_X(z + c1_g)` and `E_Y(z + c2_g)` on centred cells. `H` subsets:
    `min(256, n)` from `default_rng(20261040 + g)`.
  - Pairs where either cluster has < 100 cells get `Q_g = I` (shift only);
    the number is recorded.
- **Unassigned cells:** mapped with the pair's global O1o E5s rotation
  (read from O1o's `Q_hats.npz`), identically in every method.
- **Methods scored:**
  - `noA`
  - `global` (O1o E5s, recomputed embedding)
  - `shift`
  - `anchor_op`
  - `oracle_shift`
  - `oracle_op`
  - `harmony` (O1o's embedding, read)
- **Scores:**
  - `mix` and `sil`, as O1o.
  - `acc_anchor`: the fraction of operator-matched cluster pairs whose
    majority coarse labels agree (K4).
  - `n_groups`, `n_rotated`, `frac_unassigned`.
- **Views:** UMAPs for K2 and F30–F36: noA / shift / anchor_op /
  harmony, coloured by copy and by coarse type, with `mix` and `sil` in
  the titles.

### Out of scope
- New fits; other clusterings, resolutions or size thresholds; non-rigid
  maps within a group; label use anywhere except the oracle and scoring;
  cross-age pairs; age-retention scoring.

## Pre-registration

### Gated numbers
Aging values (`_A`) are medians over the 4 donor pairs.

| name | definition | computed by |
|------|-----------|-------------|
| `inputs_match` | 1 if O1o inputs and the 50 fit sha256 match, else 0 | `python reports/O1p/scripts/o1p_run.py --stage check` |
| `acc_anchor_K2` | K2 anchor accuracy | `--stage anchor` |
| `acc_anchor_A` | aging median anchor accuracy | same |
| `mix_shift_K2` | K2 mix, shift | `--stage estimate` |
| `mix_anchor_K2` | K2 mix, anchor_op | same |
| `sil_anchor_K2` | K2 sil, anchor_op | same |
| `mix_noA_A` | aging median mix, noA | same |
| `mix_shift_A` | aging median mix, shift | same |
| `mix_anchor_A` | aging median mix, anchor_op | same |
| `mix_oracle_A` | aging median mix, oracle_op | same |
| `mix_harmony_A` | aging median mix, Harmony | same |
| `sil_noA_A` | aging median sil, noA | same |
| `sil_shift_A` | aging median sil, shift | same |
| `sil_anchor_A` | aging median sil, anchor_op | same |
| `sil_oracle_A` | aging median sil, oracle_op | same |
| `n_pass_A` | aging pairs passing rule 1 (0–4) | `--stage report` |
| `n_shift_pass_A` | aging pairs where shift alone meets rule 1's first two conditions | same |
| `o1p_outcome` | `"anchored_works"`, `"shift_suffices"` or `"fails"` | same |

### Decision rule
0. **Validity → REVISE:** `inputs_match = 0`.
1. An aging pair **passes** if all three hold:
   - `mix_anchor ≥ 0.80`;
   - `sil_anchor ≥ 0.9 × sil_noA`;
   - the operators add over shift: `mix_anchor ≥ mix_shift + 0.05` **or**
     `sil_anchor ≥ sil_shift + 0.05`.
2. Outcome:
   - **anchored_works** if `n_pass_A ≥ 3`;
   - **shift_suffices** if `n_pass_A < 3` and `n_shift_pass_A ≥ 3` (the
     anchors work, the operator rotations add nothing);
   - **fails** otherwise.

Reading (not the outcome):
- `oracle_op` vs `anchor_op`: how much matching errors cost.
- vs `harmony`: the remaining gap.
- K2: whether per-group correction keeps O1o's global result.

### What would change my mind
If operator-matched shifts plus per-group operator rotations cannot mix
donor pairs (mix < 0.80) even with oracle matching, then donor effects
are not rigid within cell groups either. scOPAlign's correction would then
be limited to technical batches, with operators used to validate rather
than to correct for biological batches.

## Procedure

1. `--stage check` — O1o inputs and fits; `inputs.json`.
2. `--stage anchor` — clusters and matches (operator and oracle);
   `clusters.csv`, `matches.csv`.
3. `--stage estimate` — per-group shifts and rotations; `groups.csv`
   (pair, g, n1, n2, rotated, L_hat, err_vs_I); `pairs.csv` (pair, method,
   mix, sil, n_groups, frac_unassigned); `Q_groups.npz`.
4. `--stage report` — `fig_mix.png`, `fig_umap_K2.png`,
   `fig_umap_F30F36.png`, `report.json`, checksums, lockfile.

`--stage all`; detached, stdin `/dev/null`, resumable, checkpoint per pair.
Mac (**Manus Desktop, not the cloud**), worktree
`~/Documents/scopalign-manus-O1p`, branch `task/O1p` from `main`; commit
and push after each stage. O1o's worktree and the aging file are
read-only.

## Deliverables

- `report.json`
- `requirements.lock`
- `clusters.csv`
- `matches.csv`
- `groups.csv`
- `pairs.csv`
- `fig_mix.png`
- `fig_umap_K2.png`
- `fig_umap_F30F36.png`
- `scripts/o1p_run.py`

`Q_groups.npz` and `inputs.json` go in `artifacts`.

## Checker expectations

- `gated` keys: exactly the eighteen above (or a subset with `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
