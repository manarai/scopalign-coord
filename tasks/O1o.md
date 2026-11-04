# O1o — Real batches at 10 PCs: chemistry pair and donor pairs

status: pre-registered (approved by Tom 2026-11-04)
blocked_by: none (O1n ACCEPT)
seed: 20261014
planner_commit: 1cb8507
kind: follow-up (O1n: operator rotation works on same-run halves at d = 5–10)
motivated_by: Tom 2026-11-04 — "proceed" (test a real batch at 10 PCs). O1n aligned two PCA frames of one run. A real batch also shifts and reshapes cells, so a rotation may not be enough. The planner's exploratory aging run (previews/aging) showed donor effects that differ by cell type.

## Purpose

Apply O1n's E5s at **d = 10** to real batches, with each copy in its own
PCA frame. The questions:

1. Do the operators find a rotation that mixes the two copies while
   keeping cell types apart?
2. Do they beat the positions-only rule?
3. How close do they get to the gene-defined rotation (T02c) and to
   Harmony?

There is no known answer here. The gene rotation `Q_load` is the
reference, and Harmony is the external benchmark.

## Scope

### In scope
- **Pairs:**
  - **K2 — chemistry:** matched v2 (A) vs depth-matched v3-thin (B), K2's
    cells and coarse `cell_type` (as T02c/T05a). Latents
    `reports/K2/latent_pca30_a.npy` / `_b.npy`, first 10 columns. Target
    `Q_load,10 = U Vᵀ` of `svd(W_B[:10] W_A[:10]ᵀ)` (T02b `loadings.npz`).
  - **Aging donors — four within-age, same-sex pairs:** (F30, F36),
    (M28, M36), (F66, F77), (M60, M66), from
    `/Users/terooatt/Documents/ChatGPT/Bio559R_studio0/aging_pbmc_8donors_prepared.h5ad`
    (sha256 `43e39835…dd4`; read-only).
    - Genes: `sc.pp.highly_variable_genes(n_top_genes=2000,
      flavor="seurat", batch_key="sample_id")` on the stored log-normalised
      `X`, over all cells. Record `genes.txt`.
    - Cells: 4000 per donor, `default_rng(20261500 + i).choice(n_d, 4000,
      replace=False)`, with i = the donor's index in sorted donor order.
    - Per donor: `PCA(30, svd_solver="full", random_state=20261014)` on its
      cells × genes (centred on its own mean); `z = scores[:, :10]`,
      `W = components_[:10]`.
    - Pair target: `Q_load,10` = `U Vᵀ` of `svd(W_2 W_1ᵀ)`.
    - Coarse labels from `cell_type`, by prefix: Mono (CD14-MC*,
      Intermed-MC, CD16-MC), NK, CD4 T, CD8 T, B (Naive-B*, Memory-B*,
      ABC, Plasma), DC (cDC, pDC), Mk, RBC (excluded from `sil`).
- **Fits:** O1n's config with input dimension 10 and 5 seeds per copy
  (20261014, 20261016–19); ensemble mean. 2 (K2) + 8 donors (each fitted
  once and shared by its pair) = 10 copies × 5 = **50 fits**.
  `symmetry_check` per fit.
- **Estimators per pair** (copy 1 is rotated into copy 2's frame):
  - `noA`: no rotation.
  - `skew`: O1n's positions-only sign rule.
  - `E5s`: O1n exactly (all 1024 sign patterns, refine the 8 lowest).
  - `gene`: `Q_load,10` (uses gene loadings; reference only).
  - `harmony`: a different embedding. For each pair, stack both copies'
    log matrices (same genes), take a joint `PCA(10, svd_solver="full",
    random_state=20261014)`, then `harmonypy.run_harmony(..., "copy",
    random_state=20261014, max_iter_harmony=20)`.
- **Scores per pair and method:**
  - `mix`: 30 nearest other-copy cells, ÷ 0.5, as O1n.
  - `sil`: coarse types with ≥ 20 cells in both copies.
  - `err_gene = ‖Q̂ − Q_load,10‖_F / √20` for skew and E5s.
  - `L(Q̂)`, `L(Q_load,10)`, `screen_rank_best`.
- **Views:** UMAPs (O1g's settings) for K2 and for pair (F30, F36): noA /
  E5s / gene / harmony, coloured by copy and by coarse type, with `mix`
  and `sil` in the titles.

### Out of scope
- Other d, pairs, donors, cell counts, seeds or starts; cross-age pairs
  (age would be confounded with donor); age-retention scoring.
- Per-group (anchored) rotations; correcting PCs beyond 10.
- Any re-tuning after seeing a score.

## Pre-registration

### Gated numbers
For K2, then for the aging pairs as the median over the 4 pairs (`_A`).

| name | definition | computed by |
|------|-----------|-------------|
| `symmetry_max` | max `median_asym` over 50 fits | `python reports/O1o/scripts/o1o_run.py --stage fits` |
| `mix_noA_K2` | K2 mix, no rotation | `--stage estimate` |
| `mix_E5s_K2` | K2 mix, E5s | same |
| `mix_skew_K2` | K2 mix, skew | same |
| `mix_gene_K2` | K2 mix, gene rotation | same |
| `mix_harmony_K2` | K2 mix, Harmony | same |
| `sil_noA_K2` | K2 sil, no rotation | same |
| `sil_E5s_K2` | K2 sil, E5s | same |
| `err_gene_K2` | K2 `err_gene`, E5s | same |
| `mix_noA_A` | aging median mix, no rotation | same |
| `mix_E5s_A` | aging median mix, E5s | same |
| `mix_skew_A` | aging median mix, skew | same |
| `mix_gene_A` | aging median mix, gene | same |
| `mix_harmony_A` | aging median mix, Harmony | same |
| `sil_noA_A` | aging median sil, no rotation | same |
| `sil_E5s_A` | aging median sil, E5s | same |
| `err_gene_A` | aging median `err_gene`, E5s | same |
| `n_pass_A` | aging pairs passing rule 1 (0–4) | `--stage report` |
| `o1o_outcome` | `"works_real"`, `"chemistry_only"`, `"donor_only"` or `"fails"` | same |

### Decision rule
0. **Validity → REVISE:** `symmetry_max > 1e-4`.
1. A pair **passes** if all three hold:
   - `mix_E5s ≥ 0.80`;
   - `sil_E5s ≥ 0.9 × sil_noA`;
   - `mix_E5s ≥ mix_skew + 0.10` (beats positions).
2. Outcome:
   - **works_real** if K2 passes and `n_pass_A ≥ 3`;
   - **chemistry_only** if K2 passes and `n_pass_A < 3`;
   - **donor_only** if K2 fails and `n_pass_A ≥ 3`;
   - **fails** otherwise.

Reading (not the outcome):
- `mix_E5s` vs `mix_gene`: do operators reach what the genes give?
- vs `mix_harmony`: the remaining gap to a non-rigid method.
- `err_gene` near 0 means operators and genes agree on the frame.
- `L(Q_load) < L(Q̂)` would mean the search missed.

### What would change my mind
If operators align same-run halves (O1n) but no real batch at the same d,
then real batch effects are not rotations even in 10 PCs. A global
operator rotation is then not a batch corrector, and any real use needs
per-group rotations (anchoring) or a non-rigid step.

## Procedure

1. `--stage check` — inputs sha256, `genes.txt`, cell choices, latents,
   `Q_load` per pair; `inputs.json`.
2. `--stage fits` — 50 fits; `fit_symmetry.csv`.
3. `--stage estimate` — 5 pairs; `pairs.csv` (pair, method, mix, sil,
   err_gene, L_hat, L_gene, screen_rank_best); `Q_hats.npz`.
4. `--stage report` — `fig_mix.png` (mix and sil per pair and method, line
   at 0.80), `fig_umap_K2.png`, `fig_umap_F30F36.png`, `report.json`,
   checksums, lockfile.

`--stage all`; detached, stdin `/dev/null`, resumable, checkpoint per pair.
Mac (**Manus Desktop, not the cloud**), worktree
`~/Documents/scopalign-manus-O1o`, branch `task/O1o` from `main`; commit
and push after each stage. Other worktrees and the aging file are
read-only.

## Deliverables

- `report.json`
- `requirements.lock`
- `genes.txt`
- `pairs.csv`
- `fit_symmetry.csv`
- `fig_mix.png`
- `fig_umap_K2.png`
- `fig_umap_F30F36.png`
- `scripts/o1o_run.py`

`Q_hats.npz` and `inputs.json` go in `artifacts`.

## Checker expectations

- `gated` keys: exactly the nineteen above (or a subset with `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
