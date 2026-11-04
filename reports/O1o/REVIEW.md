# O1o — Review

Reviewer: planner (Claude), 2026-11-04. Branch `task/O1o` @ b313ad6, Mac worktree `~/Documents/scopalign-manus-O1o`.

## Checker
`python checker/check_report.py tasks/O1o.md reports/O1o/report.json --prev reports/O1n/report.json` → **PASS** (0 fail, 0 warn).

## Recomputation
Light recompute from `pairs.csv`: all aging medians, `n_pass_A` and the pass rule recomputed exactly. Not recomputed: fits, LBFGS, Harmony.

| pair | noA | skew | **E5s** | gene | Harmony | sil noA → E5s | err_gene | passes |
|---|---|---|---|---|---|---|---|---|
| K2 (chemistry) | 0.526 | 0.514 | **0.885** | 0.880 | 0.986 | 0.568 → 0.574 | 0.085 | **yes** |
| F30–F36 | 0.047 | 0.020 | 0.070 | 0.120 | 0.785 | 0.180 → 0.072 | 0.975 | no |
| M28–M36 | 0.072 | 0.032 | 0.178 | 0.265 | 0.849 | 0.211 → 0.154 | 0.978 | no |
| F66–F77 | 0.047 | 0.042 | 0.045 | 0.129 | 0.668 | 0.227 → 0.111 | 1.078 | no |
| M60–M66 | 0.051 | 0.107 | 0.131 | 0.236 | 0.824 | 0.176 → 0.275 | 0.748 | no |

symmetry_max 1.9e-7. Deviation (declared): local caches and in-process executor. No effect on values.

## Decision rule outcome
Validity: pass. K2 passes (0.885 ≥ 0.80; sil kept; beats skew by 0.37). Aging: 0 of 4 pass.
→ **chemistry_only**.

## Verdict
**ACCEPT** — outcome **chemistry_only**.
- On a real technical batch (10x v2 vs v3, depth-matched), operator-estimated rotation in 10 PCs aligns the two datasets as well as the gene-defined rotation: mix 0.53 → 0.89 (gene 0.88), with cell types kept. It agrees with the gene frame (err 0.085), uses no genes, labels or pairing, and beats positions (0.51). Harmony mixes better (0.99).
- On donor pairs, no rotation works: the gene rotation reaches only 0.12–0.27 and Harmony 0.67–0.85. The operators find a rotation with a lower loss than the gene answer in every pair, but it misplaces cell types (sil falls in 3 of 4 pairs; UMAP of F30–F36 shows NK and monocytes laid onto different types).

So donor batches here are not rotations of each copy's own top-10 PCA frame. Fitting the fields better will not fix that.

## Quarantined (not evidence)
- Planner seeds:
  - Each donor's own top-10 PCA subspace may span different genes, so a 10×10 rotation cannot map one onto the other. A shared gene subspace (e.g. a joint PCA followed by a per-donor operator rotation) would separate "wrong subspace" from "not a rotation".
  - In donor pairs, the loss prefers a wrong rotation (L 1.2–1.4 < L_gene 1.8–2.3): operator fields differ by donor beyond rotation, consistent with donor effects that vary by cell type (planner's exploratory aging run).
  - The K2 no-rotation mix (0.53) differs from T02c's 0.64: likely a different neighbour set or cell set, not checked.
