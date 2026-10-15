# T00c review

reviewed_commit: 91a5449 (task/T00c)
date: 2026-10-14
verdict: ACCEPT — decision-rule outcome: FA

## Checker

`check_report.py tasks/T00c.md reports/T00c/report.json --prev reports/T01/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/T00c/`.

## Recomputation

Planner's own script (independent gene filter, HVG, sorting, normalisation,
refits, tangent rank, FA noise count):

| gated | report | recomputed | match |
|---|---|---|---|
| `n_genes_detected_both` | 12890 | 12890 | yes |
| `n_hvg_zero_v2` | 0 | 0 | yes |
| `n_hvg_zero_v3` | 0 | 0 | yes |
| `n_fa_tiny_noise_v2` | 0 | 0 | yes |
| `n_fa_tiny_noise_v3` | 0 | 0 | yes |
| `rank_pca_v2` | 21 | 21 | yes |
| `rank_pca_v3` | 19 | 19 | yes |
| `rank_fa_v2` | 11 | 11 | yes |
| `rank_fa_v3` | 14 | 14 | yes |
| `rank_ldvae_v2` | 19 | 19 (seeds 20/19/19) | yes |
| `rank_ldvae_v3` | 13 | 13 (seeds 13/13/13) | yes |
| `frac_fa_at_ceiling` | 0.0 | 0.0 | yes |
| `frac_ldvae_at_ceiling` | 0.0 | 0.0 | yes |

`genes.txt` matches the planner's list exactly. FA latents reproduce to
≤ 9e-7; LDVAE seed-20261014 latents reproduce bit-for-bit. (LDVAE medians
equal T00b's; confirmed real, not a stale file.)

Figure `fig_rank_dist.png`: per-cell rank histograms, v2 vs v3, one panel
per embedding, as specified.

## Decision rule outcome

1. `n_hvg_zero_v2` = `n_hvg_zero_v3` = 0 → gene filter applied.
2. FA: medians 11 and 14 (both ≥ 10), gap 3 (≤ 5), ceiling 0.0 (≤ 0.20) →
   **admissible → FA is the positional layer for T01c.**

(LDVAE, not needed: medians 19/13, gap 6 → would not be admissible.)

## Deviations

None listed. None found.

## Verdict

**ACCEPT.** T01c reads `reports/T00c/latent_fa_v2.npy` / `latent_fa_v3.npy`.

## Planner reading (not findings)

With genes measured in both chemistries, FA's tangent rank rises from 6
(T00b) to 11/14 and its near-zero noise variances disappear (239/49 → 0/0).
This supports, but does not prove, that T01's uninformative operators came
from the gene list. T01c tests that directly.

## Quarantined (not evidence)

- "The LDVAE histogram pools the three prescribed seed-specific rank arrays
  for each dataset." — describes the figure; consistent with the task.
