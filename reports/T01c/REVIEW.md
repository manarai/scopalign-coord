# T01c review

reviewed_commit: 13fb1bc (task/T01c)
date: 2026-10-14
verdict: ACCEPT (run) — decision-rule outcome: INCONCLUSIVE (rule 3)

## Checker

`check_report.py tasks/T01c.md reports/T01c/report.json --prev reports/T00c/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/T01c/`.

## Recomputation

Planner's own script (same as for T01, inputs switched to the T00c FA
latents): drifts refit with the pre-registered config; Jacobians
bit-identical to `ops_v2.npz` / `ops_v3.npz` (max |ΔJ| = 0.0); input sha256
446a5414… / 691d66b9… match the T00c deliverables.

| gated | report | recomputed | match |
|---|---|---|---|
| `n_cell_types_tested` | 4 | 4 | yes |
| `symmetry_v2` | 8.69e-08 | 8.69e-08 | yes |
| `symmetry_v3` | 1.24e-07 | 1.24e-07 | yes |
| `spectral_w2_pooled` | 506.25 | 506.25 | yes |
| `spectral_w2_null` | 9.436 | 9.436 | yes |
| `spectral_ratio` | 53.65 | 53.65 | yes |
| `basin_w1_pooled` | 5.163 | 5.163 | yes |
| `basin_w1_null` | 2.243 | 2.243 | yes |
| `basin_ratio` | 2.302 | 2.302 | yes |
| `class_tv_pooled` | 0.3618 | 0.3618 | yes |
| `class_tv_null` | 0.05559 | 0.05559 | yes |
| `class_tv_ratio` | 6.508 | 6.508 | yes |
| `spectral_sep_v2` | 2.958 | 2.958 | yes |
| `spectral_sep_v3` | 13.77 | 13.77 | yes |
| `basin_sep_v2` | 1.221 | 1.221 | yes |
| `basin_sep_v3` | 3.530 | 3.530 | yes |
| `basin_index_rho` (secondary) | 0.20 | 0.20 | yes |
| `basin_class_agreement` (secondary) | 0.25 | 0.25 | yes |

Differences ≤ 1e-7 relative (float rounding).

Figures: `fig_spectra.png` (per-type eigenvalue densities v2 vs v3, null
inset with pooled value) and `fig_basin.png` (per-type median basin_index,
identity line, ρ) are as specified.

## Decision rule outcome

1. Symmetry 8.7e-08 / 1.2e-07 ≤ 1e-4 → valid run.
2. 4 types ≥ 3 → continue.
3. `basin_sep_v2` = 1.22 < 2 (the other three separations pass: 2.96,
   13.77, 3.53) → **INCONCLUSIVE.** Rules 4–6 are not reached.

## Deviations

None listed. None found.

## Verdict

**ACCEPT** the run. Outcome per the pre-registered rule: **inconclusive**.
T02 stays blocked.

## Planner reading (not findings)

- The gene filter fixed what it was meant to fix: spectra now separate cell
  types in both chemistries (sep 2.96, 13.77; T01 had 0.73, 1.24).
- The between-chemistry spectral difference is large: `spectral_ratio` 53.7.
  Without the v2 basin guard, rule 4 would have called **axiom fails**.
- Most of it is a global scale: median |λ| is 30.3 in v2 vs 4.8 in v3
  (≈ 6.3×); range −1727…+195 vs −550…+161. v2's wells are uniformly
  much steeper. Because `tol` is 0.1 × each dataset's median |λ|, the
  classes then differ too: v2 has no holders (747 drift / 245 fork), v3 has
  247 holders (377 drift / 489 fork).
- An orthogonal A cannot change eigenvalues, so if chemistry rescales the
  whole spectrum, the axiom as written ("wells are not reshaped") fails even
  if the wells' *shape* agrees. Whether shape agrees after a single
  per-dataset scale is not tested by any pre-registered number.
- `basin_sep_v2` < 2 means v2's basin_index barely separates the four types
  (per-type medians −0.91…−1.91), so v2 adds little on that axis.

## Quarantined (not evidence)

- Observations: input latent sha256 values (confirmed by recomputation).
- `deviations_requested`: none.
