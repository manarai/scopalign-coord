# T01 review

reviewed_commit: 2673fe9 (task/T01)
date: 2026-10-14
verdict: ACCEPT (run) — decision-rule outcome: INCONCLUSIVE (rule 3)

## Checker

`check_report.py tasks/T01.md reports/T01/report.json --prev reports/T00b/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/T01/`.

## Recomputation

Planner's own script: drifts refit from the committed T00b FA latents with
the pre-registered config, own Hessians/eigenvalues/tol/classification,
own draws on the pre-registered generator and order. Refit Jacobians are
bit-identical to `ops_v2.npz` / `ops_v3.npz` (max |ΔJ| = 0.0). Input
latent sha256 match the T00b deliverables.

| gated | report | recomputed | match |
|---|---|---|---|
| `n_cell_types_tested` | 4 | 4 | yes |
| `symmetry_v2` | 1.156e-07 | 1.156e-07 | yes |
| `symmetry_v3` | 1.118e-07 | 1.118e-07 | yes |
| `spectral_w2_pooled` | 1233.12 | 1233.12 | yes |
| `spectral_w2_null` | 61.871 | 61.871 | yes |
| `spectral_ratio` | 19.93 | 19.93 | yes |
| `basin_w1_pooled` | 3.388 | 3.388 | yes |
| `basin_w1_null` | 0.5659 | 0.5659 | yes |
| `basin_ratio` | 5.988 | 5.988 | yes |
| `class_tv_pooled` | 0.00603 | 0.00603 | yes |
| `class_tv_null` | 0.000581 | 0.000581 | yes |
| `class_tv_ratio` | 10.38 | 10.38 | yes |
| `spectral_sep_v2` | 0.733 | 0.733 | yes |
| `spectral_sep_v3` | 1.237 | 1.237 | yes |
| `basin_sep_v2` | 10.11 | 10.11 | yes |
| `basin_sep_v3` | 1.036 | 1.036 | yes |
| `basin_index_rho` (secondary) | −0.80 | −0.80 | yes |
| `basin_class_agreement` (secondary) | 1.0 | 1.0 | yes |

Differences are float rounding (≤ 1e-6 relative).

Figures: `fig_spectra.png` shows per-type eigenvalue densities v2 vs v3 for
the four tested types with the null histogram inset and pooled value
marked; `fig_basin.png` shows per-type median basin_index v2 vs v3 with
identity line and ρ. Both as specified.

## Decision rule outcome

1. Symmetry 1.2e-07 / 1.1e-07 ≤ 1e-4 → valid run.
2. 4 types ≥ 3 → continue.
3. `spectral_sep_v2` = 0.73, `spectral_sep_v3` = 1.24, `basin_sep_v3` = 1.04
   are < 2 → **INCONCLUSIVE: operators do not distinguish cell types, so
   agreement is uninformative.** Rules 4–6 are not reached.

## Deviations

None listed. None found.

## Verdict

**ACCEPT** the run. Outcome per the pre-registered rule: **inconclusive**.
T02 stays blocked.

## Planner reading (not findings)

- Had rule 3 not existed, rule 4 would have fired (`spectral_ratio` 19.9,
  `basin_ratio` 6.0 > 3). Within each chemistry, cell types look alike;
  between chemistries, everything differs. The operators are dominated by
  something global to each dataset, not by cell-type wells.
- Every v2 cell (992/992) and 1106/1113 v3 cells are classed `drift`;
  `class_tv` is near zero everywhere, so the class statistic carries little.
- Eigenvalue scale: median |λ| 1906 (v2) vs 1093 (v3), range down to −6269.
  Leading eigenvalues ≈ +2.5 sit far inside tol (191 / 109). The spectrum is
  dominated by very stiff directions. This matches T00b's FA median tangent
  rank of 6: ~24 of 30 FA directions have almost no local spread, so the
  learned density is extremely steep there, and how steep depends on
  chemistry (see T00b review: 135 HVGs absent in v2, 239 near-zero FA noise
  variances in v2 vs 49 in v3).
- So the inconclusive outcome cannot separate "the axiom fails" from "the
  FA positional layer made operators measure off-manifold stiffness". The
  T00b caveats are the leading explanation, not yet tested.

## Quarantined (not evidence)

- Observations: input latent sha256 values (confirmed by recomputation).
- `deviations_requested`: none.
