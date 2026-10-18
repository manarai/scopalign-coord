# O1i — Operator-estimated rotation restricted to the operator-informative subspace (E5r)

status: pre-registered (approved by Tom 2026-10-18)
blocked_by: none (O1g, O1h ACCEPT)
seed: 20261014
planner_commit: 8afa001
kind: follow-up (O1h reading "partial"; planner judgement approved in principle by Tom 2026-10-18)
motivated_by: reports/O1h/REVIEW.md — E5b recovers the top-3 axes to 0.05 and the error grows smoothly with k (0.12 at 5, 0.60 at 30); only 1 of 20 starts reached it

## Purpose

Estimate the batch rotation **only where the operators carry information**,
and leave every other direction alone. Three changes to E5b:

1. **Subspace from the operators themselves** (rule 13: nothing from
   positions, labels or genes): each copy's top-`k` eigenvectors of a
   frame-covariant operator statistic.
2. **Loss and rotation inside that subspace only** (a `k × k` rotation),
   so the weak, noisy directions no longer pull the estimate.
3. **Starts that cover the sign ambiguity** (O1f: curvature fixes axes up
   to sign; flow decides): all `2^k` sign patterns plus random starts.

The complement is mapped "as unchanged as possible" (no operator
information used there). This is what lets scOPAlign later separate batch
(the rotation), noise (the refit null) and biology (operator mismatch left
after rotation), because weak directions are no longer scrambled.

## Scope

### In scope
- **Data, Test R / Test F constructions, fields:** as `tasks/O1g.md`
  (`n = 2759` only, `r ∈ (0, 1, 2)`, `Q* = ortho_group.rvs(30,
  random_state=20261031 + r)`, 20-seed ensemble fields). **No new drift
  fits:** reuse O1g's fits read-only (import `o1g_run.py` by absolute path
  as O1h did; record sha256 of every fit loaded). Write nothing under
  `reports/O1g/`.
- **Operator statistic** per copy (cells `c` of that copy, ensemble fields,
  float64): `S = Ū / tr Ū + J̃ / tr J̃` with `Ū = mean_c ū(c) ū(c)ᵀ`,
  `J̃ = mean_{c ∈ subset} J̄(c) J̄(c)ᵀ` (the copy's fixed 256-cell `H` subset).
  `V` = top-`k` eigenvectors of `S` (`eigh`, descending), `U⊥` = the other
  `30 − k`. Record each copy's eigenvalue gap `λ_k / λ_{k+1}`.
- **Rotation:** `Q(R) = V_Y R V_Xᵀ + U⊥_Y W U⊥_Xᵀ`, `R ∈ O(k)`,
  `W = polar(U⊥_Yᵀ U⊥_X)` (fixed; closest-to-identity map of the complements,
  no operator information). `Q(R)` is orthogonal.
- **Restricted loss E5r** (all as O1g's E5, with `P_Y = V_Y V_Yᵀ`,
  `P_X = V_X V_Xᵀ`):
  ```
  D_k = [ Σ_i ‖P_Y(ū^Y(Q x_i) − Q ū^X(x_i))‖² / Σ_i ‖P_X ū^X(x_i)‖²
        + Σ_j ‖P_X(ū^X(Qᵀ y_j) − Qᵀ ū^Y(y_j))‖² / Σ_j ‖P_Y ū^Y(y_j)‖² ] / 2
  H_k = same with operators, ‖P_Y(J̄^Y(Q x) − Q J̄^X Qᵀ)P_Y‖_F² etc.
  L_k = D_k + H_k
  ```
- **Parametrisation / optimiser:** `R = R₀ · matrix_exp(A − Aᵀ)`, `A` `k × k`,
  LBFGS as O1g (`max_iter = 500`).
- **`k`:** gated at **`k = 5`**; `k = 3` and `k = 8` run and reported.
- **Starts** (per `k`), screen then refine: candidates `R₀ = diag(s)` for
  all `2^k` sign vectors `s`, plus 20 random
  `ortho_group.rvs(k, random_state=20265000 + i)`, `i = 0…19`. **Screen:**
  evaluate `L_k(Q(R₀))` once for every candidate (no optimisation).
  **Refine:** LBFGS from the 8 lowest (ties by candidate order: signs in
  binary order, then random). Workers ≤ 3 (memory).
- **Errors** (Test R, as O1h): `err_k(Q̂) = ‖(Q̂ − Q*)[:, top-k of x]‖_F / sqrt(2k)`
  with x's variance-ordered axes (for scoring only), plus
  `err_30`. **Subspace check:** `sub_k` = largest principal-angle sine
  between `span(V_Y)` and `Q* span(V_X)` (does the operator subspace
  itself transfer?).
- **Robustness:** `n_good_starts` = refined starts with `L ≤ 1.01 × best L`.
- **Test F:** K0 halves, functional (mix, sil, frame ratio, no rotation vs
  `Q̂` vs `Q_load`), as O1g; UMAP views as O1g (four rows: no rotation,
  E5b from O1g, E5r k = 5, known answer).

### Out of scope
- New fits; positions, labels or gene loadings in `V` or `Q̂`; other `k`,
  statistics, thresholds or seeds.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_top5_2759` | median over r of `err_5(Q̂)`, k = 5, n = 2759 | `python reports/O1i/scripts/o1i_run.py --stage testR` |
| `err_R_top3_2759` | k = 3 estimator, `err_3`, n = 2759 | same |
| `err_R_top8_2759` | k = 8 estimator, `err_8`, n = 2759 | same |
| `sub5_max` | max over the 3 cases of `sub_5` | same |
| `frac_robust` | fraction of the 3 k = 5 cases with `n_good_starts ≥ 3` | same |
| `mix_F` | Test F, k = 5 | `--stage testF` |
| `mix_F_noA` | no rotation | same |
| `sil_F` | k = 5 | same |
| `sil_F_noA` | no rotation | same |
| `frame_ratio_F` | k = 5 | same |
| `symmetry_max` | as O1g (fits loaded) | same |
| `e5r_outcome` | `"works"` / `"fails"` | `--stage report` |

### Decision rule
0. Validity: `symmetry_max > 1e-4` → REVISE.
1. **Test R passes** if `err_R_top5_2759 ≤ 0.10` **and** `frac_robust ≥ 2/3`.
2. **Test F passes** if `mix_F ≥ 0.80`, `frame_ratio_F ≤ 1.5`,
   `sil_F ≥ 0.9 × sil_F_noA`.
3. E5r **works** if both pass; it is then the O1-series estimator for O2
   (rule 14 amended: E5r inserted after O1g), with `k = 5` and the claim
   limited to the operator subspace. Otherwise **fails** → O1c.
- Reading (not the outcome): `sub5_max > 0.3` means the operator subspace
  itself does not transfer between copies (the statistic `S` is too noisy);
  `err_R_top3_2759 ≤ 0.10 < err_R_top5_2759` means `k = 3` is the usable size.

### What would change my mind
If restricting the loss does not beat O1h's unrestricted 0.125 on the top
5 axes, the weak directions were not what limited the estimate; the limit
is the operators' information on those 5 axes, and the usable subspace is
smaller (k ≈ 3).

## Procedure

1. `--stage testR` — for k ∈ (3, 5, 8), 3 cases: `S`, `V`, gaps, screen,
   refine; `screen.csv` (k, r, candidate, L0); `testR.csv` (k, n, r, err_k, err_30, sub_k, L_hat, L_star_k, n_good_starts,
   best_start); `Q_hats.npz` (best `Q̂` per k and case).
2. `--stage testF` — k = 5; `testF.csv`, `Q_hat_F.npy`; `views.csv`,
   `fig_umap.png` (4 rows × batch / cell type, mix·sil·frame in titles).
3. `--stage report` — `fig_topk.png` (err_k vs k for E5r and O1h's E5b),
   `report.json`, checksums, lockfile.

`--stage all`; detached, stdin `/dev/null`; worktree
`~/Documents/scopalign-manus-O1i`, branch `task/O1i` from `main`.

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `views.csv`
- `fig_umap.png`
- `fig_topk.png`
- `scripts/o1i_run.py`

`Q_hats.npz`, `screen.csv` and `fits_used.csv` go in `artifacts`.

## Checker expectations

- `gated` keys: exactly the twelve above (or a subset with `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
