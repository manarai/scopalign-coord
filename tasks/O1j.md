# O1j — Positions choose the axes, operators choose the rotation (E5p)

status: pre-registered (approved by Tom 2026-10-18, with rule 13a)
blocked_by: none
seed: 20261014
planner_commit: b1ce09b
kind: follow-up (O1i failed on its subspace step; Tom 2026-10-18: "subspace from each copy's position covariance, rotation inside it from operators")
motivated_by: reports/O1i/REVIEW.md — operator statistics have flat spectra and their subspaces do not transfer (sub_5 0.84–0.98); each copy's position covariance does (sub_3 0.07, sub_5 0.11). O1h: unrestricted E5b already recovers the top-3 axes to 0.05.

## Purpose

Split the job by what each source can identify (rule 13a):

- **Axes (positions):** each copy's own principal axes, from its own cell
  covariance. No pairing, labels or genes. This says *which* directions
  matter and lines them up only up to sign and up to mixing among axes of
  similar variance.
- **Orientation (operators):** everything the axes leave open — the
  rotation inside the top-`k` block (where axes of similar variance mix)
  and the sign of every remaining axis. These come from the operators and
  drift (E5 loss); O1f showed curvature cannot set signs and flow can.

The correction is then applied to the latent space as before (rule 13).

## Scope

### In scope
- **Data, Test R / Test F, fields:** as `tasks/O1i.md` (n = 2759,
  r ∈ (0, 1, 2), `Q*` seeds as O1g, O1g's 20-seed fits reused read-only,
  sha256 recorded, nothing written under `reports/O1g/`).
- **Axes:** per copy, `C = cov(cells)` (`numpy.cov`, rows = cells),
  `eigh`, descending: `U_X`, `U_Y` (30 × 30). Record eigenvalue gaps.
- **Rotation family:** `Q(R, s) = U_Y · blockdiag(R, diag(s)) · U_Xᵀ`,
  `R ∈ O(k)`, `s ∈ {±1}^(30−k)`. Orthogonal by construction.
- **Step 1 — block rotation `R`:** E5r loss restricted to the top-`k` axes
  exactly as O1i (`P_X = U_X[:, :k] U_X[:, :k]ᵀ`, same for Y), with
  `s = +1` during this step; screen-then-refine starts as O1i (`2^k` sign
  patterns + 20 random `ortho_group.rvs(k, random_state=20265000 + i)`;
  LBFGS from the 8 lowest; workers ≤ 3).
- **Step 2 — complement signs `s`** (flow only, cheap): with `R` fixed, for
  `j = k+1 … 30` in order, set `s_j` to the sign giving the lower full
  flow term `D(Q)` (O1g's `D`, all cells, both directions); sweep again
  until no sign changes (max 5 sweeps). Record the sweeps used.
- **`k`:** gated at **k = 5**; k = 3 reported.
- **Errors (Test R):** `err_k` on x's top-`k` variance axes and `err_30`,
  as O1h; `sub_k` for the position subspace (as O1i); `n_good_starts` as
  O1i.
- **Test F (functional, as O1g):** k = 5; plus the reading variant
  `s = +1` (no sign step) reported, not gated. UMAP views, 4 rows: no
  rotation, O1g's E5b, E5p (k = 5), known answer; mix · sil · frame in titles.

### Out of scope
- New fits; labels, gene loadings or cross-copy pairing anywhere in the
  estimate; positions in anything beyond each copy's own covariance; other
  `k`, sweeps, starts, thresholds or seeds.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_top5_2759` | median over r of `err_5(Q̂)`, k = 5 | `python reports/O1j/scripts/o1j_run.py --stage testR` |
| `err_R_top3_2759` | k = 3 estimator, `err_3` | same |
| `err_R_top10_2759` | k = 5 estimator, `err_10` (block + signs) | same |
| `err_R_full_2759` | k = 5 estimator, `err_30` | same |
| `sub5_max` | max over cases of the position `sub_5` | same |
| `frac_robust` | fraction of the 3 k = 5 cases with `n_good_starts ≥ 3` | same |
| `mix_F` | Test F, k = 5 | `--stage testF` |
| `mix_F_noA` | no rotation | same |
| `sil_F` | k = 5 | same |
| `sil_F_noA` | no rotation | same |
| `frame_ratio_F` | k = 5 | same |
| `symmetry_max` | as O1g (fits loaded) | same |
| `e5p_outcome` | `"works"` / `"fails"` | `--stage report` |

### Decision rule
0. Validity: `symmetry_max > 1e-4` → REVISE.
1. **Test R passes** if `err_R_top5_2759 ≤ 0.10` **and** `frac_robust ≥ 2/3`.
2. **Test F passes** if `mix_F ≥ 0.80`, `frame_ratio_F ≤ 1.5`,
   `sil_F ≥ 0.9 × sil_F_noA`.
3. E5p **works** if both pass → it is the O1-series estimator for O2
   (rule 14: inserted after O1i), with the claim stated as in rule 13a.
   Otherwise **fails** → O1c.
- Reading (not the outcome): `err_R_top10_2759` and `err_R_full_2759`
  show how far the sign step carries beyond the block (axes with close
  variances can mix, so `err_30` is expected to stay well above 0.10).
  Test F `s = +1` vs gated: the sign step's contribution to alignment.

### What would change my mind
If the operators cannot orient a 5-axis block that positions have already
nearly located (`err_R_top5 > 0.10` with `sub5_max ≈ 0.1`), they carry too
little information to correct even the major axes, and the
operator-central plan needs a different operator summary (O1c) or more
cells.

## Procedure

1. `--stage testR` — k ∈ (3, 5), 3 cases: axes, gaps, screen, refine, sign
   sweeps; `screen.csv`; `testR.csv` (k, r, err_k, err_10, err_30, sub_k,
   L_hat, L_star_k, n_good_starts, best_start, sweeps); `Q_hats.npz`.
   One JSON line per finished case in `run.log`.
2. `--stage testF` — `testF.csv` (incl. the `s = +1` variant), `Q_hat_F.npy`,
   `views.csv`, `fig_umap.png`.
3. `--stage report` — `fig_topk.png` (err_k vs k: E5p, O1h's E5b, O1i's
   E5r), `report.json`, checksums, lockfile.

`--stage all`; detached, stdin `/dev/null`; worktree
`~/Documents/scopalign-manus-O1j`, branch `task/O1j` from `main`; commit
`reports/O1j` and push the branch.

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `views.csv`
- `fig_umap.png`
- `fig_topk.png`
- `scripts/o1j_run.py`

`Q_hats.npz`, `screen.csv` and `fits_used.csv` go in `artifacts`.

## Checker expectations

- `gated` keys: exactly the thirteen above (or a subset with `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
