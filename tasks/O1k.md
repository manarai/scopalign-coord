# O1k — Positive control: does the operator estimator recover a known rotation on a synthetic landscape?

status: pre-registered (approved by Tom 2026-11-01)
blocked_by: none (O1g ACCEPT; uses O1g's code and one O1g cached case read-only)
seed: 20261014
planner_commit: 88626c9
kind: diagnostic (positive control; not an estimator candidate, rule 14 unchanged)
motivated_by: Tom 2026-10-29 — "I can't figure out if the problem is the way we are testing it or if it is the operator". Every O1 test so far used PBMC, where neither the true landscape nor its information content is known.

## Purpose

Run the estimator that failed on PBMC (O1g's E5b: flow + curvature field
equivariance, no pairing) on a **synthetic landscape whose true energy is
known**, step by step:

- **Rung A — exact fields:** the true energy's own flow and operators, no
  fitting. If E5b cannot recover the rotation here, the estimator, its
  search or our test is at fault (**test_fault**), and the O1-series
  failures say nothing against operators.
- **Rung B — fitted fields:** scJDO2 drift fits on cells sampled from the
  same landscape, exactly as on PBMC. Failing here but not in A means
  fitting noise destroys the information (**fit_limited**).
- **Anisotropy sweep:** the wells range from round to clearly elongated.
  The rotation is set, outside the directions fixed by the well centres,
  only by the wells' shapes, so this shows how elongated the wells must be.
  Passing only at strong elongation (**anisotropy_limited**) means operators
  can align only landscapes more structured than PBMC; passing in PBMC's
  range (**pbmc_range_ok**) means the PBMC failures come from the real data
  (the two copies' fitted fields disagree beyond a rotation), not from
  operators lacking information.

## Scope

### In scope
- **Landscape (fixed, d = 10):** Gaussian mixture, `K = 4` wells, equal
  weights. `p(x) = ¼ Σ_k N(x; μ_k, Σ_k)`; true energy `E(x) = −log p(x)`
  (β = 1, so `u = −∇E` and `J = −∇²E`, the quantities ESM fits estimate).
  - Means: `μ_k` = `(m_k, 0, …, 0)`, `m_k ∈ ℝ³`, `m_k ~ N(0, 3² I₃)` from
    `numpy.random.default_rng(20264200)`, k = 0…3 in order, one
    `standard_normal(3)` call per k. The means span only the first 3
    coordinates, so the 7 remaining directions ("shape subspace",
    coordinates 4–10 of X's frame) can be oriented only by well shapes.
  - Shapes: `Σ_k = R_k diag(σ²) R_kᵀ`, `R_k = ortho_group.rvs(10,
    random_state=20264210 + k)`, `σ_j² = γ^(5.5 − j)`, j = 1…10 (geometric
    mean 1; adjacent variance ratio `γ`).
  - **Anisotropy `γ` ∈ (1.0, 1.15, 1.33, 1.6)** (condition number γ⁹:
    1, 3.5, 13, 69). PBMC's adjacent operator-eigenvalue ratios were
    1.03–1.33. **γ = 1.0 is the negative control:** round wells with
    centres in 3 coordinates are unchanged by any rotation of the shape
    subspace, so no estimator can recover it; recovering it would mean the
    test leaks the answer.
- **Cases:** `n ∈ (1000, 3000, 10000)` cells per copy × `r ∈ (0, 1, 2)` ×
  4 γ = 36 cases per rung.
  - `Q*_r = ortho_group.rvs(10, random_state=20264100 + r)` (the same for
    every γ and n).
  - Cells: `rng = default_rng(20264300 + r)` per case, draws in this order:
    X labels `rng.integers(0, 4, n)`, X normals `rng.standard_normal((n,
    10))`, then the same two for Y'. `x = μ_label + R_label (σ ⊙ ε)`;
    `Y = Y' · Q*ᵀ` (rows; Y is a fresh sample seen in a rotated frame).
- **Rung A fields:** a torch float64 module `ExactEnergy` returning `E(x)`
  per cell (logsumexp of the four Gaussian log-densities). X's model:
  `E(x)`; Y's model: `E(Q*ᵀ y)`. Ensembles of one model per copy. Fields
  through O1g's `fields`/`ensemble` unchanged.
- **Rung B fields:** drift fits on X and on Y, O1g's `cfg(seed)`
  (hidden 96, depth 3, 1500 epochs, batch 256, ESM, Hutchinson, lr 5e-4),
  input dimension 10, **5 seeds per copy** (20261014, 20261016–19); O1g's
  fit worker, cached, resumable, parallel (`fit_workers()`). Ensembles =
  mean over the 5 seeds (O1g found 20 seeds no better than 5).
- **Estimator:** O1g's E5b — `Problem` (loss `D + H`, `H` on 256-cell
  subsets from `default_rng(20261040 + r)`), `Closure`, LBFGS settings and
  `estimate` — copied into `o1k_run.py` with `D` as a parameter, nothing
  else changed. **Starts:** 20 — `Q₀ = I` and `ortho_group.rvs(10,
  random_state=20264000 + i)`, i = 1…19; `start_overlap` against the three
  `Q*_r` as O1g. Start workers ≤ 3.
- **Copy check:** on O1g's cached Test R case (n = 2759, r = 0; fits read
  from `~/Documents/scopalign-manus-O1g`, read-only, sha256 recorded),
  `L(I)` and `L(Q*_0)` from the copied loss (D = 30) vs O1g's own script
  imported by path → `loss_copy_check` = max absolute difference.
- **Errors per case:** `err = ‖Q̂ − Q*‖_F / √20`;
  `err_shape = ‖(Q̂ − Q*) P_S‖_F / √14`, `P_S` = projector on coordinates
  4–10; `err_pos` likewise on coordinates 1–3 (`/ √6`); `L(Q̂)`, `L(Q*)`,
  `best_start`, number of starts with final `L ≤ 1.01 × L(Q̂)`.
- **Field fit error (rung B):** at each copy's own cells,
  `fe_u = ‖ū − u_exact‖_F / ‖u_exact‖_F` over all cells and
  `fe_J = ‖J̄ − J_exact‖_F / ‖J_exact‖_F` over the 256-cell `H` subset.
- **Views (illustration, not gated):** case γ = 1.33, n = 3000, r = 0:
  three stacked embeddings `[x · Q̂ᵀ ; y]` for no rotation (`Q̂ = I`),
  rung B's `Q̂`, and `Q*`; UMAP as O1g's views; titles with `mix` (30
  nearest other cells from the other copy, ÷ 0.5) and `sil` (well
  labels). Frame ratio is not computed (no refit null defined here).

### Out of scope
- d = 30, other K, mean layouts, γ, n, replicates, seeds, starts, loss
  weights or subset sizes; 20-seed ensembles; any other estimator
  (E5r, E5p, O1c–O1e); non-Gaussian wells; PBMC data beyond the copy
  check; any batch effect other than a rotation.

## Pre-registration

### Gated numbers
Medians are over the 3 replicates; "A" = rung A, "B" = rung B.

| name | definition | computed by |
|------|-----------|-------------|
| `n_cases` | cases completed over both rungs (must be 72) | `python reports/O1k/scripts/o1k_run.py --stage estimate` |
| `start_overlap` | must be 0 | same |
| `loss_copy_check` | copied loss vs O1g's, max abs difference | `--stage check` |
| `symmetry_max` | max `median_asym` over all rung-B fits | `--stage fits` |
| `negctrl_err_shape_A` | A, γ = 1.0, n = 3000: median `err_shape` | `--stage estimate` |
| `errA_max` | A: max over γ ∈ (1.15, 1.33, 1.6) and all n of median `err` | same |
| `errB_g115_n3000` | B, γ = 1.15, n = 3000: median `err` | same |
| `errB_g133_n3000` | B, γ = 1.33, n = 3000 | same |
| `errB_g160_n3000` | B, γ = 1.6, n = 3000 | same |
| `errB_g115_n10000` | B, γ = 1.15, n = 10000 | same |
| `errB_g133_n10000` | B, γ = 1.33, n = 10000 | same |
| `errB_g160_n10000` | B, γ = 1.6, n = 10000 | same |
| `frac_hat_below_star_B` | fraction of the 27 B cases with γ > 1 and `L(Q̂) ≤ L(Q*)` | same |
| `fe_u_3000` | median `fe_u` over B cases with n = 3000 (all γ, r, both copies) | same |
| `fe_J_3000` | median `fe_J`, same cases | same |
| `o1k_outcome` | `"test_fault"`, `"fit_limited"`, `"anisotropy_limited"` or `"pbmc_range_ok"` | `--stage report` |

### Decision rule
0. Validity → REVISE: `n_cases ≠ 72`, `start_overlap ≠ 0`,
   `loss_copy_check > 1e-10`, `symmetry_max > 1e-4`, or
   `negctrl_err_shape_A < 0.30` (the negative control was recovered: the
   test leaks the answer).
1. **test_fault** if `errA_max > 0.10`.
2. Otherwise, with `γ_min` = the smallest γ in (1.15, 1.33, 1.6) whose
   `errB_g<γ>_n3000 ≤ 0.10`:
   - **pbmc_range_ok** if `γ_min ∈ (1.15, 1.33)`;
   - **anisotropy_limited** if `γ_min = 1.6`;
   - **fit_limited** if there is none.
- Reading (not the outcome): the n = 10000 numbers show whether more cells
  rescue a fit-limited or anisotropy-limited result; `frac_hat_below_star_B`
  near 1 with large `err` means fitting noise moved the minimum, near 0
  means the search missed it; `fe_u_3000` and `fe_J_3000` say how wrong the
  fitted fields are.

### What would change my mind
If the estimator recovers the rotation from exact fields but not from
fitted ones even with clearly elongated wells and 10 000 cells, then
scJDO2's fitted operators, as used here, do not carry a rotation, and
operator-estimated alignment should stop; operators remain descriptors
only.

## Procedure

1. `--stage check` — copy check; `start_overlap`; `landscape.npz` (μ, R_k,
   σ² per γ, Q*_r).
2. `--stage fits` — 36 cases × 2 copies × 5 seeds = 360 fits (cached,
   resumable, parallel); `fit_symmetry.csv`.
3. `--stage estimate` — rung A then rung B, every case; `cases.csv` (rung,
   gamma, n, r, err, err_shape, err_pos, L_hat, L_star, best_start,
   n_near_best, fe_u_X, fe_u_Y, fe_J_X, fe_J_Y); `Q_hats.npz`. One JSON line
   per finished case in `run.log`; checkpoint per case.
4. `--stage report` — `fig_ladder.png` (median `err` vs γ: one panel per
   n, rungs A and B, replicate points, lines at 0.10; second row
   `err_shape`), `fig_umap.png` (views), `report.json`, checksums, lockfile.

`--stage all`; detached, stdin `/dev/null`; Mac, worktree
`~/Documents/scopalign-manus-O1k`, branch `task/O1k` from `main`; commit
`reports/O1k` and push the branch after each stage.

## Deliverables

- `report.json`
- `requirements.lock`
- `cases.csv`
- `fit_symmetry.csv`
- `fig_ladder.png`
- `fig_umap.png`
- `scripts/o1k_run.py`

`Q_hats.npz` and `landscape.npz` go in `artifacts`.

## Checker expectations

- `gated` keys: exactly the sixteen above (or a subset with `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
