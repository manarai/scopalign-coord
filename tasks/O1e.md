# O1e — Operator-estimated rotation, E8: fused Gromov–Wasserstein matching + operator/drift ICP

status: pre-registered (written 2026-10-15 before any O1-series result; may run in parallel, PROJECT.md rule 14)
blocked_by: none
seed: 20261014
planner_commit: 85b10e7
kind: planned (operator-central alignment, PROJECT.md rules 13–14)

## Purpose

O1's matching estimators pair cells only by their own operator spectra.
Fused Gromov–Wasserstein (FGW) also uses how operators relate **between**
cells within each copy, which does not depend on the frame, so the initial
pairing should be better. E8 is O1's E3 with FGW pairing in place of plain
optimal transport. Same known-answer tests and pass rule as O1.

## Scope

### In scope
- **Everything as `tasks/O1.md` E3** (data, Test R / Test F, drift config,
  five seeds, seed-averaged `M̄`, `ū`, matching groups, loss `H + D`,
  initial rotation by Procrustes on drift vectors of the initial pairs,
  ICP up to 10 iterations, LBFGS, error measures, `err_random`), with one
  change: the **initial matching**.
- **Initial matching (FGW):** per group, feature cost
  `C_ij = ‖λ(M̄_i^X) − λ(M̄_j^Y)‖²` (ascending eigenvalues), divided by its
  median; structure matrices `D^X_ii' = ‖M̄_i^X − M̄_i'^X‖_F`,
  `D^Y_jj' = ‖M̄_j^Y − M̄_j'^Y‖_F` (frame-free within a copy), each divided
  by its median; `P = ot.gromov.fused_gromov_wasserstein(C, D^X, D^Y,
  uniform, uniform, loss_fun="square_loss", alpha=0.5)`; partner of each X
  cell = `argmax_j P_ij`. Groups larger than 1500 cells in either copy are
  first subsampled to 1500 for the FGW plan
  (`numpy.random.default_rng(20261050)`, sorted indices); unsampled X cells
  are paired by plain OT on `C` as in O1.
- **Re-matching in ICP iterations:** as E3 (plain OT on the rotated cost).
- **Views (illustration; not gated):** as O1b.

### Out of scope
- Other `alpha`, subsample sizes, normalisations, thresholds or seeds;
  labels beyond matching groups, gene loadings or positions in the estimate.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_E8_500` | E8, Test R, n = 500 (median of 3) | `python reports/O1e/scripts/o1e_run.py --stage testR` |
| `err_R_E8_1000` | n = 1000 | same |
| `err_R_E8_2000` | n = 2000 | same |
| `err_R_E8_2759` | n = 2759 | same |
| `err_F_E8` | E8, Test F | `--stage testF` |
| `err_random` | as O1 | `--stage testR` |
| `symmetry_max` | as O1 | same |
| `latent_check_max` | as O1 | `--stage testF` |
| `e8_outcome` | `"works"` or `"fails"` | `--stage report` |
| `min_n_e8` | smallest n with `err_R_E8_<n> ≤ 0.10`; 0 if none | same |

### Decision rule
0. Validity: `symmetry_max > 1e-4` or `latent_check_max > 1e-4` → REVISE.
1. E8 **works** if `err_F_E8 ≤ 0.25` and some `err_R_E8_<n> ≤ 0.10`;
   `min_n_e8` = smallest such n. Otherwise **fails**.
Selection across the O1 series follows PROJECT.md rule 14.

### What would change my mind
If better pairing does not rescue the matching family, pairing was not the
bottleneck for E1–E3.

## Procedure

1. `o1e_run.py --stage testR` — 120 fits (kept locally), E8 on 12 cases;
   `testR.csv` (n, replicate, err, n_iter, final_loss).
2. `o1e_run.py --stage testF` — 10 fits, E8; `testF.csv`, `Q_hat_F.npy`;
   views → `views.csv`, `fig_umap.png`.
3. `o1e_run.py --stage report` — `fig_recovery.png`, `report.json`,
   checksums, lockfile.

`o1e_run.py` must accept `--stage all`.

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `views.csv`
- `fig_umap.png`
- `fig_recovery.png`
- `scripts/o1e_run.py`

## Checker expectations

- `gated` keys: exactly the ten above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
