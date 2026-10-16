# O1d — Operator-estimated rotation, E6: one shared energy with a learnable rotation

status: pre-registered (written 2026-10-15 before any O1-series result; may run in parallel, PROJECT.md rule 14)
blocked_by: none
seed: 20261014
planner_commit: 85b10e7
kind: planned (operator-central alignment, PROJECT.md rules 13–14)

## Purpose

Make the rotation part of the operator fit itself. One energy network is
trained by exact score matching on the union of copy Y and copy X passed
through a learnable rotation. The score-matching objective is lowest when
both copies describe one landscape, so the rotation that minimises it is
the operator-derived alignment. Same known-answer tests and pass rule as O1.

## Scope

### In scope
- **Data, Test R / Test F constructions, error measures, `err_random`:**
  exactly `tasks/O1.md`. Latents float32 → float64 inside training.
- **Model:** `scjdo2.energy_drift._EnergyMLP(dim=30, hidden=96, depth=3)`
  (the network `EnergyDriftAdapter` uses), float64, CPU; weights
  initialised with `torch.manual_seed(20261014)` before construction.
- **Rotation:** `Q = Q₀ · matrix_exp(S − Sᵀ)`, `S` a free 30×30 parameter
  initialised at 0.
- **Training data each step:** union of all Y cells (as is) and all X cells
  mapped by `x ↦ Q x` (recomputed with the current `Q`); mini-batches of 256
  drawn from the union by `torch.randperm` under
  `torch.Generator().manual_seed(20261014)`; 1500 epochs (one epoch = one
  pass over the union).
- **Loss (exact score matching, Hutchinson trace, as scJDO2's `esm`):** for
  a batch `b`, `s = −∇E(b)`;
  `L = mean_b [ ½‖s‖² + tr(∇s) ]`, with `tr(∇s)` estimated by Hutchinson
  with 4 Rademacher probes per sample (drawn from the same generator).
  Gradients flow to both the network weights and `S`.
- **Optimiser:** `torch.optim.Adam` with two parameter groups — network
  `lr = 5e-4, weight_decay = 1e-4`; `S`: `lr = 1e-3, weight_decay = 0`.
- **Starts:** 5 — `Q₀ = I` and `ortho_group.rvs(30, random_state=20261021 + i)`
  for `i = 1…4`. Each start trains from the same initial weights. `Q̂` = the
  start with the lowest **final full-data ESM loss** (evaluated on the whole
  union with 4 probes under `manual_seed(20261099)`).
- **Views (illustration; not gated):** as O1b.

### Out of scope
- Pairing cells, labels, gene loadings or positions in the estimate.
- Other architectures, epochs, learning rates, probe counts, start counts,
  seeds or thresholds; any modification of the scJDO2 checkout (the loop is
  written in `reports/O1d/scripts/`, importing `_EnergyMLP` only).

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `err_R_E6_500` | E6, Test R, n = 500 (median of 3) | `python reports/O1d/scripts/o1d_run.py --stage testR` |
| `err_R_E6_1000` | n = 1000 | same |
| `err_R_E6_2000` | n = 2000 | same |
| `err_R_E6_2759` | n = 2759 | same |
| `err_F_E6` | E6, Test F | `--stage testF` |
| `n_starts_near_best_F` | Test F starts with final ESM loss ≤ best + 1% of \|best\| | same |
| `err_random` | as O1 | `--stage testR` |
| `latent_check_max` | as O1 | `--stage testF` |
| `e6_outcome` | `"works"` or `"fails"` | `--stage report` |
| `min_n_e6` | smallest n with `err_R_E6_<n> ≤ 0.10`; 0 if none | same |

### Decision rule
0. Validity: `latent_check_max > 1e-4` → REVISE.
1. E6 **works** if `err_F_E6 ≤ 0.25` and some `err_R_E6_<n> ≤ 0.10`;
   `min_n_e6` = smallest such n. Otherwise **fails**.
Selection across the O1 series follows PROJECT.md rule 14.

### What would change my mind
If a shared score-matching fit does not pull the copies into one frame,
the energy landscape alone does not fix the rotation, and operator
alignment needs explicit field equations (O1b) or summaries (O1c).

## Procedure

1. `o1d_run.py --stage testR` — 12 cases × 5 starts trainings; `testR.csv`
   (n, replicate, start, final_loss, err).
2. `o1d_run.py --stage testF` — 5 trainings; `testF.csv`, `Q_hat_F.npy`;
   views → `views.csv`, `fig_umap.png`.
3. `o1d_run.py --stage report` — `fig_recovery.png`, `report.json`,
   checksums, lockfile.

`o1d_run.py` must accept `--stage all`.

## Deliverables

- `report.json`
- `requirements.lock`
- `testR.csv`
- `testF.csv`
- `Q_hat_F.npy`
- `views.csv`
- `fig_umap.png`
- `fig_recovery.png`
- `scripts/o1d_run.py`

## Checker expectations

- `gated` keys: exactly the ten above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
