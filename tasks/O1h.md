# O1h — Diagnostic: does E5b recover the rotation in the high-variance directions?

status: pre-registered (2026-10-17; routine follow-up, Tom's standing instruction)
blocked_by: O1g's run must have exited (PID 8200; memory) — do not run concurrently
seed: 20261014
planner_commit: dcf935d
kind: follow-up (diagnostic; no estimator change, not an O1-series selection candidate)
motivated_by: O1g partial Test R (5/6 cases): err 0.58–0.73 with `L(Q̂) ≈ L(Q*)` — the objective is flat. Planner check on O1b's Test F estimate: error ≤ 0.12 in the top-5 axes but there `Q_load ≈ I`, so Test F cannot tell. O1g did not save per-case `Q̂`.

## Purpose

Decide whether E5b's leftover error is confined to the low-variance
directions. If the estimate is right on the top axes (which carry ~80 % of
the variance), the next estimator restricts the rotation to those
directions; if not, the O1 series moves on to O1c (PROJECT.md rule 14).

## Scope

### In scope
- **One Test R case:** `n = 2759`, `r = 0`, exactly as O1g (`x`, `y`, `Q*`
  from O1g's construction; O1g's 20-seed fields; 20 clean starts; E5 loss;
  LBFGS settings). **No new drift fits:** import O1g's script read-only from
  `~/Documents/scopalign-manus-O1g/reports/O1g/scripts/o1g_run.py`
  (importlib, by absolute path) and use its cached fits for case
  `R_n2759_r0` (O1g cache + O1b fits it reused). Record the sha256 of every
  fit file loaded. Do not write anything under `reports/O1g/`.
- **Save** every start's final `Q` (20 × 30 × 30), its final `L`, and `Q*`.
- **Axes:** the latent axes are PCA components in variance order; check
  that the variances of `x` are non-increasing along the axes (record
  `axes_sorted`; if not, order axes by the variance of `x` and record it).
- **Top-k error** for a rotation `Q` and `k ∈ (1, 2, 3, 5, 8, 10, 15, 20, 30)`:
  `err_k(Q) = ‖(Q − Q*)[:, :k]‖_F / sqrt(2k)` (k = 30 equals O1g's `err`).
  Computed for `Q̂` (lowest `L`), every start, `Q = I`, and the median over
  200 random rotations `ortho_group.rvs(30, random_state=20264000 + i)`.
- **Cell displacement:** `disp(Q) = sqrt(mean_i ‖(Q − Q*) x_i‖² / (2 mean_i ‖x_i‖²))`
  for `Q̂`, `I` and the random median.

### Out of scope
- Other cases, cell counts, starts, loss changes, refits; any change to O1g's
  files or outcome.

## Pre-registration

### Gated numbers
| name | definition | computed by |
|------|-----------|-------------|
| `repro_err_abs` | \|`err_30(Q̂)` − 0.6037296396289642\| (O1g's value for this case) | `python reports/O1h/scripts/o1h_run.py --stage estimate` |
| `repro_L_abs` | \|`L(Q̂)` − 0.14564707632646573\| | same |
| `axes_sorted` | true / false | same |
| `err_top3` | `err_3(Q̂)` | `--stage measure` |
| `err_top5` | `err_5(Q̂)` | same |
| `err_top10` | `err_10(Q̂)` | same |
| `err_top5_random` | median `err_5` over the 200 random rotations | same |
| `n_starts_top5_ok` | starts with `err_5 ≤ 0.10` (of 20) | same |
| `disp_hat` | `disp(Q̂)` | same |
| `disp_random` | median `disp` over the random rotations | same |
| `reading` | `"top_recovered"`, `"partial"` or `"not_recovered"` | `--stage report` |

### Decision rule
0. Validity: `repro_err_abs > 1e-6` or `repro_L_abs > 1e-6` → **REVISE**
   (the case did not reproduce).
1. `reading = "top_recovered"` if `err_top5 ≤ 0.10`;
   `"not_recovered"` if `err_top5 > 0.25`; otherwise `"partial"`.
2. Consequence (planner, not executor): `top_recovered` → draft the
   restricted-subspace estimator (rotation on the top axes only);
   `not_recovered` → O1c next per rule 14; `partial` → planner judgement,
   shown to Tom.

### What would change my mind
If even the top-3 axes are not recovered (`err_top3 > 0.25`), the operator
fields do not fix the rotation anywhere, and restricting the subspace will
not help.

## Procedure

1. `o1h_run.py --stage estimate` — runs E5b on the case (start workers ≤ 3,
   as O1g); `starts.csv` (start, L, err_30), `Q_hats.npz` (`Q` (20,30,30),
   `L` (20,), `Q_star`, `best`), `fits_used.csv` (path, sha256).
2. `o1h_run.py --stage measure` — `topk.csv` (k, err_hat, err_identity,
   err_random_median, n_starts_le_010); `fig_topk.png` (`err_k` vs k: Q̂ bold,
   the 20 starts thin, identity and random median dashed, line at 0.10).
3. `o1h_run.py --stage report` — `report.json`, checksums, lockfile.

`o1h_run.py` must accept `--stage all`. Run detached with stdin from
`/dev/null` (MANUS.md "Long runs"), in worktree
`~/Documents/scopalign-manus-O1h` on branch `task/O1h` from `main`.

## Deliverables

- `report.json`
- `requirements.lock`
- `starts.csv`
- `topk.csv`
- `fits_used.csv`
- `fig_topk.png`
- `scripts/o1h_run.py`

`Q_hats.npz` goes in `artifacts` (with sha256).

## Checker expectations

- `gated` keys: exactly the eleven above (or a subset with
  `"status": "stopped"`).
- `known_failures`: must equal `KNOWN_FAILURES` in PROJECT.md.
