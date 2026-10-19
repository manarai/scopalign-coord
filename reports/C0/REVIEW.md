# C0 review

reviewed_commit: 296fc13 (task/C0)
date: 2026-10-18
verdict: **ACCEPT** (after a planner tolerance correction, approved by Tom) — the Manus cloud sandbox may run T04
recomputation: none needed (identity checks); relative latent error computed from the reported max and the committed K0 latents.

## Results

| check | value | rule |
|---|---|---|
| raw sha256 | 3 / 3 | = 3 ✔ |
| content fingerprints | 8 / 8 | = 8 ✔ |
| scJDO2 | 59cbbf3e… | starts with 59cbbf3 ✔ |
| `latent_check_max` | 1.46e-3 | ≤ 1e-4 ✘ (as registered) |
| `drift_repro_max_abs` | 2.1e-6 | ≤ 1e-4 → cloud and Mac drift fits interchangeable |

Platform Linux x86_64, torch 2.14.1+cpu.

## Planner spec error — latent tolerance

The registered rule gives REVISE. The 1e-4 absolute tolerance was set for
same-machine reruns and is too tight for a float32 PCA (5518 × 2000, full
SVD) computed with a different BLAS/LAPACK. The difference is 1.46e-3 on
latents with max |z| ≈ 31: relative 5e-5, float32 round-off accumulated
through the SVD — not a data or pipeline difference (all 8 content
fingerprints match exactly, and the drift fit on the committed K2 latents
reproduces to 2e-6).

**Correction (Tom, 2026-10-18):** for cross-platform checks the latent
criterion is relative: `max|Δz| / max|z| ≤ 1e-4`. Here 4.7e-5 → passes.
Recorded as a planner error, not an executor deviation.

## Deviations (executor, declared)

Python 3.11 via uv; lock installed with `--no-deps` (lock's internal
dynamo/anndata conflict); `torch==2.14.1+cpu` (CUDA wheel unusable on CPU);
T00 script run with its default command. None can move a gated number; the
drift reproduction (2e-6) shows the environment is equivalent.

## Consequence

T04 may run in this sandbox. Its numbers are internally consistent (all
latents built in the cloud). Drift fits are interchangeable with the Mac's
(2e-6); latents agree to float32 round-off. Rule kept: a T04 number is not
compared with a Mac number without saying so.
