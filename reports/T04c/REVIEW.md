# T04c review

reviewed_commit: 8bb57c5 (task/T04c; Manus cloud, T04's sandbox)
date: 2026-10-20
verdict: **ACCEPT** — outcome: **weak** (as registered)
recomputation: light — per-seed values in `contrasts.csv` checked against
the gated numbers (contrast = median of per-seed Δ = 0.0883; per-seed nulls
0.043–0.053; seed_agreement 2/5); figure not needed for the verdict.

## Checks

Checker PASS. `inputs_match` 1 (B fits identical to T04's), `symmetry_max`
1.9e-7, `type_ratio_S1` 15.5 (≥ 2: operators readily tell monocytes from
T cells). No deviations.

## Results (inside batch B; ratio = observed ÷ relabelling null)

| comparison | S1 | S0 (control) |
|---|---|---|
| operators: responders vs non-responders | **1.83** (per seed 1.64–2.56; 2 of 5 ≥ 2) | 0.88 (0.83–0.94) |
| positions: responders vs non-responders | **18.3** | 0.93 |
| operators: monocytes vs T cells | 15.5 | — |

Rule: 1.5 < 1.83 < 2.0 and seed agreement 0.4 < 0.8 → **weak**.

## Reading

The interferon response is a real but small operator signal: about twice
refit-level noise inside one batch, with a clean control (S0 ≈ 0.9). The
same response is ~18× noise in positions. The response **moves** cells far
more than it **reshapes** their landscape. Operators are therefore a weak
detector for this kind of biology — ten times less sensitive than the
positions they are fitted on — and could not have flagged it across
batches, where the batch effect alone gives ~2.4 (T04).

## Consequence

Per Tom's order: T04b next (rigid correction vs Harmony vs scVI across
effect sizes and batch strengths). Operators are not part of the method;
whether to keep the within-batch operator contrast as a secondary read-out
in T04b is Tom's call (it costs drift fits per condition).
