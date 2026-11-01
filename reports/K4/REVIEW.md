# K4 — Review

Reviewer: planner (Claude), 2026-11-01. Branch `task/K4` @ 8cb24dc (report data from d90203e), Manus cloud sandbox.

## Checker
`python checker/check_report.py tasks/K4.md reports/K4/report.json --prev reports/K3/report.json` → **PASS** (0 fail, 7 warn). The warnings are cache artifacts (cluster assignments, match state) that live in the sandbox cache and were not committed; they are not deliverables.

## Recomputation
Light recompute, from `fingerprints.npz` on the Mac:
- Re-standardised the raw 30-d operator fingerprints (pooled over both copies) and recomputed the Euclidean distances: max |Δ| vs saved distances = 0.
- Re-ran `linear_sum_assignment` on the operator and expression distances and rescored against the cluster labels: all six accuracies match.
- Regenerated the permutation null with `default_rng(20261400)`, pairs in order P1, P2, P3: identical to the saved nulls (max |Δ| 0); chance and p-value match.
- Not recomputed: Leiden clustering, centroid operators and expression means (cluster assignments are in the sandbox cache).

| gated | value |
|---|---|
| n_pairs_P1 / P2 | 19 / 17 |
| acc_op_P1 / P2 / P3 | 1.000 / 0.824 / 0.824 |
| acc_expr_P1 / P2 / P3 | 1.000 / 1.000 / 1.000 |
| chance_P2 | 0.219 |
| p_op_P2 | 0.001 (minimum possible) |
| symmetry_max | 2.02e-7 |

Figure: `fig_anchors.png` shows the three operator distance matrices with matches marked and the accuracy bars against chance, as asked.
Deviations: none.

## Decision rule outcome
Validity: 19 and 17 pairs ≥ 8; symmetry 2.0e-7 → valid.
Rule 1: p 0.001 ≤ 0.01 and 0.824 ≥ 0.219 + 0.20 → not no_anchor.
Rule 2: 0.824 < 1.000 + 0.10 → not anchors_better.
Rule 3: acc_op_P1 1.00 ≥ 0.80 and acc_op_P2 0.82 ≥ 0.80 → **anchors_equal**.

Wording note (planner): the label "anchors_equal" means "operator anchors pass the 0.80 bar but are not better than expression". Under the gene-level batch they are **worse** than expression (0.82 vs 1.00). The label was chosen before the data and is kept as registered.

## Verdict
**ACCEPT** — outcome **anchors_equal**. Frame-free operator fingerprints match cell clusters across batches without labels: perfectly across two halves of one run, and 14 of 17 under the gene-level batch (chance 0.22). The interferon response (P3) did not reduce accuracy further (0.82). With shared genes, expression means match perfectly in all three pairs. So operator anchors work, but on this data they add nothing over expression.

## Quarantined (not evidence)
- Executor observations: fingerprints came from saved float64 operators; labels are majority coarse types used only for scoring.
- Planner notes (seeds only):
  - The 3 P2 errors are T↔ILC and Monocyte↔DC swaps (a myeloid pair and a lymphoid pair). In the distance matrices, lymphoid clusters (T, B, ILC, Mk, pDC) form one dark block, and operators separate myeloid from lymphoid much better than within either.
  - Coarse labels make "correct" lenient: 6 T-cell and 6 monocyte clusters can match any same-type cluster.
  - Expression's perfect score relies on shared genes and a mild batch. Where operators could add value (no shared genes, or a strong batch) is untested.
