# T04 review

reviewed_commit: b9e7c74 (task/T04; Manus cloud sandbox, after C0)
date: 2026-10-20
verdict: **ACCEPT** — outcome: **not_supported** (as registered)
recomputation: light — verdict re-derived from `scores.csv` and `flags.csv`
(exact); figures inspected.

## Checks

Checker PASS after T04 amendment 1 (planner format error: per-method gated
rows written with a `<m>` placeholder; spelled out, no change of meaning).
Validity: `n_bio_genes` 22, `bio_keep_oracle` 0.939 (≥ 0.3), `symmetry_max`
1.9e-7. Deviations (declared, cannot move a number): per-start
checkpointing; 3 spawn workers × 2 threads (benchmark L agreement 3.9e-12).

## Results

| method | S1 mix_bg | S1 retention | S0 bio_false | removes batch | keeps biology | no false biology |
|---|---|---|---|---|---|---|
| none | 0.34 | 0.84 | 0.24 | ✘ | ✔ | ✘ |
| joint | 0.27 | 1.02 | 0.32 | ✘ | ✔ | ✘ |
| harmony | 0.83 | 0.63 | −0.02 | ✔ | ✘ | ✔ |
| scvi | 0.93 | 0.63 | 0.09 | ✔ | ✘ | ✔ |
| gene_rot | 0.83 | **0.91** | 0.04 | ✔ | ✔ | ✔ |
| **scopalign (E5)** | **0.07** | 0.68 | 0.24 | ✘ | ✘ | ✘ |

Operator flag (S1 responders / max other group / S0 max):
none 5.02 / 4.95 / 4.98; gene_rot 2.19 / 2.19 / 2.42; scopalign 1.50 /
2.32 / 4.06 → **works for no method**.

Rule: scOPAlign does not keep biology (0.68 < 0.8) → **not_supported**.

## Reading

1. **Rigid correction preserves biology; cell-moving methods erode it.**
   The gene-loading rotation is the only method that removes the batch and
   keeps the interferon response (0.91 vs 0.63 for Harmony and scVI) with
   no false biology. This supports the "rigid can't erase biology" half of
   the rationale.
2. **The operator-estimated rotation fails under a realistic batch.** With
   a gene-level (non-rotational) batch effect, E5's objective stays high at
   its best (L 1.30 S0, 0.99 S1, vs ~0.13 on the technical halves) and the
   estimate rotates the halves apart (mix 0.07, below no correction). The
   operator fields are not related by a rotation when the batch is not one,
   so E5 has no correct target.
3. **Operator mismatch does not separate biology from batch.** Under the
   reference rotation, B's monocytes are flagged equally whether they
   responded or not (2.19 vs 2.19), and in S0 with no biology at all (2.42).
   The gene-level batch effect itself changes monocyte operators beyond
   refit noise; the interferon response adds nothing visible on top.

## Consequence (planner; Tom decides)

On this pre-registered test the hypothesis is not supported in its
operator-central form: operators neither estimate the correction nor flag
the biology. What is supported is the value of a rigid correction.
