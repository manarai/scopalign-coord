# O1n — Review

Reviewer: planner (Claude), 2026-11-04. Branch `task/O1n` @ 211e0d8, Mac worktree `~/Documents/scopalign-manus-O1n`.

## Checker
`python checker/check_report.py tasks/O1n.md reports/O1n/report.json --prev reports/O1m/report.json` → **PASS** (0 fail, 0 warn).

## Recomputation
Light recompute: every gated value and the rule recomputed exactly from `testR.csv`, `testF.csv` and `rho_kde.csv`. An independent planner run at d = 2 (exploratory, before issue; same construction) reproduced the d = 2 numbers to all digits. Not recomputed: d = 3/5/10 fits, LBFGS.

| d | err_R E5s | err_R skew | mix_F E5s | mix_F skew | mix_F noA | passes |
|---|---|---|---|---|---|---|
| 2 | 0.005 | 0.010 | 1.001 | 0.999 | 0.999 | no (does not beat positions) |
| 3 | 0.002 | 0.013 | 0.994 | 0.992 | 0.992 | no (does not beat positions) |
| 5 | 0.006 | 0.092 | 0.991 | 0.981 | 0.981 | **yes** |
| 10 | 0.016 | 0.467 | 0.979 | 0.605 | 0.344 | **yes** |

Other gated values: symmetry_max 1.8e-7; latent_check_rel 0; rho_kde_d2 0.805; rho_kde_d3 0.151; d_max_op 10.
Test F at d = 10 matches the gene-loading answer (`Q_load` mix 0.979, sil 0.488 vs E5s 0.488). L(Q̂) ≤ L(Q*) in all 12 Test R cases, with err ≤ 0.02.
Deviation (declared): report-local caches and O1l's in-process executor (≤ 3 workers). No effect on values.

## Decision rule outcome
Validity: pass.
Rule 1: d = 5 and d = 10 pass (err ≤ 0.10, mix ≥ 0.80, and beat the skew baseline by ≥ 0.05 in err); d = 2 and 3 do not beat positions.
→ **op_works_lowd**, `d_max_op` = 10.

## Verdict
**ACCEPT** — outcome **op_works_lowd**. In 5 and 10 PCs of real PBMC, operator-estimated rotation (E5s: axis-sign starts plus full refinement) recovers a known rotation (err 0.006 / 0.016). It also aligns the two real halves as well as the gene-defined answer (d = 10: mix 0.34 → 0.98), with no genes, labels or pairing. A positions-only sign rule fails at d = 10 (err 0.47, mix 0.61). In 2–3 PCs, positions alone are enough.

This is the first operator-estimated alignment of real cells that passes. Together with O1m (fails at d = 30) and O1l (synthetic, works at d = 10), the working range is **5–10 dims**, breaking down somewhere between 10 and 30.

Scope of the claim (rule 13a): axes from positions, orientation from operators. Same-run halves only: the frame difference is between two PCAs of one dataset, not a biological or technical batch.

## Quarantined (not evidence)
- rho_kde_d3 = 0.15: the fitted energies track a KDE poorly at d = 3, yet the rotation is recovered. Equivariance seems to need fields that agree between copies, not fields that are accurate. Seed only.
- Untested: a real batch at d = 10 (chemistry or donors); the break point between 10 and 30; correcting the remaining PCs beyond d.
