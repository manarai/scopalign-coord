# T01x review (exploratory)

reviewed_commit: a5d797a (task/T01x)
date: 2026-10-15
verdict: REJECT — one gated number does not match its definition

## Checker

`check_report.py tasks/T01x.md reports/T01x/report.json --prev reports/T01c/report.json`
→ PASS (0 fail, 0 warn). Branch changes only `reports/T01x/`. Input ops
sha256 (ab5238dc… / 7b1fb01d…) match the committed T01c files.

## Recomputation

Planner's own script, reading `reports/T01c/ops_*.npz` as committed:

| gated | report | recomputed | match |
|---|---|---|---|
| `scale_v2` | 30.2874 | 30.2874 | yes |
| `scale_v3` | 4.80332 | 4.80332 | yes |
| `scale_ratio` | 6.3055 | 6.3055 | yes |
| `spectral_w2_pooled_s` | 11.778 | 11.778 | yes |
| `spectral_w2_null_s` | 0.75438 | 0.75438 | yes |
| `spectral_ratio_s` | 15.613 | 15.613 | yes |
| `basin_w1_pooled_s` | 1.0558 | 1.0558 | yes |
| `basin_w1_null_s` | 0.19483 | 0.19483 | yes |
| `basin_ratio_s` | 5.4188 | 5.4188 | yes |
| `class_tv_ratio_s` | 6.508188 | 6.508188 | yes |
| `max_sep_change` | 4.4e-16 | **4.8e-05** | **no** |

`max_sep_change` is defined as the largest |scaled separation − **T01c's
reported value**|. The executor's script (lines 190–211) compares the scaled
separations with an unscaled recomputation inside the same run, which is a
different quantity. Unlisted deviation from a gated definition → REJECT.

Planner error found in the same check: computed as defined, the value is
4.8e-05 because T01c's eigenvalues are stored as float32. The pre-registered
validity threshold (< 1e-6) is therefore too tight, so even a correct run
would have hit reading rule 1 (REVISE). Fixing that requires changing the
task, hence a new task id rather than a revision.

Figure `fig_spectra_scaled.png`: per-type scaled eigenvalue densities v2 vs
v3 with null inset; as specified.

## Reading-rule outcome

Not applied (REJECT). For the record only, not a finding: the ten matching
numbers give `spectral_ratio_s` 15.6 and `basin_ratio_s` 5.4, both > 3,
which under rule 3 would read "shape differs".

## Verdict

**REJECT.** Nothing from this run enters FINDINGS. Branch `task/T01x` kept
until Tom confirms deletion. Proposed replacement: T01y, identical to T01x
except the validity check (relative tolerance 1e-4 against T01c's reported
values, stated explicitly).

## Quarantined (not evidence)

- Observations: input ops sha256 values (confirmed).
