# Stage 15 — rung-2 geometry on ensemble poses

## Why this stage exists

Stage 11 concluded that receptor rigidity was not what limited rungs 2-3, but it compared
**docking scores** only. Rung-2 descriptors are geometric, and two poses with the same score
can sit differently in the gorge. Recording that as a limitation would have invited exactly
the question a reviewer should ask; this answers it with data.

The ensemble poses had been saved (`docking/ens_poses/`, 2,832 files) — a fact this stage
established late, after the opposite had been asserted. Without them the question would have
been permanently closed.

## Method

Descriptors are computed exactly as `scripts/rung2_descriptors.py` computes them — 4.5 Å
contact cutoff, 6.0 Å cation-π, top-ranked pose, median aggregation — so the numbers are
comparable to `08_rung2_descriptors.csv`. **Each receptor's selections are rebuilt**: CAS,
PAS and gate sidechains, the Trp86 and Trp286 ring centroids, every protein heavy atom, and
the fourteen lining residues. All eight were verified first to carry the right residues with
the right identities and six-atom rings; an empty selection yields wrong descriptors without
raising anything.

Four descriptor sets, one model suite, the frozen CV folds:

| set | what it is |
|---|---|
| `4ey7` | `08_rung2_descriptors.csv` — median over 3 seeds, one rigid receptor |
| `best` | descriptors from whichever of eight receptors scored the ligand lowest |
| `median` | median of each descriptor across all eight receptors |
| `med_fix4ey7` | **control** — all eight poses, every one measured against 4EY7 only |

The published rung-2 figure of R² = 0.167 also used the axis-projection bridging features,
which have no ensemble counterpart. The `4ey7` baseline here excludes them, so all four sets
carry the same features and the comparison is matched rather than inherited.

The control exists because the first result looked like a finding and was not. `median` beat
`4ey7` while `best` did far worse — the pattern noise reduction predicts, not the pattern
"the ensemble finds better poses" predicts. If one receptor genuinely bound a ligand better,
taking *that* receptor's descriptors would win, not averaging all eight. Vina is stochastic,
so one pose is one noisy sample. Holding the measuring receptor fixed at 4EY7 isolates the
eight-fold averaging from the receptor conformations.

## Result — no demonstrable improvement

Framing B, SVM (RBF) fixed a priori, 95% CI from 4,000 bootstrap resamples of the 145
evaluation **pairs** (both antisymmetric rows of a pair kept together, frozen folds untouched):

| set | R² | 95% CI |
|---|---:|---|
| `4ey7` | 0.160 | [+0.026, +0.271] |
| `best` | 0.005 | [−0.186, +0.168] |
| `med_fix4ey7` | 0.207 | [+0.058, +0.339] |
| `median` | **0.274** | [+0.148, +0.384] |

Paired differences, **none significant**:

| comparison | Δ R² | 95% CI | p |
|---|---:|---|---:|
| denoising alone (8 samples vs 3) | +0.047 | [−0.120, +0.212] | 0.56 |
| receptor geometry alone | +0.068 | [−0.037, +0.176] | 0.21 |
| both together | +0.114 | [−0.024, +0.257] | 0.10 |
| best-scoring pose vs 4EY7 | −0.155 | [−0.338, +0.019] | 0.08 |

`median` at R² = 0.274 is significantly above zero, but **not significantly better than
4EY7 at 0.160**. With 145 pairs this study cannot resolve a 0.11 difference in R² from
sampling error.

The split into "0.047 from denoising, 0.068 from receptor geometry" is **not a finding** and
is not claimed as one. Both are point estimates whose intervals contain zero; presenting them
as separable mechanisms would be reading structure into noise.

## What this means for the ladder

| rung | framing B R² |
|---|---:|
| 1 — 2D fingerprints | +0.200 |
| 2 — pose geometry, 4EY7 | +0.160 |
| 2 — pose geometry, ensemble median | +0.274 (not significantly above 4EY7) |
| 3 — single-snapshot MM/GBSA | −0.055 |

The best geometry eight receptor conformations can produce **matches 2D fingerprints**. That
is not a rescue. Stage 11's conclusion stands, now on geometric as well as energetic grounds:
receptor rigidity is not what limited rungs 2-3, which leaves dynamics as the only untested
structural explanation.

One receptor dominates the `best` selection — 6WUY for 216 of 354 ligands, 4EY4 for none.
That looks like a scoring-function preference for a more open gorge rather than a discovery
about binding, and it is consistent with `best` performing worst of the four sets.

## Limitation that remains

**The size confound is untested for these descriptors.** Stage 10 found that after size
adjustment, 0 of 19 rung-2 descriptors separated cliffs from controls, and these models
include `n_heavy` and contact counts that scale with ligand size. Because the ensemble
difference is not significant anyway, repeating the size adjustment would not change the
conclusion — but the caveat belongs in the manuscript.

Produced by `scripts/rung2_ensemble.py` (descriptors), `rung2_ensemble_control.py` (fixed
receptor), `rung2_ensemble_model.py` (model suite), `rung2_ensemble_bootstrap.py` (intervals).
Artefacts: `15_rung2_ensemble_per_pose.csv`, `15_rung2_ensemble_best.csv`,
`15_rung2_ensemble_median.csv`, `15b_rung2_ensemble_median_fixed4ey7.csv`,
`15_rung2_ensemble_benchmark.csv`.
