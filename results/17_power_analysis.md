# Stage 17 — what effect this study can detect, and what caps it

Run before committing money to rung 4, because the answer decides what a null result is
allowed to mean. Both questions are answered from data already in hand.

## Noise ceiling: not the limiting factor

pChEMBL aggregates assays across labs, so part of Δp is measurement noise no descriptor can
explain. 336 of 5,609 compounds (6.0%) carry a repeat measurement; their between-assay
spread is median 0.064, mean 0.184, p90 0.560 log units.

| σ taken as | ceiling R² |
|---|---:|
| median, 0.064 | 0.995 |
| mean, 0.184 | 0.958 |
| p90, 0.560 (pessimistic) | **0.616** |

Even pessimistically the ceiling is 0.62, far above rung 1's 0.200. **About 80% of the
variance in Δp is unexplained and is not measurement noise** — there is real signal that no
descriptor family has found yet.

Caveat stated rather than buried: only 6% of compounds have a repeat, the spread is skewed,
and that subset may not represent the rest.

## Power: the study sees large effects, not moderate ones

Synthetic descriptor matrices pushed through the **same** SVM and the **same** frozen 5-fold
CV every rung was scored with, compared against a synthetic baseline tuned to rung 1's 0.200
by paired bootstrap over the 145 dev pairs (150 simulations, 800 resamples).

| candidate R² | mean gain over baseline | power | CI width |
|---:|---:|---:|---:|
| 0.125 | −0.100 | 1% | 0.33 |
| **0.335** | +0.114 | **36%** | 0.30 |
| **0.479** | +0.256 | **90%** | 0.28 |

**Rung 4 must reach roughly R² = 0.45–0.48 to be separated convincingly from rung 1.** A
moderate improvement to 0.335 would be detected about a third of the time. Interpolating,
power at R² = 0.40 is roughly 50–60%.

The CI width sits near 0.30 across the whole range, which is the more portable statement:
**this study can separate two R² values about 0.30 apart, and not less.**

### Two earlier versions of this analysis were wrong

Recorded because the conclusions were quoted before they were checked.

1. The first simulation generated *predictions* rather than descriptors and left them
   uncalibrated. A signal with correlation √r² scores R² = 2√r² − 1, so a true 0.20 appeared
   as −0.104 — a systematic penalty mistaken for sampling error. It made the study look far
   weaker than it is.
2. The corrected version still bypassed model fitting and treated rung 1's 0.200 as a
   constant. Both sharpen the comparison: real out-of-fold predictions carry the variance of
   fitting on 4/5 of the data, and a baseline that is itself an estimate is harder to beat.
   That version claimed **100% power at R² = 0.40**; the proper analysis puts it near 55%.

What is reported above runs descriptors through the real pipeline and re-simulates the
baseline every time. It is still an upper bound: the synthetic descriptors are Gaussian and
uncorrelated with ligand size, while real MD descriptors will not be — and stage 10 showed
the size confound already destroyed four rung-2/3 "findings".

## Consequence for the manuscript

A null at rung 4 is reportable, but only as: **"an improvement to R² ≥ 0.45 can be ruled
out"**, not "dynamics do not explain activity cliffs". The likeliest outcome given rungs 1-3
(0.200, 0.160, −0.055) and ensemble rung 2 (0.274, not significant) is a value near 0.30-0.35
that is higher than rung 1 but not significantly so — the same ambiguity stage 15 already
produced. The full design was chosen with that risk understood.

Produced by `scripts/power_and_ceiling.py` and `scripts/power_proper.py`.
