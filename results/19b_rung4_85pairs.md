# Stage 19b — rung 4 at 85 of 145 dev pairs: the preliminary advantage does not survive

Read on 2026-10-08, 20 ns window, 85 of 145 dev pairs (156 of 217 phase-1 runs extracted).
The reading threshold of 80 pairs was set in advance, in stage 19, together with what would
be watched at it: not whether the headline R² rose, but whether the lower bound of the gap
rose above zero.

## The numbers

| model set | R² | algorithm |
|---|---|---|
| rung 4, replicate-filtered (12 features) | **+0.317** | k-NN (k=5) |
| rung 4, all descriptors (33 features) | +0.226 | SVM (RBF) |
| rung 1 2D, recomputed on exactly these pairs (2265 features) | **+0.321** | Gradient Boosting |

Gap, rung 4 − rung 1: **−0.004**.

Bootstrap over pairs, 300 replicates, best-of-six re-selected inside every replicate:

| quantity | median | 95% CI |
|---|---|---|
| rung 4 | +0.298 | [+0.058, +0.514] |
| rung 1 2D | +0.267 | [−0.078, +0.539] |
| **gap** | **+0.032** | **[−0.271, +0.393]** |

## What this says

At 59 pairs rung 4 read +0.478 and looked like a result. At 85 pairs it reads +0.317, and a
2D baseline recomputed on the same pairs reads +0.321. The advantage is gone. The honest
statement is not that it is "nearly significant": the interval on the gap is 0.66 wide and
sits almost symmetrically around zero.

**This was close to unobservable.** `bootstrap_gap`, the only function that answers the
question this stage exists to ask, sat below the `if __name__ == "__main__":` block. Python
had not reached its `def` while main ran, and nothing called it from anywhere in the file.
Through the entire preliminary reading the project reported point estimates with no interval
and no one noticed, because the number it would have produced was never missing from the
output -- it had never been in it. A function placed after the block that uses it is not a
late definition, it is an absent one.

## The design cannot rescue this by waiting

Bootstrap width scales about as 1/sqrt(n_pairs). Projecting the measured width forward:

| pairs | projected width | interval around the current median | contains zero |
|---|---|---|---|
| 85 | 0.66 | [−0.300, +0.364] | yes |
| 100 | 0.61 | [−0.274, +0.338] | yes |
| 120 | 0.56 | [−0.248, +0.311] | yes |
| 145 (all) | 0.51 | [−0.223, +0.286] | yes |

For the lower bound to clear zero at the full 145 pairs, the observed gap would have to be
about **+0.254**. It is **+0.032**, an order of magnitude smaller. This agrees with stage 17,
which found the study separates R² values about 0.30 apart; that figure was derived before
any rung-4 run finished and it is holding.

So the remaining 60 pairs will not change the conclusion unless the point estimate moves by
a factor of eight, and nothing in the data so far suggests it will. **As designed, this study
cannot show that MD descriptors beat a matched 2D baseline on framing B.** That is a finding
about the design, not a reason to keep sampling until an interval cooperates.

## One reason the 85 pairs are not simply "145 with more noise"

MD throughput depends on system size, so the ligands that finish first are the smaller ones.
Measured over the 217 phase-1 ligands: median 32 heavy atoms among the 156 extracted against
35 among the 61 outstanding, Mann-Whitney p = 0.009. The covered subset is therefore biased
toward small ligands, and the direction that bias pushes the gap is not known. The reading at
145 pairs will differ from this one systematically, not just by sampling error -- which is an
argument for finishing the runs, but not an argument for expecting a different verdict.

## What is still worth doing

- **Finish the runs.** 61 phase-1 runs remain; all have hosts. The 145-pair reading is the
  one that goes in the paper, biased subsets cannot be reported as if complete, and the
  size-bias measurement above only becomes answerable once the large ligands are in.
- **Score the other framings.** Framing B is the primary endpoint, and it is the hardest:
  a signed difference between two ligands. Framings A and C have not been scored at rung 4.
  A null on B is not a null on all three, and it must not be reported as one.
- **Do not add replicates for this.** Stage 18 showed replicate noise is not the binding
  constraint (descriptor reliability median 0.979, ceiling R² 0.69-1.00). The gap is small
  because the descriptors do not separate these pairs, not because the measurement is noisy.
- **Write the null up.** "Dynamic descriptors did not outperform 2D on signed pair
  differences, with an interval wide enough to state the limit of what the design could
  detect" is a publishable and useful claim, and it is the one the data supports.
