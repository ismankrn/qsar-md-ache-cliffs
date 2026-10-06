# Stage 19 — rung 4 scored, preliminary: the signal is there, the advantage is not yet proven

Read at 52 of 145 dev pairs, on the 20 ns window. Both numbers below move as pairs arrive;
this is recorded now because the threshold set in advance (48 pairs) was crossed, not
because the analysis is finished.

## What was measured

`scripts/rung4_model.py`, built to be comparable with stage 15 by construction: the same
frozen folds (seed 20260902), the same six algorithms, the same antisymmetric pair rows for
framing B, the same preprocessing.

Three choices that the project had already paid to learn:

**One window for every run.** Production was 50 ns until 2026-10-02 and 25 ns after, so
runs carry 45 or 20 ns of analysed trajectory. Every JSON holds the ablation windows, so
the table is built on the common 20 ns window and no ligand is scored on more evidence than
another.

**The replicate filter.** Of ~33 descriptors, the pilot cleared 12 at this window. The
filtered set is the primary result and the unfiltered one is reported beside it.

**A matched baseline.** Rung 1's published +0.200 is over 145 pairs. Comparing a subset
result to it would compare two different problems, so rung 1 is recomputed on exactly the
pairs rung 4 can see.

## Result

| | R² (point) | R² bootstrap median [95% CI] |
|---|---:|---|
| rung 4, filtered (12 features) | **+0.510** | **+0.478 [+0.234, +0.721]** |
| rung 4, all 33 features | +0.355 | — |
| rung 1 2D, same 52 pairs | +0.059 | +0.077 [−0.494, +0.524] |
| **gap** | +0.451 | **+0.406 [−0.027, +0.940]** |

Best algorithm per row; SVM (RBF) for both rung-4 sets, Random Forest for rung 1. The
bootstrap resamples **pairs**, not rows — the two antisymmetric rows of a pair are one
observation — and re-selects the best of six algorithms inside every replicate, because
picking the winner is part of the procedure and fixing it would report a CI for something
nobody ran.

## What can be said, and what cannot

**MD descriptors predict the signed pair difference.** R² = +0.478 with a lower bound of
+0.234: the interval excludes zero on 52 pairs with algorithm selection paid for. That is
the first evidence in this project that rung 4 carries real information.

**They are not yet shown to beat 2D.** The gap's lower bound is −0.027 — it touches zero.
What prevents the claim is not rung 4 but rung 1, whose estimate on this subset spans
−0.494 to +0.524. On 52 pairs, 2D can look like almost anything.

**The filter helps rather than merely tidies.** 12 descriptors beat 33, +0.510 against
+0.355. Dropping the 15 seed-dominated descriptors improves the model, which is what the
replicate pilot predicted and the first independent check of it.

## Why the point estimate will probably fall

The 52 complete pairs are not a random sample of 145 — they are the ligands the slices
happened to reach first. Best-of-six selection inflates a maximum when the spread between
algorithms is wide, and it is: +0.510 for SVM against −0.332 for gradient boosting on the
same features. The bootstrap already discounts this (+0.478 against the +0.510 point
estimate) but cannot correct for a non-random subset.

What is less likely to move is the *direction*. Rank correlation is more robust at small n
than R², and rung 4's Spearman rho runs +0.53 to +0.75 across the sensible algorithms
against +0.45 at best for rung 1.

## When this becomes decidable

Stage 17 calibrated the study to separate R² values about 0.30 apart at 145 pairs. The
observed gap is larger than that, so if it survives it should reach significance well
before the full set — the binding constraint now is rung 1's variance, which tightens with
pairs faster than the gap moves. Re-run at ~80 and ~120 pairs.
