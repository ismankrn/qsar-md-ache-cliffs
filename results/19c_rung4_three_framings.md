# Stage 19c — rung 4 on all three framings: a null on the pair endpoints, a significant loss on potency

Read 2026-10-08, 20 ns window, 158 dev ligands and 86 covered dev pairs. Rung 4 had never
been scored on anything but framing B before this date.

## Results

| framing | endpoint | rung 4 (filtered) | rung 1 2D, same units | gap | gap 95% CI |
|---|---|---|---|---|---|
| **A** | R² per ligand (158) | +0.147 | **+0.413** | **−0.266** | **[−0.508, −0.071]** |
| B | R² per signed pair (86) | +0.317 | +0.321 | −0.004 | [−0.271, +0.393] |
| C | AUC per pair (48/38) | 0.813 | **0.861** | −0.048 | [−0.172, +0.104] |

Classifier families, the symmetric pair construction `[da+db, |da−db|]`, and the metrics are
copied from `rung2_model.py`, so rung 4 is comparable with rung 2 by construction.

## The finding is no longer only a null

**Framing A's interval excludes zero.** This is the first interval in the whole project that
does. It does not say rung 4 is unproven; it says rung 4 is **worse**, with the difference
resolved by the design. For absolute potency, 25 ns of dynamics carries less information
about pChEMBL than a count fingerprint plus RDKit 2D descriptors, and the study can tell.

The two pair endpoints remain genuinely undecided: B's interval is 0.66 wide and C's 0.28
wide, both straddling zero. They are nulls in the honest sense -- not evidence of equality,
just intervals too wide to separate anything.

So the summary across framings is: **significantly worse on potency, indistinguishable on
both pair endpoints.** Not "MD does not help" -- that would overstate B and C -- and not
"MD is promising" -- that would ignore A.

## The 2D baseline is not a subset artefact here

Stage 19b showed the covered subset is biased toward small ligands (median 32 heavy atoms
against 35 outstanding, Mann-Whitney p = 0.009), which is why every baseline is recomputed on
the covered units rather than quoted from the published full-dev number. For framing A that
correction turns out to be small: `05_rung1_benchmark.csv` gives rung 1 framing A a best of
+0.386 over all 217 dev ligands, and the same feature set used here (count ECFP4 + RDKit2D)
reads +0.380 there against +0.413 on the 158-ligand subset. An inflation of +0.033.

That matters for the conclusion. Framing A's result is not an artefact of which ligands
finished first: the matched baseline agrees with the independently computed full-dev value,
and rung 4 is far below both.

## Caveats that belong next to the numbers

**Feature counts are not matched** -- 12 filtered MD descriptors (33 unfiltered) against 2265
2D features. This is a comparison between descriptor sets as they would actually be used, not
a dimensionality-controlled experiment. It is the comparison a practitioner faces, but a
reviewer is entitled to ask for the controlled version, and we have not run it.

**Bootstrap medians are not performance.** Resampling units and then running the frozen folds
over the resample puts duplicates on both sides of the train/test split. Framing A shows the
effect plainly: the rung-4 median is +0.084 where the unresampled fit gives +0.147. Quote the
point estimates for performance and the intervals only for the gap.

**`SVC(probability=True)` is deprecated** in scikit-learn 1.9 and goes away in 1.11. It is in
`CLF` here and in `rung2_model.py`. Replacing it with `CalibratedClassifierCV(SVC(),
ensemble=False)` changes framing C's numbers, so it has to be done for both rungs at once or
the comparability this stage was built for is lost.

## What this closes

The hope recorded in stage 19b -- that framing B might be a null only because it is the
hardest endpoint, with A or C showing signal -- is closed. A is worse and resolved; C is
undecided and pointing the same way. No framing supports rung 4 over 2D.
