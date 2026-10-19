# Stage 16d — duration ablation: how much trajectory does rung 4 actually need?

This exists to answer a question before it is asked. The predictable objection to a
descriptor set computed from molecular dynamics is that it costs too much to be useful, and
the honest form of that objection is a real one: a descriptor that needs 50 ns per compound
cannot be applied to a screening library.

The objection is answerable with a curve rather than an argument, but only while the
trajectories exist. The full `prod.xtc` is deleted as soon as descriptors are extracted, and
the reduced archive is 100 ps — too coarse for the rate descriptors, and it carries no water.
So the ablation has to be computed during production or not at all. It was added on
2026-09-28, before any run of the corrected pipeline had completed, so all 30 runs of the
first instance's slice carry it.

## What is computed

Every descriptor is recomputed on the first 5, 10, 20 and 30 ns of analysed trajectory,
alongside the full 45 ns, and stored under a `windows` key in each run's JSON. Windows are in
ns of *analysed* trajectory; the simulation each implies is `EQUIL_NS` (5 ns) longer, and that
is the number to quote as a cost.

| window | simulated | GPU hours | cost per ligand |
|---:|---:|---:|---:|
| 5 ns | 10 ns | 1.4 | $0.22 |
| 10 ns | 15 ns | 2.1 | $0.32 |
| 20 ns | 25 ns | 3.4 | $0.54 |
| 30 ns | 35 ns | 4.8 | $0.75 |
| 45 ns (full) | 50 ns | 6.9 | $1.08 |

At 174 ns/day measured on a rented RTX 3090 at about $0.156/h. Descriptor extraction is
roughly one minute per run and does not enter the comparison.

## How it is computed, and why that matters

The refactor moved the per-frame accumulators from running counters to series, so a window is
a *slice* — `summarise(n)` is called once with the whole trajectory and once per window, and
slicing is the only difference between them. A window therefore cannot be computed by a
different rule than the full result it will be compared against, which is the failure mode an
ablation is most exposed to.

The cost of the ablation itself is nothing measurable: 33 s against 35 s for the same 3,172
frames, within noise, because the trajectory is read once either way.

**Verified bit-identical.** The refactored extraction was run against the deployed version on
the same 3,172-frame snapshot of a live production trajectory: same keys, **zero differing
values**. An ablation whose full-window result had drifted from the production numbers would
be worse than none, so this was checked rather than assumed.

## How to use it

Score the ladder on each window separately and report R² against simulated time. Three
outcomes, all publishable:

- **The curve plateaus early.** Then the practical cost of rung 4 is the plateau, not 50 ns,
  and the objection is answered quantitatively — $0.32 per ligand rather than $1.08.
- **The curve is still climbing at 45 ns.** Then 50 ns was too short, which is a finding about
  the timescale of whatever separates activity cliffs, and it says the next study needs longer
  trajectories rather than more of them.
- **Different descriptors converge at different rates.** The most likely outcome. Occupancies
  and fluctuations should settle early; rate descriptors such as `gate_flips_per_ns` and
  `bridge_runs_per_ns` need long enough to observe several events, and `bridge_longest_ns` is
  bounded by the window by construction and must be read with that in mind.

One caution against over-reading a single run: a ligand that never leaves its pose will look
converged at 5 ns for every descriptor, and says nothing about one that does. The ablation is
only interpretable across the 217 runs, on the same frozen folds as the main analysis.

---

## Decision (2026-10-02): production is 25 ns from here on

Twenty nanoseconds of analysed trajectory, five discarded, so 25 ns simulated. The evidence
is in the replicate pilot rather than here: 12 of the 15 usable descriptors survive at the
20 ns window, and seed noise does not grow as the window shortens — `occ_bridge` has
`sd_rep` 0.006 from 5 to 30 ns against 0.010 at 45. The three descriptors lost are
`frac_lig_buried_mean`, `occ_124` and `gate_chi1_sd`.

What makes this reversible rather than a bet: the tpr is left at its original 25,000,000
steps and the length is imposed with `mdrun -nsteps`, so extending a run needs no
`convert-tpr` at all. Drop the override, resume with `-cpi`, and it continues to the 50 ns
the tpr already specifies. Checkpoints are kept for exactly this.

What is NOT concluded: that 20 ns is enough to *predict*. The pilot measured reproducibility,
not predictive value. With 6 complete cliff pairs of the 145 needed, a correlation estimate
spans [-0.52, +0.93] and settles nothing. The window comparison on the real endpoint becomes
readable at roughly 125 ligands / 48 pairs.

Applied to: vast1 and vast2 from 2026-10-02 (`NSTEPS=12500000`), and BRIN jobs 553833/553834
through `md/NSTEPS_PROD`, which `run_md.sh` reads at runtime — those jobs were already PENDING
and their exported environment could not be changed without discarding a queue position that
took 82 hours to earn.

Not applied to: the 28 phase-1 runs and 12 replicates already finished at 50 ns. They need no
adjustment. Every JSON carries the 5/10/20/30 ns windows, so the primary table is built on the
common 20 ns window for every run regardless of how long it was simulated, and the 50 ns
subset remains available as a robustness check.

---

## Re-run (2026-10-19): recomputed against the re-run replicate filter

The filter this ablation holds fixed was itself recomputed in stage 36, once all six
replicate sets and all 217 ligands had finished, so the curve was regenerated. The
numbers above that describe the mechanics, the cost table and the bit-identical check
are unaffected, as none of them depends on which descriptors the filter selects.

The outcome of the ablation did not change. With the descriptor set held at the 20 ns
selection, framing A is flat across 5, 10 and 20 ns, framing B peaks at 10 ns, and
framing C rises by 0.049 over the same range, every spread being inside the resolution
this study can claim. The 2026-10-02 decision to simulate 25 ns stands, and it was taken
on the replicate pilot rather than on this curve in any case.

Note that the paragraph above recording that decision quotes the filter as it stood then,
i.e. 12 of 15 usable descriptors surviving at 20 ns. That was the state of the evidence
on which the decision was made and is left as written; the current filter passes 11 at
20 ns, which stage 36 records.
