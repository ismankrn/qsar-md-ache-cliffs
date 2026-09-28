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
