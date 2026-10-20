# Stage 21 — what a partial subset does to an expensive descriptor set's apparent gain

Expensive descriptors arrive incrementally, and the tempting comparison while they are still
arriving is the subset in hand against the baseline's published figure over the full set. This
project ran that comparison for real: at 59 of 145 pairs rung 4 read +0.478 against rung 1's
published +0.200 and looked like a result. It was not.

## The reconstruction is faithful, not simulated

Completion order comes from each run's JSON mtime, written on the machine that produced it and
preserved through `rsync -a`, so the prefixes are the ligands this project actually held at
each point. Two checks that it is the real history: the 60-pair prefix reproduces **+0.477**
against the +0.478 reported at the time, and the 85-pair prefix **+0.314** against +0.317.

## The table

| subset | ligands | pairs | rung 4 | rung 1 matched | gain vs published | gain vs matched |
|---|---:|---:|---:|---:|---:|---:|
| arrival | 111 | 60 | +0.477 | +0.295 | **+0.277** | **+0.183** |
| random (40 draws) | 111 | — | +0.093 | +0.271 | −0.107 | −0.179 |
| arrival | 153 | 85 | +0.314 | +0.284 | +0.114 | +0.030 |
| random (40 draws) | 153 | — | +0.153 | +0.231 | −0.047 | −0.078 |
| arrival | 196 | 120 | +0.045 | +0.177 | −0.155 | −0.132 |
| random (40 draws) | 196 | — | +0.118 | +0.199 | −0.082 | −0.081 |
| arrival | 217 | 145 | +0.124 | +0.200 | −0.076 | −0.077 |

## Two findings, and the second is the one that matters

**It is not a small-sample effect.** Random subsets of the same size put rung 4 at +0.093
where arrival order puts it at +0.477. What matters is *which* ligands arrive first, not how
many have arrived.

**Recomputing the baseline is necessary and not sufficient.** Matching rung 1 to the same 60
pairs removes 0.094 of the 0.277, and +0.183 survives. The rest goes only when the set is
finished. A study that corrects its baseline onto the available subset and reports the result
is still reporting an inflated number; the correction is not a licence to read a partial set.

## Why the early subset is different — three answers, not one

The convenient claim is ligand size. The first 111 ligands to finish do carry a median of 31
heavy atoms against 35 for the rest (Mann-Whitney p = 0.0005), and the bias decays exactly as
a throughput mechanism predicts: p = 0.0005 at 111 ligands, 0.030 at 153, 0.54 at 196. Random
subsets show no such bias (p = 0.68, 0.42, 0.58), which is the control working.

But two other things travel with arrival order. The early ligands hold **every multi-replicate
run in the study** (mean 1.11 replicates against 1.00 later), and replicate averaging lowers
descriptor noise. They also carry **28.8 ns of analysed trajectory on average against 20.0**,
because the campaign began on a 50 ns design. Smaller systems run faster, replicates were
piloted first, and the long runs came first; all three move together and n = 217 cannot
separate them.

So the supported claim is that the early subset was systematically different in at least three
measured ways and that comparing it against a published full-set baseline manufactured a gain
of +0.277. Not that ligand size caused it.

## Why this generalises

Nothing here is specific to AChE or to molecular dynamics. Any descriptor that is expensive
enough to compute incrementally — MD, FEP, QM, crystallography — delivers its early subset in
an order set by cost, and cost correlates with molecular properties. The reporting practice
that turns this into a false positive, comparing a partial new descriptor set against a
published full-set baseline, is common and cheap to avoid: recompute the baseline on the same
units, report the subset size, and treat anything short of the full set as provisional.

Regenerated 2026-10-20 after the replicate pilot was re-run; the numbers above are the
Oct 12 values and are superseded by `21_subset_inflation.csv`.

Script: `scripts/table_subset_inflation.py`. Outputs `21_subset_inflation.csv` and
`21_subset_inflation.tex`.
