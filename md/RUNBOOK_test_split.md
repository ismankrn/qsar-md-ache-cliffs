# Simulating the 54 sealed test ligands

The development split was simulated; the test split was not. Closing that gap lets the
held-out data validate the whole ladder, rung 4 included, instead of only the four rungs
whose descriptors already cover it.

Nothing needs redocking. All 54 test ligands already have poses under `docking/poses/`
from stage 07, and the receptor is the same capped `4ey7_rec_md.pdb`.

## Scale, from your own measurements

The 97 runs in `md/logs_vast/` took a median of **2.30 h** each on one RTX 3090, which is
**261 ns/day**. Fifty-four runs of 25 ns, one per ligand and no replicates:

| GPUs | wall clock | runs per GPU |
|---|---|---|
| 1 | 124 h (5.2 d) | 54 |
| 2 | 62 h (2.6 d) | 27 |
| **3** | **41 h (1.7 d)** | **18 each, which divides exactly** |
| 4 | 31 h (1.3 d) | 14, 14, 13, 13 |

Total 124 GPU-hours. Reduced archives come back at about 15 MB each, so 0.9 GB for all
54; the full trajectory during a run is around 2 GB and `vast_worker.sh` deletes it once
the descriptors are out.

## Instance specification

What the previous runs used, and what to ask for again:

- **RTX 3090**, one per instance. The timings above are measured on this card.
- **CUDA 12.4 or later.**
- **Disk 250 GB.** A run needs its own 2 GB trajectory plus the system files while it
  works, and the worker deletes as it goes, so this is comfortable rather than tight.
- **4 CPU threads per GPU**, which is what `NTOMP` was set to.
- Three instances, 18 runs each.

Rented on 2026-10-23 and reachable as `vast-1`, `vast-2` and `vast-3`, all RTX 3090
with 250 GB free and 20, 64 and 32 cores. They authenticate with `~/.ssh/vast` and
nothing else, since a BRIN or campus key on a rented box would hand the cluster account
to whoever compromises it.

## Sequence

Each block below says where it runs. Nothing is built on the instances except the
equilibration; they receive systems and return descriptors.

### On the laptop, in the project environment

The base environment has no OpenMM, which is what fails first if you forget.

```bash
conda activate qsar-md
python scripts/build_md_systems.py --manifest 13b_md_ligands_test.csv
python scripts/make_md_runs.py --manifest 13b_md_ligands_test.csv --out 14m_md_runs_test.csv
```

### On the laptop, staging

`MANIFEST` is not optional. Left unset it defaults to `14_md_runs_phase1.csv`, the 217
development runs, and you would stage the wrong campaign without being told.

```bash
MANIFEST=results/14m_md_runs_test.csv bash md/stage_vast.sh vast-1 vast-2 vast-3
```

The script cuts the manifest into one disjoint slice per host by stride and writes each
to `/root/qsar/slice.txt` on its instance. It does its own slicing, so
`make_md_runs.py --slices` is only for looking at the division beforehand and is not a
step.

### On each instance

```bash
cd /root/qsar && SLICE=slice.txt bash vast_worker.sh
```

### On the laptop, afterwards

```bash
bash md/pull_vast.sh
conda run -n qsar-md python scripts/rung4_descriptors.py
```

## Two things to hold to

**The replicate filter does not get recomputed.** It was fixed on the development pilot
before any replicate had been scored, and refitting it now on data that includes the test
split would let the held-out ligands influence which descriptors they are then used to
judge. The 11 descriptors that passed at the 20 ns window stay as they are.

**The models are not retrained or reselected.** Each rung keeps the descriptor arm and
the algorithm family that won on the development split, refits on all 217 development
ligands, and is scored once on the test split. Selection stays on development; the test
split sees each rung exactly once.

## Then

With `results/16_rung4_per_run/` covering all 271 ligands, the test-split evaluation
becomes a single scoring pass over 54 ligands and 38 pairs, and the manuscript gains the
subsection it is missing. Expect intervals about twice as wide as the development ones,
since 38 pairs against 145 widens them by roughly the square root of the ratio.
