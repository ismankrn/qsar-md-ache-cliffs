# Stage 16b — two defects in the rung-4 descriptor extraction, and their repair

Found on 2026-09-28, after two of 217 production runs had completed. The question that
exposed it was simply "are you going to plot RMSD and RMSF?" — the plot was never made, and
that is exactly why the defect survived into production.

## What was wrong

`rung4_descriptors.py` read the raw production trajectory and compared coordinates directly
against a frozen reference frame, with no superposition, no minimum-image convention and no
reassembly of molecules across the periodic boundary. The docstring recorded "no realignment"
as a deliberate choice. It is not a defensible one: over 45 ns a solvated protein rotates and
diffuses freely, so a fixed lab-frame reference measures the receptor's tumbling, not the
ligand's behaviour.

Two independent mechanisms, both present:

**1. Free rotation of the complex.** Anything compared against a frozen reference absorbed the
protein's rotation. This alone is enough to destroy the descriptor.

**2. Periodic wrapping.** The capped receptor is three moleculetypes plus the ligand, and mdrun
wraps each independently. In `CHEMBL104553_r1` the ligand sat a whole box vector from its
pocket in 124 of 501 archived frames. Worse, the first receptor molecule was itself cut in
half by a box face, spanning 48.3 Å more than it does in the reference PDB — so placing whole
molecules into a common image was not sufficient either; the molecules had to be made whole
along their bonds first.

The tell-tale was internal contradiction in the very first JSON:

```
rmsf_gorge = 38.4 Å      buried sidechains, which in reality fluctuate ~1 Å
rmsd_mean  = 44.7 Å      for a ligand contacting Trp86 in 99.9% of frames
```

## What was affected, and what was not

The defect splits the descriptor set exactly along one line: whether the quantity is measured
*between* atoms in the same frame, or *against* a frozen reference.

**This section was wrong when first written; the correction follows it.** The original claim
was that only quantities compared against a frozen reference could be affected — seven of
roughly thirty — because everything else measures distances or angles between atoms in the
same frame and is therefore invariant to the protein's rotation.

That reasoning is sound for rotation and silent about the second mechanism. Distances between
atoms in the same frame are rotation-invariant but they are **not** periodic-image-safe, and an
angle among four atoms of one residue is meaningless if a box face happens to cut between them.

It was believed because it held on the only two trajectories then available. Both had already
been deleted, so the test was made on their reduced archives — and `trjconv -pbc mol` had
already made those whole, so the archives could not show an effect that had been present in the
raw data. Measuring it needed a raw trajectory with its water still in it, which existed again
only while the next run was in flight.

### The measurement (28.95 ns, 29,287 waters, 96,239 atoms, identical frames)

| descriptor | old | new | change |
|---|---:|---:|---:|
| `rmsd_mean` | 26.70 | 3.14 | −23.55 |
| `rmsf_gorge` | 37.09 | 0.93 | −36.16 |
| `rmsf_lig` | 35.60 | 1.18 | −34.43 |
| `depth_sd` | 14.31 | 0.60 | −13.70 |
| `lig_dihedral_sd` | 21.72 | 6.57 | **−15.15** |
| `gate_chi1_sd` | 11.18 | 8.97 | **−2.21** |
| `gate_flips_per_ns` | 0.414 | 0.000 | **−0.414** |
| `n_bridge_water` | 1.969 | 2.121 | **+0.152** |
| `occ_447` | 0.506 | 0.585 | **+0.079** |
| `occ_121` | 0.818 | 0.876 | **+0.058** |
| `occ_74`, `occ_341`, `occ_124`, `occ_337` | | | +0.031, +0.027, +0.002, +0.000 |
| `frac_lig_buried_mean` | 0.977 | 0.986 | +0.008 |

Nineteen descriptors differ, not seven. Thirteen are bit-identical, among them `occ_cas`,
`occ_pas`, `occ_bridge`, `bridge_runs_per_ns` and `bridge_longest_ns`.

Three things the original table got wrong:

- **Contact occupancies are affected.** Six of the fourteen lining residues moved, every one of
  them *upward*. A contact whose partner was written an image away is a contact that was
  missed, so the error is a systematic undercount and never an overcount. Trp86 and Trp286
  happen to be unaffected here, which is the only reason the bridging quantities are identical
  — luck of where the box faces fall relative to those residues, not a property that holds by
  construction.
- **`n_bridge_water` is affected, by about 8%.** A water hydrogen-bonded to both ligand and
  protein but written into a neighbouring image is simply not counted. This is the one
  descriptor that no other test could have reached: the reduced archives contain no water.
- **The dihedral spreads are affected, `lig_dihedral_sd` by a factor of three.** A molecule
  straddling a box face has its halves written a box vector apart, so a torsion computed across
  the cut is meaningless. The ligand is its own moleculetype and straddles like any other.

The defensible statement is therefore much narrower than the original one: **only the
quantities that came back bit-identical are demonstrably safe, and only on this trajectory.**
Which descriptors a given run corrupts depends on where its box faces fall relative to its
molecules, so it varies from run to run and cannot be predicted from another run. That is the
reason the two completed runs are being redone rather than patched.

## The repair

1. Molecules are made whole along their bonds once, on the reference frame (`unwrap_mols`).
2. Every later frame is carried forward by rounding each atom's displacement from the previous
   frame to whole box vectors (`nojump_into`). 10 ps between frames against a ~99 Å box leaves
   four orders of magnitude of margin, so motion can never be mistaken for a jump.
3. Molecules are placed into a common periodic image (`make_whole`).
4. `rmsd`/`rmsf` superpose the receptor's Cα onto the reference each frame (`fit_rotation`,
   Kabsch) and carry the ligand along with it.
5. `depth` takes its CAS→PAS axis from the frame it measures, which makes it invariant to
   rotation by construction — no superposition needed at all.
6. `box=` is passed to every distance call, so waters and the three receptor molecules are
   handled by minimum image.

Two guards were added, and both have already earned their place. The extent check against the
reference PDB is what revealed the split molecule — the reduced-archive validation could not
have found it, because `trjconv -pbc mol` had already made those molecules whole. And
`ca_rmsd_mean`, now reported in every JSON, is the one number that cannot look plausible while
the receptor is broken: physical values are 1–3 Å, a broken receptor gives tens.

### Cost

Calling MDAnalysis' bond-based `unwrap` every frame costs 397 ms against 4 ms to read one —
30 minutes per run, 108 GPU-hours over the full design. The nojump scheme reproduces it to
4×10⁻⁷ (largest relative difference, in a circular standard deviation, i.e. summation order)
and runs 28× faster: 6 s against 166 s for 406 frames.

## Validation

Three independent tests, all against the code as deployed.

**A/B on the reduced archives of both finished runs, old code vs new, identical frames.**

| | `CHEMBL104553_r1` | | `CHEMBL1084093_r1` | |
| | old | new | old | new |
|---|---:|---:|---:|---:|
| `ca_rmsd_mean` | — | 1.35 | — | 1.21 |
| `rmsd_mean` | 39.71 | **2.49** | 11.45 | **1.36** |
| `rmsd_max` | 96.83 | **4.00** | 22.82 | **2.39** |
| `rmsf_lig` | 35.48 | **1.88** | 11.00 | **0.86** |
| `rmsf_gorge` | 40.40 | **0.87** | 11.58 | **0.71** |
| `depth_sd` | 5.24 | **0.61** | 7.73 | **0.33** |
| `escaped` | 1 | **0** | 1 | **0** |
| `occ_cas` | 0.44 | 1.00 | 1.00 | 1.00 |
| `frac_lig_buried_mean` | 0.70 | 0.98 | 0.9876 | 0.9876 |

The second run is the cleaner experiment. Every contact descriptor is *identical* between old
and new — that trajectory had no image jumps at all — yet `rmsd_mean` still falls from 11.45 to
1.36 Å and `rmsf_gorge` from 11.58 to 0.71 Å. Rotation is isolated from PBC there, and it alone
accounts for an order of magnitude. The first run adds the PBC mechanism on top.

The unchanged contact descriptors in that run are also the control: they confirm the repair
does not perturb what was already correct.

**Full solvated system with `.tpr` topology**, on a partial trajectory of `CHEMBL118968_r1`
(343 frames, 3.43 ns). This is the only test that exercises water bridges, ligand bonds and the
real 96,239-atom system: `ca_rmsd_mean` 1.04 Å, `rmsd_mean` 1.09 Å, `rmsf_gorge` 0.91 Å,
`escaped` 0, `n_bridge_water` 2.29, `lig_dihedral_sd` 6.70°. All physical.

**Cheap scheme against expensive scheme**, same 406 frames: 37 descriptors, largest relative
difference 4.03×10⁻⁷.

## Consequences for the two completed runs

Their production JSONs carry the bad values, and the water A/B above shows that is nineteen
descriptors rather than the seven first believed. The structural ones can be recomputed from
the reduced archives — validated above — but the archives are 100 ps and hold no water, so
`gate_flips_per_ns`, `bridge_runs_per_ns`, `n_bridge_water` and `lig_dihedral_sd` cannot be
recovered from them (`gate_flips_per_ns` reads 0.00 at 100 ps against 1.02 at 10 ps). Both
ligands are simply run again instead, at about 14 GPU-hours and roughly $2.20.

Neither of them can be cleared by inspection. It was briefly claimed here that
`CHEMBL1084093_r1` was the safer of the two because its contact occupancies came back
unchanged; that was read off six residues whose occupancy was saturated at exactly 1.000, where
an undercount is the only error possible and none had occurred. It says nothing about the other
eight, and nothing at all about its water bridges, which the archive cannot hold. Which
descriptors a run loses depends on where its own box faces fell, and for a deleted trajectory
that is no longer answerable.

## The QC that should have existed first (stage 16c)

`scripts/rung4_qc_plots.py` produces one page per run from the reduced archive — receptor Cα
RMSD, ligand RMSD in the receptor's frame, depth along the gorge axis, per-atom and
per-residue fluctuation, and a CAS/PAS/bridge occupancy timeline — plus a fleet-wide summary
and a table of every run that trips a threshold.

It imports `unwrap_mols`, `nojump_into`, `make_whole` and `fit_rotation` from
`rung4_descriptors.py` rather than reimplementing them. A QC plot carrying its own copy of
that geometry would be able to agree with a broken extraction, which is the one failure a QC
plot must not have. (`rung4_descriptors.py` gained a `__main__` guard so it can be imported;
nothing else about its behaviour as a script changed.)

Two of its thresholds are physical impossibilities rather than preferences — `ca_rmsd_mean`
above 8 Å and `rmsf_gorge` above 5 Å — and exist specifically to catch a regression of this
defect. Run against the two pre-16b JSONs they fire correctly, on all four counts:

```
CHEMBL104553_r1 : EXTRACTION: gorge RMSF 38.4 A is non-physical;
                  EXTRACTION: no ca_rmsd_mean, written by a pre-16b version of the script;
                  ligand left its pose, RMSD 44.7 A; ligand left the gorge at some point
CHEMBL1084093_r1: EXTRACTION: gorge RMSF 11.6 A is non-physical; ...
```

The panels themselves, recomputed from the archives with the repaired geometry, agree with
`analyse()` to two decimals on an independent accumulation path: Cα RMSD 1.35 and 1.21 Å,
ligand RMSD 2.49 and 1.36 Å.

## Requeued, and quarantined in the meantime

Both ligands were put back on the rented instance's queue on 2026-09-28: their JSONs and
reduced archives were deleted there and their names appended to the end of `slice.txt`, so
they are re-run at full 10 ps fidelity after the other 28. Appending to a file that a running
`while read ... done < file` loop is part-way through was verified on that instance's own bash
before doing it, not assumed.

Their defective results are kept here but renamed to `*.json.stale` and `*.xtc.stale`. The
rename is the point: a scoring script globbing `16_rung4_per_run/*.json` would otherwise ingest
two wrong rows among 217 without anything stopping it, which is the same class of silent
failure this record exists to document. `--update` in `pull_vast.sh` will drop the corrected
files in alongside them; the `.stale` copies can be deleted once that has happened.

Any future scoring of rung 4 should require `ca_rmsd_mean` to be present and below 3 Å per run.
Its absence identifies a pre-16b extraction exactly, and no other field does.
