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

| | Affected | Why |
|---|---|---|
| `rmsd_mean`, `rmsd_max` | yes | compared to a frozen reference |
| `rmsf_lig`, `rmsf_gorge` | yes | fluctuation about a lab-frame mean |
| `depth_mean`, `depth_sd` | yes | axis and origin frozen at frame 0 |
| `escaped` | yes | threshold on the same displacement |
| `occ_*`, `occ_cas/pas/bridge` | no | inter-atomic distances, rotation-invariant |
| `bridge_runs_per_ns`, `bridge_longest_ns` | no | derived from those contacts |
| `gate_chi1_sd`, `gate_flips_per_ns`, `lig_dihedral_sd` | no | internal dihedrals |
| `n_bridge_water`, `frac_lig_buried_mean` | no | inter-atomic distances |

Seven of roughly thirty descriptors. The contact and dihedral descriptors — which carry the
CAS–PAS bridging hypothesis, the actual scientific content of rung 4 — were never affected.

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

Their production JSONs carry the seven bad values. The structural ones can be recomputed from
the reduced archives — validated above — but the archives are 100 ps and hold no water, so
`gate_flips_per_ns`, `bridge_runs_per_ns`, `n_bridge_water` and `lig_dihedral_sd` cannot be
recovered from them (`gate_flips_per_ns` reads 0.00 at 100 ps against 1.02 at 10 ps). Mixing
fidelities across 217 runs is not worth the saving: both ligands should simply be run again,
which costs about 14 GPU-hours, roughly $2.20.

`CHEMBL1084093_r1` is the less urgent of the two — its contact descriptors are provably
unaffected. `CHEMBL104553_r1` had real image jumps, so its `occ_*` cannot be trusted to better
than a few hundredths.

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
