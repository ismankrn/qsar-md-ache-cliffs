# Stage 13 — building the 217 MD systems

One system per **ligand**, not per trajectory. The 81 cliff-pair ligands get three
replicates each, but a replicate reuses the same system with a different velocity seed —
that is what makes the replicates a test of sampling rather than of structure. 217 systems
carry 379 trajectories.

| | |
|---|---:|
| Ligands (all appearing in a dev-split pair) | 217 |
| — in a cliff pair → 3 replicates | 81 |
| — control only → 1 replicate | 136 |
| Trajectories | 379 |
| Total simulation | 18,950 ns |
| Atoms per system | 96,220 – 96,298 |
| Box | 100.1 Å cubic |
| Build time | 8 min on 24 cores, 0 failures |

Force field follows rung 3 exactly: protein `amber14-all.xml`, ligand **OpenFF 2.3.0 Sage**
via `SMIRNOFFTemplateGenerator`, water TIP3P, 0.15 M NaCl. Rung 4 must match rung 3 or the
rung-to-rung comparison measures force fields instead of descriptors.

Receptor is `docking/4ey7_rec_md.pdb` (stage 12b). Systems are created with
`constraints=None, rigidWater=False` so ParmEd writes every bond explicitly; GROMACS then
applies `constraints = h-bonds` at grompp time and takes rigid water from the `[ settles ]`
block. **FLEXIBLE must not be defined in production** — see
`results/12_gmx_export_validation.md`.

## Guards that are checked per system, not assumed

| guard | result |
|---|---|
| Net charge | max \|q\| = 2.8e-13 e across all 217 |
| Pose fits the capped receptor | **0 heavy–heavy clashes**; closest 2.38 Å |
| Heavy atoms unmoved by H relaxation | fails hard above 1e-6 Å |
| Force field identical everywhere | 1 variant |
| Receptor identical everywhere | 8,173 atoms |

The pose-clearance check initially looked alarming — 42 systems with contacts under 2.0 Å —
but that guard compared receptor heavy atoms against *every* ligand atom. Splitting it
(`results/13b_pose_clearance.csv`) showed heavy–heavy contacts of 2.38–3.52 Å and
heavy–hydrogen of 1.67–2.66 Å: those are hydrogen bonds, not clashes. **No pose fails.**

## The defect this stage found: overlapping ligand hydrogens

`Chem.AddHs(pose, addCoords=True)` places hydrogens at idealised geometry without checking
for overlap. For CHEMBL3632626 that left H50 and H74 **0.880 Å** apart and 132,204 kJ/mol
in the ligand's own nonbonded term — 137,000 kJ/mol above every other system. Same class of
failure as the rung-3 hydrogen problem, inside the ligand this time.

It surfaced only from looking at the *spread* of energies across systems; every system had
built successfully and every grompp passed.

Fixed by relaxing the ligand hydrogens with the heavy atoms frozen at mass 0, so the docked
pose is preserved exactly (checked: heavy-atom shift below 1e-6 Å).

| | before | after |
|---|---:|---:|
| Ligand energy range | −53 … **132,503** | −98 … **1,244** |
| Ligands improved by >1,000 kJ/mol | — | **35 of 217** |
| System energy spread | 147,151 | **50,841** |
| System energy s.d. | 11,214 | **5,890** |

**35 of 217 ligands — 16% — carried more than 1,000 kJ/mol of hydrogen strain.** Fixing
only the one outlier that surfaced would have left 34 others to go into production with the
same defect.

## One outlier remains, and it is solvation

CHEMBL4204913 sits 6.9 s.d. above the median at −806,983 kJ/mol. Its ligand is clean after
relaxation and its receptor contacts are normal, but it has **31 water–water pairs closer
than 1.5 Å** (O···O is normally ~2.8 Å) from `addSolvent`'s lattice placement.

**Checked, not assumed** (`results/14b_npt_audit.csv`, all 379 runs after equilibration):
minimisation and equilibration absorb it completely. CHEMBL4204913 ends up **0.18 s.d.** from
the median, and across every run the potential-energy spread falls from 50,841 kJ/mol to
2,210 (s.d. 5,890 → 396), with temperature inside 0.34 K of 300 and density inside 2.35
kg/m³. Zero physics deviations, zero outliers beyond 5 s.d.

Produced by `scripts/build_md_systems.py`; manifest `results/13_md_ligands.csv`;
per-system provenance in `md/systems/<CID>/build.json`.
