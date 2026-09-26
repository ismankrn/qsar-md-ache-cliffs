"""Stage 12 — prove the OpenFF/OpenMM -> GROMACS export is faithful.

Production runs GROMACS on BRIN (no GPUs there), but the project parameterises
ligands with OpenFF Sage through OpenMM, exactly as rung 3 does. A bad export
does not fail loudly: it produces trajectories that look plausible and are wrong.
That is the same failure mode as the rung 3 atom-ordering bug (+8,884 kcal/mol vdW).

So: build one production-like solvated complex, compute its single-point energy in
OpenMM, export it to GROMACS, and compare term by term.

Several settings differ by default between the two engines and would masquerade as
export errors if left alone. They are pinned here and mirrored in the .mdp:
  - dispersion correction  OpenMM on  / GROMACS off   -> both off
  - vdW cutoff treatment   OpenMM hard cutoff / GROMACS potential-shift -> both raw
  - PME grid and alpha     chosen independently       -> GROMACS told to use OpenMM's
  - constraints            irrelevant to potential energy, but off in both
"""
import json
from pathlib import Path

import numpy as np
from openmm import app, unit
import openmm
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
import parmed
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DOCK, RES = ROOT/"docking", ROOT/"results"
OUT = ROOT/"md"/"export_check"
OUT.mkdir(parents=True, exist_ok=True)

CUTOFF = 1.0*unit.nanometer
PADDING = 1.0*unit.nanometer
KJ = unit.kilojoule_per_mole

# ---------------------------------------------------------------- ligand pose
score = pd.read_csv(RES/"07_docking_scores.csv")
CID = sorted(score.molecule_chembl_id.unique())[0]
s = score[score.molecule_chembl_id == CID].sort_values("score_1")
seed = int(s.iloc[len(s)//2].seed)
print(f"ligan uji : {CID}  (seed {seed})", flush=True)

pm = PDBQTMolecule.from_file(str(DOCK/"poses"/f"{CID}_seed{seed}.pdbqt"), skip_typing=True)
pose = RDKitMolCreate.from_pdbqt_mol(pm)[0]
pose = Chem.AddHs(pose, addCoords=True)
off = Molecule.from_rdkit(pose, allow_undefined_stereo=True)
lig_top = off.to_topology().to_openmm()
lig_xyz = np.array(pose.GetConformer(0).GetPositions())
print(f"atom ligan: {lig_top.getNumAtoms()}", flush=True)

# ---------------------------------------------------------------- force field
gen = SMIRNOFFTemplateGenerator(molecules=off, forcefield="openff-2.3.0")
ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
ff.registerTemplateGenerator(gen.generator)

rec = app.PDBFile(str(DOCK/"4ey7_rec_clean.pdb"))
mod = app.Modeller(rec.topology, rec.positions)
mod.add(lig_top, lig_xyz*unit.angstrom)
print(f"kompleks  : {mod.topology.getNumAtoms()} atom, melarutkan...", flush=True)

mod.addSolvent(ff, model="tip3p", padding=PADDING,
               ionicStrength=0.15*unit.molar, neutralize=True)
n_atoms = mod.topology.getNumAtoms()
print(f"tersolvasi: {n_atoms} atom", flush=True)

system = ff.createSystem(mod.topology,
                         nonbondedMethod=app.PME,
                         nonbondedCutoff=CUTOFF,
                         constraints=None,
                         rigidWater=False)

# pin the settings that otherwise differ silently
for f in system.getForces():
    if isinstance(f, openmm.NonbondedForce):
        f.setUseDispersionCorrection(False)
        f.setUseSwitchingFunction(False)

# one force group per force so the energy can be decomposed
for i, f in enumerate(system.getForces()):
    f.setForceGroup(i)
groups = {i: type(f).__name__ for i, f in enumerate(system.getForces())}

ctx = openmm.Context(system, openmm.VerletIntegrator(1*unit.femtosecond),
                     openmm.Platform.getPlatformByName("CPU"))
ctx.setPositions(mod.positions)

# PME parameters OpenMM actually used, so GROMACS can be told to match
nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
alpha, nx, ny, nz = nb.getPMEParametersInContext(ctx)
# OpenMM returns this either as a Quantity or a bare float depending on version
alpha = alpha.value_in_unit(unit.nanometer**-1) if unit.is_quantity(alpha) else float(alpha)
print(f"PME OpenMM: alpha={alpha:.6f} /nm, grid={nx}x{ny}x{nz}", flush=True)

# Virtual sites would need matching treatment on the GROMACS side; the SMIRNOFF
# generator warned that the force field defines them, so check whether any were
# actually placed on this system.
n_vsite = sum(1 for k in range(system.getNumParticles()) if system.isVirtualSite(k))
print(f"virtual site: {n_vsite}", flush=True)

energies = {}
for i, name in groups.items():
    e = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(KJ)
    energies[name] = energies.get(name, 0.0) + e
total = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)

print("\n=== energi OpenMM (kJ/mol) ===", flush=True)
for k, v in energies.items():
    print(f"  {k:28s} {v:18.4f}")
print(f"  {'TOTAL':28s} {total:18.4f}", flush=True)

# ---------------------------------------------------------------- export
print("\nmengekspor ke format GROMACS via ParmEd...", flush=True)
state = ctx.getState(getPositions=True)
struct = parmed.openmm.load_topology(mod.topology, system,
                                     xyz=state.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
bv = state.getPeriodicBoxVectors().value_in_unit(unit.angstrom)
lengths, angles = parmed.geometry.box_vectors_to_lengths_and_angles(
    *[np.array(v) for v in bv])
# these come back as Quantities in some parmed versions and floats in others
def _plain(x):
    return float(x.value_in_unit(x.unit)) if hasattr(x, "value_in_unit") else float(x)
lengths = [_plain(v) for v in lengths]
angles = [_plain(v) for v in angles]
struct.box = [*lengths, *angles]
print(f"kotak     : {lengths[0]:.2f} x {lengths[1]:.2f} x {lengths[2]:.2f} A, "
      f"sudut {angles[0]:.1f}/{angles[1]:.1f}/{angles[2]:.1f}", flush=True)
struct.save(str(OUT/"system.top"), overwrite=True)
struct.save(str(OUT/"system.gro"), overwrite=True)
print(f"  {OUT/'system.top'}", flush=True)
print(f"  {OUT/'system.gro'}", flush=True)

# .gro stores coordinates to 3 decimals in nm, i.e. 0.01 A. Bond energies are
# quadratic in bond length, so that rounding raises them systematically -- it is
# not symmetric noise. Ship the exact coordinates in a binary .trr as well and
# rerun from it, otherwise the rounding is indistinguishable from a bad export.
import MDAnalysis as mda
xyz = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
u = mda.Universe(str(OUT/"system.gro"))
u.atoms.positions = xyz
u.dimensions = [*lengths, *angles]
with mda.Writer(str(OUT/"system.trr"), n_atoms=u.atoms.n_atoms) as W:
    W.write(u.atoms)
np.save(OUT/"positions_A.npy", xyz)
print(f"  {OUT/'system.trr'}   (koordinat presisi penuh)", flush=True)

# Decisive check: recompute the OpenMM energy from the coordinates GROMACS will
# actually read. .trr is float32, and the GROMACS build is mixed precision, so
# some disagreement is precision rather than parameters. Feeding the same rounded
# coordinates back into OpenMM separates the two: whatever difference survives
# here is a real parameter difference.
u2 = mda.Universe(str(OUT/"system.gro"), str(OUT/"system.trr"))
xyz_rt = np.array(u2.atoms.positions, dtype=float)
ctx.setPositions(xyz_rt*unit.angstrom)
rt = {}
for i, name in groups.items():
    e = ctx.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(KJ)
    rt[name] = rt.get(name, 0.0) + e
rt_total = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
print("\n=== energi OpenMM pada koordinat hasil baca-ulang .trr (kJ/mol) ===", flush=True)
for k, v in rt.items():
    print(f"  {k:28s} {v:18.4f}")
print(f"  {'TOTAL':28s} {rt_total:18.4f}", flush=True)
print(f"  pergeseran koordinat maks: {np.abs(xyz_rt - xyz).max():.3e} A", flush=True)

# What actually decides whether the two engines produce the same dynamics is the
# force on each atom, not the potential energy. A constant energy offset moves
# nothing. Save OpenMM's forces at the exported coordinates so they can be compared
# against GROMACS forces from the same .trr.
ctx.setPositions(xyz*unit.angstrom)
F = ctx.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
    unit.kilojoule_per_mole/unit.nanometer)
np.save(OUT/"forces_openmm_kJ_per_mol_nm.npy", np.asarray(F))
fmag = np.linalg.norm(np.asarray(F), axis=1)
print(f"  gaya OpenMM: rms {fmag.mean():.2f}, maks {fmag.max():.2f} kJ/mol/nm", flush=True)

# How much of any OpenMM-vs-GROMACS force difference is simply PME discretisation?
# Recompute OpenMM's own forces on a much finer reciprocal grid. Whatever OpenMM's
# answer moves by is the size of its own PME error -- a floor below which agreement
# between two independent PME implementations cannot be expected.
fine = int(nx*2)
nb.setPMEParameters(alpha, fine, fine, fine)
ctx.reinitialize(preserveState=True)
F_fine = ctx.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
    unit.kilojoule_per_mole/unit.nanometer)
dF = np.linalg.norm(np.asarray(F_fine) - np.asarray(F), axis=1)
print(f"  PME grid {nx}^3 -> {fine}^3 : pergeseran gaya OpenMM sendiri "
      f"rms {np.sqrt((dF**2).mean()):.4f}, maks {dF.max():.4f} kJ/mol/nm", flush=True)
np.save(OUT/"forces_openmm_finePME.npy", np.asarray(F_fine))
nb.setPMEParameters(alpha, nx, ny, nz)
ctx.reinitialize(preserveState=True)

# ---------------------------------------------------------------- mdp mirror
# ewald-rtol is what GROMACS uses to derive alpha; invert OpenMM's alpha for it.
from math import erfc
def rtol_for(alpha, rc):
    return erfc(alpha*rc)
rtol = rtol_for(alpha, CUTOFF.value_in_unit(unit.nanometer))
(OUT/"rerun.mdp").write_text(f"""; single-point rerun, mirroring the OpenMM system exactly
; ParmEd writes water twice: flexible bonds/angles under #ifdef FLEXIBLE, and a
; rigid [ settles ] block otherwise. GROMACS takes settles by default, which would
; compare rigid water against OpenMM's rigidWater=False -- a different physical
; system, not an export error. For PRODUCTION leave this undefined: rigid water
; plus h-bond constraints is what allows a 2 fs timestep.
define                  = -DFLEXIBLE
integrator              = md
nsteps                  = 0
cutoff-scheme           = Verlet
nstlist                 = 1
pbc                     = xyz
verlet-buffer-tolerance = -1
rlist                   = {CUTOFF.value_in_unit(unit.nanometer)}
coulombtype             = PME
rcoulomb                = {CUTOFF.value_in_unit(unit.nanometer)}
coulomb-modifier        = none
ewald-rtol              = {rtol:.8e}
fourier-nx              = {nx}
fourier-ny              = {ny}
fourier-nz              = {nz}
pme-order               = 5
vdwtype                 = cutoff
vdw-modifier            = none
rvdw                    = {CUTOFF.value_in_unit(unit.nanometer)}
DispCorr                = no
constraints             = none
continuation            = yes
nstfout                 = 1
nstxout                 = 1
nstenergy               = 1
""")
(OUT/"openmm_energies.json").write_text(json.dumps(
    {"ligand": CID, "seed": seed, "n_atoms": n_atoms,
     "alpha_per_nm": alpha, "grid": [nx, ny, nz], "ewald_rtol": rtol,
     "terms": energies, "total": total,
     "terms_trr_roundtrip": rt, "total_trr_roundtrip": rt_total}, indent=2))
print(f"\nmdp dan energi acuan ditulis ke {OUT}", flush=True)
