#!/usr/bin/env python
"""End-to-end smoke test: SMILES -> OpenFF parameters -> solvated box -> MD.

This is the whole OpenMM pipeline in one file. If this runs, your environment
can do everything the activity-cliff study needs; the production version just
swaps the ligand-in-water box for a ligand-in-protein complex.

Usage:  python smoke_test.py [--platform CPU|OpenCL] [--steps N]
"""
import argparse, time

from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import openmm
from openmm import app, unit

p = argparse.ArgumentParser()
p.add_argument("--smiles", default="CC(=O)Oc1ccccc1C(=O)O", help="ligand SMILES")
p.add_argument("--platform", default="CPU")
p.add_argument("--steps", type=int, default=2000)
p.add_argument("--padding", type=float, default=0.9, help="solvent padding, nm")
p.add_argument("--ligand-ff", default="openff-2.3.0", help="SMIRNOFF force field")
args = p.parse_args()

t0 = time.time()

# 1. Ligand: SMILES -> 3D -> OpenFF parameters (no manual topology editing)
mol = Molecule.from_smiles(args.smiles)
mol.generate_conformers(n_conformers=1)
print(f"[1] ligand: {mol.n_atoms} atoms, formula {mol.to_smiles(explicit_hydrogens=False)}")

smirnoff = SMIRNOFFTemplateGenerator(molecules=mol, forcefield=args.ligand_ff)
ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
ff.registerTemplateGenerator(smirnoff.generator)
print(f"[2] force field: amber14 + tip3p + ligand={args.ligand_ff}")

# 2. Build and solvate
modeller = app.Modeller(mol.to_topology().to_openmm(), mol.conformers[0].to_openmm())
modeller.addSolvent(ff, padding=args.padding * unit.nanometer, model="tip3p")
n_atoms = modeller.topology.getNumAtoms()
print(f"[3] solvated system: {n_atoms:,} atoms")

# 3. System + integrator
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=0.9 * unit.nanometer,
    constraints=app.HBonds,
)
integrator = openmm.LangevinMiddleIntegrator(
    300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
)
platform = openmm.Platform.getPlatformByName(args.platform)
sim = app.Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)
print(f"[4] platform: {sim.context.getPlatform().getName()}")

# 4. Minimise + run
e0 = sim.context.getState(energy=True).getPotentialEnergy()
sim.minimizeEnergy()
e1 = sim.context.getState(energy=True).getPotentialEnergy()
print(f"[5] minimised: {e0.value_in_unit(unit.kilojoule_per_mole):,.0f}"
      f" -> {e1.value_in_unit(unit.kilojoule_per_mole):,.0f} kJ/mol")

sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
t1 = time.time()
sim.step(args.steps)
dt = time.time() - t1

ns = args.steps * 0.002 / 1000
nsday = ns / (dt / 86400)
print(f"[6] ran {args.steps} steps ({ns*1000:.0f} ps) in {dt:.1f}s "
      f"-> {nsday:.2f} ns/day at {n_atoms:,} atoms")
print(f"\nOK - full pipeline works. total {time.time()-t0:.1f}s")
