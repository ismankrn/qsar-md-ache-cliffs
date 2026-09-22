#!/usr/bin/env python
"""Benchmark this machine on an AChE-sized solvated system.

Builds a box of roughly the atom count a solvated acetylcholinesterase +
ligand complex needs (~65k atoms) and measures ns/day on each available
OpenMM platform. Water-box proxy: real protein systems run ~10-20% slower
because of bonded terms and a larger PME grid, so treat these as optimistic.
"""
import sys, time
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
import openmm
from openmm import app, unit

TARGET_ATOMS = 65_000
STEPS = {"CPU": 300, "OpenCL": 1000, "Reference": 0}

mol = Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
mol.generate_conformers(n_conformers=1)
smirnoff = SMIRNOFFTemplateGenerator(molecules=mol, forcefield="openff-2.3.0")
ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
ff.registerTemplateGenerator(smirnoff.generator)

print(f"building ~{TARGET_ATOMS:,}-atom box ...", file=sys.stderr)
modeller = app.Modeller(mol.to_topology().to_openmm(), mol.conformers[0].to_openmm())
modeller.addSolvent(ff, padding=4.0 * unit.nanometer, model="tip3p")
n = modeller.topology.getNumAtoms()
print(f"system: {n:,} atoms\n")

system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                         nonbondedCutoff=1.0 * unit.nanometer,
                         constraints=app.HBonds, rigidWater=True)

names = [openmm.Platform.getPlatform(i).getName()
         for i in range(openmm.Platform.getNumPlatforms())]

print(f"{'platform':<12} {'steps':>7} {'wall':>8} {'ns/day':>9}   {'30 ns takes':>14}")
print("-" * 60)
for pname in ("CPU", "OpenCL"):
    if pname not in names:
        continue
    steps = STEPS[pname]
    try:
        plat = openmm.Platform.getPlatformByName(pname)
        integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin,
                                                1 / unit.picosecond,
                                                0.002 * unit.picoseconds)
        sim = app.Simulation(modeller.topology, system, integ, plat)
        sim.context.setPositions(modeller.positions)
        sim.minimizeEnergy(maxIterations=100)
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
        sim.step(20)                       # warm up / JIT
        t = time.time(); sim.step(steps); dt = time.time() - t
        nsday = (steps * 0.002 / 1000) / (dt / 86400)
        days30 = 30 / nsday
        print(f"{pname:<12} {steps:>7} {dt:>7.1f}s {nsday:>9.2f}   "
              f"{days30:>10.1f} days")
    except Exception as e:
        print(f"{pname:<12} FAILED: {type(e).__name__}: {str(e)[:60]}")
