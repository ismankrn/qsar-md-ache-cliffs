"""Stage 13 — build the 217 solvated systems for rung-4 MD and export them for GROMACS.

One system per ligand, not per trajectory: the three replicates of a cliff ligand reuse
the same system with a different velocity seed. That is what makes them a test of sampling
rather than of structure.

Force field follows rung 3 exactly (scripts/rung3_mmgbsa.py): protein amber14SB, ligand
OpenFF Sage 2.3.0 via SMIRNOFFTemplateGenerator, water TIP3P. Rung 4 must match rung 3 or
the rung-to-rung comparison measures force fields instead of descriptors.

Receptor is docking/4ey7_rec_md.pdb -- three capped segments, no fabricated peptide bonds
(see results/12b_receptor_chain_breaks.md). Its gorge is bit-identical to the docking
receptor, so the stage-07 poses still apply; this script verifies that per ligand rather
than assuming it.

The system is created with constraints=None and rigidWater=False so ParmEd writes every
bond explicitly. GROMACS then applies `constraints = h-bonds` at grompp time and takes
water from the [ settles ] block that ParmEd emits under #else. Do NOT define FLEXIBLE in
production: rigid water plus h-bond constraints is what permits the 2 fs timestep.
See results/12_gmx_export_validation.md.

Usage:  python scripts/build_md_systems.py [--index I --total N] [--only CHEMBLxxxx]
        With --index/--total, builds the stride slice I, I+N, I+2N ... for a SLURM array.
"""
import argparse, json, os, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from openmm import app, unit, openmm
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
import parmed
import MDAnalysis as mda

ROOT = Path(__file__).resolve().parent.parent
DOCK, RES = ROOT/"docking", ROOT/"results"
OUT = ROOT/"md"/"systems"
RECEPTOR = DOCK/"4ey7_rec_md.pdb"

CUTOFF  = 1.0*unit.nanometer
PADDING = 1.0*unit.nanometer      # gives ~91k atoms, matching the benchmarked throughput
IONIC   = 0.15*unit.molar
FF_LIG  = "openff-2.3.0"
KJ      = unit.kilojoule_per_mole
CLASH_LIMIT = 1.5                 # A; below this the pose does not fit the receptor

def add_posres(top_path, n_rec, n_lig, rec, off, fc=1000.0):
    """Insert an #ifdef POSRES block into the protein and ligand moleculetypes.

    ParmEd writes no restraints, and GROMACS cannot restrain by selection at run time --
    the indices have to be in the topology, numbered within their own moleculetype.
    Heavy atoms are identified from the mass column rather than from atom names, which
    is robust to whatever naming the ligand happens to carry.
    """
    # Exclusion, not inclusion. The capped receptor has three chains, so ParmEd splits
    # the protein into system1/system2/system3 -- an inclusion list naming "system1"
    # restrained only the first segment and left residues 265-494 and 498-541 free,
    # which would have distorted the protein during restrained NVT. Anything that is not
    # solvent or a monatomic ion gets restrained, whatever it is called.
    SOLVENT = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "T3P",
               "NA", "CL", "K", "MG", "CA", "ZN", "NA+", "CL-"}
    lines = top_path.read_text().splitlines()
    out, current, atoms, in_atoms = [], None, [], False

    def flush(buf, name, entries):
        if name is not None and name.upper() not in SOLVENT and entries:
            buf.append("")
            buf.append("#ifdef POSRES")
            buf.append("[ position_restraints ]")
            buf.append(";  ai   funct        fcx        fcy        fcz")
            for nr in entries:
                buf.append(f"{nr:6d}       1 {fc:10.1f} {fc:10.1f} {fc:10.1f}")
            buf.append("#endif")
            buf.append("")

    for line in lines:
        stripped = line.strip()
        header = stripped.startswith("[") and stripped.endswith("]")
        if header:
            sect = stripped.strip("[] ").lower()
            if sect in ("moleculetype", "system", "molecules"):
                flush(out, current, atoms)
                atoms, current = [], None
            in_atoms = (sect == "atoms")
        elif current is None and out and out[-1].strip().startswith("[ moleculetype"):
            pass
        if in_atoms and not header and stripped and not stripped.startswith(";"):
            f = stripped.split()
            if len(f) >= 8:
                try:
                    if float(f[7]) > 2.0:          # heavy atom
                        atoms.append(int(f[0]))
                except ValueError:
                    pass
        out.append(line)
        # capture the moleculetype name: first data line after [ moleculetype ]
        if current is None and len(out) >= 2:
            for k in range(len(out)-1, max(-1, len(out)-6), -1):
                if out[k].strip().lower().startswith("[ moleculetype"):
                    if stripped and not stripped.startswith(";") and not header:
                        current = stripped.split()[0]
                    break
    flush(out, current, atoms)
    top_path.write_text("\n".join(out) + "\n")

ap = argparse.ArgumentParser()
ap.add_argument("--index", type=int, default=0)
ap.add_argument("--total", type=int, default=1)
ap.add_argument("--only", type=str, default=None)
ap.add_argument("--force", action="store_true")
ap.add_argument("--manifest", default="13_md_ligands.csv",
                help="ligand manifest under results/, e.g. 13b_md_ligands_test.csv")
a = ap.parse_args()

manifest = pd.read_csv(RES/a.manifest)
if a.only:
    manifest = manifest[manifest.molecule_chembl_id == a.only]
else:
    manifest = manifest.iloc[a.index::a.total]
print(f"tugas {a.index}/{a.total}: {len(manifest)} ligan", flush=True)

rec = app.PDBFile(str(RECEPTOR))
rec_heavy = np.array([i for i, at in enumerate(rec.topology.atoms())
                      if at.element is not None and at.element.symbol != "H"])
rec_xyz = np.array(rec.positions.value_in_unit(unit.angstrom))
n_rec = rec.topology.getNumAtoms()
print(f"reseptor: {n_rec} atom, {rec.topology.getNumChains()} rantai", flush=True)

OUT.mkdir(parents=True, exist_ok=True)
built, skipped, failed = 0, 0, []

for _, row in manifest.iterrows():
    cid, seed = row.molecule_chembl_id, int(row.seed)
    d = OUT/cid
    if (d/"system.top").exists() and (d/"system.gro").exists() and not a.force:
        skipped += 1
        continue
    t0 = time.time()
    try:
        d.mkdir(parents=True, exist_ok=True)

        # --- ligand from the docked pose, built FROM the pose so parameters and
        #     coordinates share one atom ordering (the rung-3 bug was exactly this)
        pm = PDBQTMolecule.from_file(str(DOCK/"poses"/f"{cid}_seed{seed}.pdbqt"),
                                     skip_typing=True)
        pose = RDKitMolCreate.from_pdbqt_mol(pm)[0]
        pose = Chem.AddHs(pose, addCoords=True)
        off = Molecule.from_rdkit(pose, allow_undefined_stereo=True)
        lig_top = off.to_topology().to_openmm()
        lig_xyz = np.array(pose.GetConformer(0).GetPositions())
        q_lig = int(round(sum(a_.formal_charge.m for a_ in off.atoms)))

        # --- relax the added hydrogens, heavy atoms frozen
        # Chem.AddHs places hydrogens at idealised geometry without checking for overlap.
        # For CHEMBL3632626 that put H50 and H74 0.880 A apart and left 132,204 kJ/mol in
        # the ligand's own nonbonded term -- 137,000 kJ/mol above every other system. The
        # docked heavy atoms are fine; only the hydrogens are misplaced. Freezing the heavy
        # atoms at mass 0 relaxes the hydrogens and leaves the pose exactly as docked.
        # Same class of failure as the rung-3 hydrogen problem, inside the ligand this time.
        gen_h = SMIRNOFFTemplateGenerator(molecules=off, forcefield=FF_LIG)
        ff_h = app.ForceField()
        ff_h.registerTemplateGenerator(gen_h.generator)
        sys_h = ff_h.createSystem(lig_top, nonbondedMethod=app.NoCutoff)
        is_h = [at.element is not None and at.element.symbol == "H"
                for at in lig_top.atoms()]
        for k in range(sys_h.getNumParticles()):
            if not is_h[k]:
                sys_h.setParticleMass(k, 0.0)
        sim_h = app.Simulation(lig_top, sys_h,
                               openmm.VerletIntegrator(1*unit.femtosecond),
                               openmm.Platform.getPlatformByName("CPU"))
        sim_h.context.setPositions(lig_xyz*unit.angstrom)
        e_h0 = sim_h.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
        sim_h.minimizeEnergy(tolerance=1.0*KJ/unit.nanometer, maxIterations=5000)
        st_h = sim_h.context.getState(getEnergy=True, getPositions=True)
        e_h1 = st_h.getPotentialEnergy().value_in_unit(KJ)
        new_xyz = np.array(st_h.getPositions(asNumpy=True).value_in_unit(unit.angstrom))
        heavy_shift = float(np.abs(new_xyz[~np.array(is_h)] -
                                   lig_xyz[~np.array(is_h)]).max())
        if heavy_shift > 1e-6:
            raise RuntimeError(f"atom berat ligan bergeser {heavy_shift:.2e} A saat relaksasi H")
        lig_xyz = new_xyz

        # --- does this pose still fit the capped receptor?
        dmin = float(np.linalg.norm(
            rec_xyz[rec_heavy][:, None, :] - lig_xyz[None, :, :], axis=2).min())
        if dmin < CLASH_LIMIT:
            raise RuntimeError(f"pose bertabrakan dengan reseptor: {dmin:.2f} A")

        # --- force field
        gen = SMIRNOFFTemplateGenerator(molecules=off, forcefield=FF_LIG)
        ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
        ff.registerTemplateGenerator(gen.generator)

        mod = app.Modeller(rec.topology, rec.positions)
        mod.add(lig_top, lig_xyz*unit.angstrom)
        mod.addSolvent(ff, model="tip3p", padding=PADDING,
                       ionicStrength=IONIC, neutralize=True)
        n_atoms = mod.topology.getNumAtoms()

        # constraints=None / rigidWater=False so every bond is written out; GROMACS
        # applies h-bond constraints and takes rigid water from [ settles ]
        system = ff.createSystem(mod.topology, nonbondedMethod=app.PME,
                                 nonbondedCutoff=CUTOFF,
                                 constraints=None, rigidWater=False)
        nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
        q_tot = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
                    for i in range(system.getNumParticles()))
        if abs(q_tot) > 0.01:
            raise RuntimeError(f"sistem tidak netral: {q_tot:+.3f} e")

        ctx = openmm.Context(system, openmm.VerletIntegrator(1*unit.femtosecond),
                             openmm.Platform.getPlatformByName("CPU"))
        ctx.setPositions(mod.positions)
        st = ctx.getState(getEnergy=True, getPositions=True)
        e0 = st.getPotentialEnergy().value_in_unit(KJ)

        # --- export
        xyz = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        struct = parmed.openmm.load_topology(mod.topology, system, xyz=xyz)
        bv = st.getPeriodicBoxVectors().value_in_unit(unit.angstrom)
        lengths, angles = parmed.geometry.box_vectors_to_lengths_and_angles(
            *[np.array(v) for v in bv])
        pl = lambda x: float(x.value_in_unit(x.unit)) if hasattr(x, "value_in_unit") else float(x)
        struct.box = [*[pl(v) for v in lengths], *[pl(v) for v in angles]]
        struct.save(str(d/"system.top"), overwrite=True)
        struct.save(str(d/"system.gro"), overwrite=True)

        # --- position restraints for equilibration, injected per moleculetype
        add_posres(d/"system.top", n_rec, lig_top.getNumAtoms(), rec, off)

        (d/"build.json").write_text(json.dumps({
            "ligand": cid, "seed": seed, "n_atoms": n_atoms,
            "n_receptor": n_rec, "n_ligand": lig_top.getNumAtoms(),
            "ligand_formal_charge": q_lig, "system_charge": q_tot,
            "box_A": [pl(v) for v in lengths],
            "min_dist_receptor_ligand_A": dmin,
            "ligand_energy_before_H_relax": e_h0,
            "ligand_energy_after_H_relax": e_h1,
            "energy_kJ_per_mol": e0,
            "forcefield": {"protein": "amber14-all.xml", "ligand": FF_LIG,
                           "water": "amber14/tip3p.xml"},
            "n_replicates": int(row.n_replicates),
        }, indent=2))
        built += 1
        print(f"  {cid}: {n_atoms} atom, q_lig {q_lig:+d}, kontak min {dmin:.2f} A, "
              f"E {e0:,.0f} kJ/mol, {time.time()-t0:.0f} s", flush=True)
    except Exception as exc:
        failed.append((cid, f"{type(exc).__name__}: {exc}"))
        print(f"  {cid}: GAGAL {type(exc).__name__}: {exc}", flush=True)

print(f"\ndibangun {built}, dilewati {skipped}, gagal {len(failed)}", flush=True)
for cid, err in failed:
    print(f"  {cid}: {err}", flush=True)
sys.exit(1 if failed else 0)
