#!/usr/bin/env python
"""Rung 3 — single-snapshot MM/GBSA on the docked poses.

    dG = [E_MM + G_polar(GB) + G_nonpolar(SASA)]_complex
         - [same]_protein - [same]_ligand

Single-trajectory throughout: protein, ligand and complex use the SAME
coordinates and nothing is re-minimised after separation, so the internal
bonded terms cancel exactly and only intermolecular interaction and
desolvation remain. Docking used a rigid receptor, so the protein term is
identical for every ligand and is computed once.

Two mistakes in the first version of this script are worth recording:

  1. Charges and Lennard-Jones parameters were taken from a system built from
     SMILES while coordinates came from the docked pose. RDKit and OpenFF
     order atoms differently, so every parameter sat on the wrong atom. The
     symptom was van der Waals energies up to +8,884 kcal/mol. The molecule is
     now built FROM the pose, so ordering is consistent by construction.

  2. The GB polar solvation term was never computed, which left a raw
     electrostatic term of about -479 kcal/mol uncancelled. For cationic
     ligands in an anionic gorge the desolvation penalty is most of the
     answer; without it this is MM, not MM/GBSA.

Not a rigorous free energy method, and not presented as one: reviewer-facing
standards reject absolute MM/PBSA values from single post-docking snapshots.
It is on the ladder because it is what the literature routinely does. Within
a matched pair the systematic error largely cancels, so the DIFFERENCE is far
more defensible than the absolute value.
"""
from __future__ import annotations
import warnings, time
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda
import openmm
from openmm import app, unit
from openff.toolkit import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from rdkit import Chem, RDLogger
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES, DOCK = ROOT/"results", ROOT/"docking"
POSES = DOCK/"poses"
TARGET = "CHEMBL220"
KCAL = unit.kilocalorie_per_mole
COUL, GAMMA = 332.0637, 0.0072
KEY_RES = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]
CPU = openmm.Platform.getPlatformByName("CPU")

rec_pdb = app.PDBFile(str(DOCK/"4ey7_rec_clean.pdb"))
rec_xyz = np.array(rec_pdb.positions.value_in_unit(unit.angstrom))
u = mda.Universe(str(DOCK/"4ey7_rec_clean.pdb"))
res_idx = {r: u.select_atoms(f"resid {r}").indices for r in KEY_RES}

ff_prot = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
sys_prot = ff_prot.createSystem(rec_pdb.topology, nonbondedMethod=app.NoCutoff)
ctx_p = openmm.Context(sys_prot, openmm.VerletIntegrator(1*unit.femtosecond), CPU)
ctx_p.setPositions(rec_pdb.positions)
E_PROT = ctx_p.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
print(f"protein reference (constant): {E_PROT:,.1f} kcal/mol", flush=True)

def nb_params(system):
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    q, s, e = [], [], []
    for i in range(nb.getNumParticles()):
        c, sg, ep = nb.getParticleParameters(i)
        q.append(c.value_in_unit(unit.elementary_charge))
        s.append(sg.value_in_unit(unit.angstrom))
        e.append(ep.value_in_unit(KCAL))
    return np.array(q), np.array(s), np.array(e)

Q_P, S_P, E_P = nb_params(sys_prot)
RAD = {"C":1.7,"N":1.55,"O":1.52,"S":1.8,"H":1.2,"F":1.47,"Cl":1.75,
       "Br":1.85,"I":1.98,"P":1.8}
R_PROT = np.array([RAD.get(a.element.symbol, 1.7) for a in rec_pdb.topology.atoms()])

def sasa(xyz, radii, probe=1.4, n=100):
    g = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*g/n); theta = np.pi*(1+5**0.5)*g
    sph = np.c_[np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)]
    R = radii + probe
    tot = 0.0
    for i in range(len(xyz)):
        pts = xyz[i] + R[i]*sph
        d = np.linalg.norm(pts[:, None] - xyz[None], axis=2)
        d[:, i] = 1e9
        tot += 4*np.pi*R[i]**2 * (d >= R[None]).all(axis=1).mean()
    return tot

def delta_sasa(lx, r_lig, cutoff=12.0):
    """Buried surface on the interface only; distant atoms cancel exactly."""
    near = (np.linalg.norm(rec_xyz[:, None] - lx[None], axis=2).min(axis=1) <= cutoff)
    px, pr = rec_xyz[near], R_PROT[near]
    return (sasa(np.vstack([px, lx]), np.concatenate([pr, r_lig]))
            - sasa(px, pr) - sasa(lx, r_lig)), int(near.sum())

def interaction(lx, ql, sl, el):
    d = np.linalg.norm(rec_xyz[:, None] - lx[None], axis=2)
    d = np.where(d < 0.5, 0.5, d)
    elec = COUL * np.outer(Q_P, ql) / d
    sig = 0.5*(S_P[:, None] + sl[None]); eps = np.sqrt(np.outer(E_P, el))
    sr6 = (sig/d)**6
    vdw = 4*eps*(sr6**2 - sr6)
    per = {f"res_{r}": float(elec[res_idx[r]].sum() + vdw[res_idx[r]].sum())
           for r in KEY_RES}
    return float(elec.sum()), float(vdw.sum()), per

score = pd.read_csv(RES/"07_docking_scores.csv")
ligands = sorted(score.molecule_chembl_id.unique())
out_csv = RES/"09_rung3_mmgbsa.csv"
done = set(pd.read_csv(out_csv).molecule_chembl_id) if out_csv.exists() else set()
print(f"ligands: {len(ligands)}, already done: {len(done)}\n", flush=True)

rows, t0, nfail = [], time.time(), 0
for i, cid in enumerate(ligands, 1):
    if cid in done:
        continue
    try:
        s = score[score.molecule_chembl_id == cid].sort_values("score_1")
        seed = int(s.iloc[len(s)//2].seed)
        pm = PDBQTMolecule.from_file(str(POSES/f"{cid}_seed{seed}.pdbqt"),
                                     skip_typing=True)
        pose = RDKitMolCreate.from_pdbqt_mol(pm)[0]
        pose = Chem.AddHs(pose, addCoords=True)

        # built FROM the pose, so parameters and coordinates share an ordering
        off = Molecule.from_rdkit(pose, allow_undefined_stereo=True)
        gen = SMIRNOFFTemplateGenerator(molecules=off, forcefield="openff-2.3.0")
        ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
        ff.registerTemplateGenerator(gen.generator)

        ltop = off.to_topology().to_openmm()
        lx = np.array(pose.GetConformer(0).GetPositions())
        sys_lig = ff.createSystem(ltop, nonbondedMethod=app.NoCutoff)
        ctx_l = openmm.Context(sys_lig, openmm.VerletIntegrator(1*unit.femtosecond), CPU)
        ctx_l.setPositions(lx*unit.angstrom)
        E_lig = ctx_l.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)

        # complex, same coordinates, carries the GB polar term
        mod = app.Modeller(rec_pdb.topology, rec_pdb.positions)
        mod.add(ltop, lx*unit.angstrom)
        sys_c = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff)
        # Relax ONLY the ligand hydrogens. AddHs places them at idealised
        # geometry, which puts some of them 1.1 A from protein atoms and sends
        # van der Waals energies positive. The docked pose itself is fine:
        # the closest heavy-heavy contact is 3.08 A. Freezing everything else
        # (mass 0) relieves the clash without moving the pose, and keeps the
        # protein term constant.
        n_prot = rec_pdb.topology.getNumAtoms()
        lig_is_h = [a.element.symbol == "H" for a in ltop.atoms()]
        for k in range(sys_c.getNumParticles()):
            free = k >= n_prot and lig_is_h[k - n_prot]
            if not free:
                sys_c.setParticleMass(k, 0.0)
        ctx_c = openmm.Context(sys_c, openmm.VerletIntegrator(1*unit.femtosecond), CPU)
        ctx_c.setPositions(mod.positions)
        openmm.LocalEnergyMinimizer.minimize(ctx_c, 1.0, 200)
        st = ctx_c.getState(getEnergy=True, getPositions=True)
        E_cplx = st.getPotentialEnergy().value_in_unit(KCAL)
        # every downstream term uses these same relaxed coordinates
        allx = np.array(st.getPositions().value_in_unit(unit.angstrom))
        lx = allx[n_prot:]
        ctx_l.setPositions(lx*unit.angstrom)
        E_lig = ctx_l.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)

        dG_mm_gb = E_cplx - E_PROT - E_lig
        ql, sl_, el = nb_params(sys_lig)
        elec, vdw, per = interaction(lx, ql, sl_, el)
        r_lig = np.array([RAD.get(a.element.symbol, 1.7) for a in ltop.atoms()])
        dsasa, n_near = delta_sasa(lx, r_lig)
        dG_np = GAMMA*dsasa

        rows.append({"molecule_chembl_id": cid, "seed_used": seed,
                     "dG_mm_gb": dG_mm_gb, "dG_nonpolar": dG_np,
                     "dG_total": dG_mm_gb + dG_np,
                     "dE_elec": elec, "dE_vdw": vdw, "dE_mm": elec+vdw,
                     "dG_polar_gb": dG_mm_gb - (elec+vdw),
                     "dSASA": dsasa, "n_interface_atoms": n_near,
                     "n_lig_atoms": len(lx), **per})
        if len(rows) >= 10:
            pd.DataFrame(rows).to_csv(out_csv, mode="a", index=False,
                                      header=not out_csv.exists())
            rows = []
        if i % 25 == 0:
            h = (time.time()-t0)/3600
            print(f"  {i}/{len(ligands)}  {h:.2f} h, ETA {h/max(i,1)*(len(ligands)-i):.2f} h",
                  flush=True)
    except Exception as e:
        nfail += 1
        print(f"  FAIL {cid}: {type(e).__name__}: {str(e)[:70]}", flush=True)

if rows:
    pd.DataFrame(rows).to_csv(out_csv, mode="a", index=False,
                              header=not out_csv.exists())
d = pd.read_csv(out_csv)
print(f"\ndone: {len(d)} ligands, {nfail} failed, {(time.time()-t0)/3600:.2f} h")
for c in ["dE_vdw", "dE_elec", "dG_polar_gb", "dG_nonpolar", "dG_total"]:
    print(f"  {c:14} median {d[c].median():>9.1f}   range {d[c].min():>9.1f} to {d[c].max():>9.1f}")
