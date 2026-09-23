#!/usr/bin/env python
"""Rung 3 — single-snapshot MM/GBSA on the docked poses.

    dG = [E_MM + G_polar(GB) + G_nonpolar(SASA)]_complex
         - [same]_protein - [same]_ligand

Single-trajectory throughout: the protein, ligand and complex all use the
SAME coordinates, and nothing is re-minimised after separation. That is what
makes the internal bonded terms cancel exactly, leaving only intermolecular
interaction and desolvation. Re-minimising the separated parts breaks the
cancellation and turns the result into the difference of two large unrelated
numbers.

Docking used a rigid receptor, so the protein term is identical for every
ligand and is computed once.

This is NOT a rigorous free energy method, and is not presented as one.
Reviewer-facing standards reject absolute MM/PBSA values from single
post-docking snapshots. It sits on the ladder because it is what the
published literature routinely does, so beating it means beating current
practice. Within a matched pair the systematic error largely cancels, which
makes the DIFFERENCE far more defensible than the absolute value.

Components and a per-residue decomposition along the gorge are reported, not
just the total: the decomposition is where a bridging ligand would show
itself, by drawing energy from Trp86 and Trp286 at once.
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
COUL = 332.0637          # kcal*A/(mol*e^2)
GAMMA = 0.0072           # kcal/mol/A^2, standard nonpolar coefficient
KEY_RES = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]

# ---- receptor: fixed for every ligand --------------------------------
rec_pdb = app.PDBFile(str(DOCK/"4ey7_rec_clean.pdb"))
rec_xyz = np.array(rec_pdb.positions.value_in_unit(unit.angstrom))
u = mda.Universe(str(DOCK/"4ey7_rec_clean.pdb"))
res_idx = {r: u.select_atoms(f"resid {r}").indices for r in KEY_RES}

ff_prot = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
sys_prot = ff_prot.createSystem(rec_pdb.topology, nonbondedMethod=app.NoCutoff,
                                constraints=None)
ctx_prot = openmm.Context(sys_prot, openmm.VerletIntegrator(1*unit.femtosecond),
                          openmm.Platform.getPlatformByName("CPU"))
ctx_prot.setPositions(rec_pdb.positions)
E_PROT = ctx_prot.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
print(f"protein reference energy (constant): {E_PROT:,.1f} kcal/mol", flush=True)

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


def interaction(lig_xyz, ql, sl, el):
    """Direct protein-ligand nonbonded interaction, decomposed."""
    d = np.linalg.norm(rec_xyz[:, None] - lig_xyz[None], axis=2)   # (nprot, nlig)
    d = np.where(d < 0.5, 0.5, d)
    elec = COUL * np.outer(Q_P, ql) / d
    sig = 0.5 * (S_P[:, None] + sl[None])
    eps = np.sqrt(np.outer(E_P, el))
    sr6 = (sig / d) ** 6
    vdw = 4 * eps * (sr6**2 - sr6)
    per_res = {f"res_{r}": float(elec[res_idx[r]].sum() + vdw[res_idx[r]].sum())
               for r in KEY_RES}
    return float(elec.sum()), float(vdw.sum()), per_res


def sasa(xyz, radii, probe=1.4, n=120):
    """Shrake-Rupley solvent accessible surface area."""
    g = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*g/n); theta = np.pi*(1+5**0.5)*g
    sph = np.c_[np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)]
    R = radii + probe
    total = 0.0
    for i in range(len(xyz)):
        pts = xyz[i] + R[i]*sph
        d = np.linalg.norm(pts[:, None] - xyz[None], axis=2)
        d[:, i] = 1e9
        total += 4*np.pi*R[i]**2 * (d >= R[None]).all(axis=1).mean()
    return total

RAD = {"C":1.7,"N":1.55,"O":1.52,"S":1.8,"H":1.2,"F":1.47,"Cl":1.75,
       "Br":1.85,"I":1.98,"P":1.8}
prot_elem = np.array([a.element.symbol for a in rec_pdb.topology.atoms()])
R_PROT = np.array([RAD.get(e, 1.7) for e in prot_elem])

score = pd.read_csv(RES/"07_docking_scores.csv")
lig_tab = pd.concat([
    pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id","smiles_md"]],
    pd.read_csv(RES/f"04_validation_ligands_{TARGET}.csv")[["molecule_chembl_id","smiles_md"]],
]).drop_duplicates("molecule_chembl_id").set_index("molecule_chembl_id")
ligands = sorted(score.molecule_chembl_id.unique())
print(f"ligands: {len(ligands)}\n", flush=True)

out_csv = RES/"09_rung3_mmgbsa.csv"
done = set()
if out_csv.exists():
    done = set(pd.read_csv(out_csv).molecule_chembl_id)
    print(f"resuming, {len(done)} already done\n", flush=True)

rows, t0 = [], time.time()
for i, cid in enumerate(ligands, 1):
    if cid in done:
        continue
    try:
        # median-seed top pose, matching rung 2
        s = score[score.molecule_chembl_id == cid].sort_values("score_1")
        seed = int(s.iloc[len(s)//2].seed)
        pm = PDBQTMolecule.from_file(str(POSES/f"{cid}_seed{seed}.pdbqt"),
                                     skip_typing=True)
        lmol = RDKitMolCreate.from_pdbqt_mol(pm)[0]

        off = Molecule.from_smiles(lig_tab.loc[cid, "smiles_md"],
                                   allow_undefined_stereo=True)
        off.generate_conformers(n_conformers=1)
        gen = SMIRNOFFTemplateGenerator(molecules=off, forcefield="openff-2.3.0")
        ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
        ff.registerTemplateGenerator(gen.generator)

        ltop = off.to_topology().to_openmm()
        sys_lig = ff.createSystem(ltop, nonbondedMethod=app.NoCutoff)
        ql, sl_, el = nb_params(sys_lig)

        lx = np.array(lmol.GetConformer(0).GetPositions())
        if len(lx) != sys_lig.getNumParticles():      # pose has no hydrogens
            lmol_h = Chem.AddHs(lmol, addCoords=True)
            lx = np.array(lmol_h.GetConformer(0).GetPositions())
        if len(lx) != sys_lig.getNumParticles():
            raise RuntimeError(f"atom count {len(lx)} vs {sys_lig.getNumParticles()}")

        # ligand alone, in its bound conformation - no re-minimisation
        ctx_l = openmm.Context(sys_lig, openmm.VerletIntegrator(1*unit.femtosecond),
                               openmm.Platform.getPlatformByName("CPU"))
        ctx_l.setPositions(lx * unit.angstrom)
        E_lig = ctx_l.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)

        elec, vdw, per_res = interaction(lx, ql, sl_, el)

        lelem = [a.element.symbol for a in ltop.atoms()]
        r_lig = np.array([RAD.get(e, 1.7) for e in lelem])
        sasa_l = sasa(lx, r_lig)
        allx = np.vstack([rec_xyz, lx]); allr = np.concatenate([R_PROT, r_lig])
        sasa_c = sasa(allx, allr)
        sasa_p = sasa(rec_xyz, R_PROT)
        dG_np = GAMMA * (sasa_c - sasa_p - sasa_l)

        rows.append({"molecule_chembl_id": cid, "seed_used": seed,
                     "E_lig_internal": E_lig,
                     "dE_elec": elec, "dE_vdw": vdw, "dE_mm": elec + vdw,
                     "dSASA": sasa_c - sasa_p - sasa_l, "dG_nonpolar": dG_np,
                     "dG_total": elec + vdw + dG_np, **per_res})
        if len(rows) % 10 == 0:
            pd.DataFrame(rows).to_csv(out_csv, mode="a", index=False,
                                      header=not out_csv.exists())
            rows = []
        if i % 25 == 0:
            el_h = (time.time()-t0)/3600
            print(f"  {i}/{len(ligands)}  {el_h:.2f} h elapsed, "
                  f"ETA {el_h/i*(len(ligands)-i):.2f} h", flush=True)
    except Exception as e:
        print(f"  FAIL {cid}: {type(e).__name__}: {str(e)[:70]}", flush=True)

if rows:
    pd.DataFrame(rows).to_csv(out_csv, mode="a", index=False,
                              header=not out_csv.exists())
d = pd.read_csv(out_csv)
print(f"\ndone: {len(d)} ligands, {(time.time()-t0)/3600:.2f} h")
print(f"  dG_total : median {d.dG_total.median():.1f}, "
      f"range {d.dG_total.min():.1f} to {d.dG_total.max():.1f} kcal/mol")
print(f"  dE_vdw   : median {d.dE_vdw.median():.1f}")
print(f"  dE_elec  : median {d.dE_elec.median():.1f}")
print(f"written: {out_csv.name}")
