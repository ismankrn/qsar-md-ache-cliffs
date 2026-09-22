#!/usr/bin/env python
"""Validasi redocking dengan replika berseed - angka siap kutip untuk manuskrip.

Vina bersifat stokastik: GNT pernah menghasilkan 5,89 A dan 1,81 A pada
pengaturan identik. Angka tunggal karena itu tidak memadai; setiap ligan
divalidasi dengan 5 seed berbeda dan dilaporkan sebagai mean +/- SD.
"""
import warnings, sys
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from vina import Vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
sys.path.insert(0, str(Path(__file__).resolve().parent))
from protonation import protonate
from validate_docking import canon, _ca

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
STRUCT, DOCK, RES = ROOT/"structures", ROOT/"docking", ROOT/"results"
REF, SEEDS, EXH = "4ey7", [1, 2, 3, 4, 5], 32

CASES = [("E20","4ey7","donepezil","self",
          "COc1cc2c(cc1OC)C(=O)[C@@H](C2)CC3CCN(CC3)Cc4ccccc4"),
         ("1YK","4m0f","territrem B","cross",
          "C[C@@]12CC[C@@]3([C@@]([C@]1(CC4=C(O2)C=C(OC4=O)c5cc(c(c(c5)OC)OC)OC)O)(C(=O)C=CC3(C)C)C)O"),
         ("GNT","4ey6","galantamin","cross",
          "C[N@@]1CC[C@@]23C=C[C@@H](C[C@@H]2Oc4c3c(ccc4OC)C1)O"),
         ("HUP","4ey5","huperzine A","cross",
          "C/C=C/1\\[C@@H]2CC3=C([C@]1(CC(=C2)C)N)C=CC(=O)N3")]

def ref_mol(comp, pdb, smi):
    mob = mda.Universe(str(STRUCT/f"{pdb}.pdb"))
    l = mob.select_atoms(f"resname {comp} and chainID A and not name H*")
    xyz = l.positions.copy()
    if pdb != REF:
        ref = mda.Universe(str(STRUCT/f"{REF}.pdb"))
        r_ca, m_ca = _ca(ref), _ca(mob)
        common = np.intersect1d(r_ca.residues.resids, m_ca.residues.resids)
        rp = r_ca.positions[np.isin(r_ca.resids, common)]
        mp = m_ca.positions[np.isin(m_ca.resids, common)]
        n = min(len(rp), len(mp)); rp, mp = rp[:n], mp[:n]
        rc, mc = rp.mean(0), mp.mean(0)
        R, _ = mda_align.rotation_matrix(mp-mc, rp-rc)
        xyz = (xyz-mc) @ R.T + rc
    t = DOCK/f"ref_{comp}.pdb"; l.write(str(t))
    m = Chem.MolFromPDBFile(str(t), removeHs=True, sanitize=False)
    m = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smi), m)
    cf = m.GetConformer()
    for k in range(m.GetNumAtoms()): cf.SetAtomPosition(k, xyz[k].tolist())
    return m

center = np.load(DOCK/f"{REF}_box_center.npy")
rows = []
print(f"reseptor {REF}, kotak 30 A, exhaustiveness {EXH}, {len(SEEDS)} seed\n")
print(f"{'ligan':>5} {'jenis':>6} {'RMSD-1 per seed':>34} {'mean+/-SD':>13} {'sukses':>7}")
print("-"*72)
for comp, pdb, label, kind, smi in CASES:
    rm = canon(ref_mol(comp, pdb, smi))
    p = protonate(smi)
    mol = Chem.AddHs(Chem.MolFromSmiles(p["smiles"]))
    par = AllChem.ETKDGv3(); par.randomSeed = 42
    AllChem.EmbedMolecule(mol, par); AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    s, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol)[0])
    f = DOCK/f"rep_{comp}.pdbqt"; f.write_text(s)
    r1s = []
    for sd in SEEDS:
        v = Vina(sf_name="vina", cpu=0, verbosity=0, seed=sd)
        v.set_receptor(str(DOCK/f"{REF}_rec.pdbqt"))
        v.compute_vina_maps(center=[float(c) for c in center], box_size=[30]*3)
        v.set_ligand_from_file(str(f))
        v.dock(exhaustiveness=EXH, n_poses=10)
        o = DOCK/f"rep_{comp}_{sd}.pdbqt"; v.write_poses(str(o), n_poses=10, overwrite=True)
        poses = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(str(o), skip_typing=True))[0]
        one = Chem.Mol(poses); one.RemoveAllConformers()
        one.AddConformer(poses.GetConformer(0), assignId=True)
        try: r1s.append(float(rdMolAlign.CalcRMS(canon(one), rm)))
        except Exception: r1s.append(np.nan)
    a = np.array(r1s, float)
    succ = int(np.sum(a < 2.0))
    print(f"{comp:>5} {kind:>6} {'  '.join(f'{x:5.2f}' for x in a):>34}"
          f" {np.nanmean(a):6.2f}+/-{np.nanstd(a):4.2f} {succ}/{len(SEEDS):>5}")
    rows.append({"ligan":comp,"nama":label,"jenis":kind,
                 "rmsd_mean":float(np.nanmean(a)),"rmsd_sd":float(np.nanstd(a)),
                 "sukses":succ,"n_seed":len(SEEDS),
                 **{f"seed_{sd}":x for sd,x in zip(SEEDS,a)}})
pd.DataFrame(rows).to_csv(RES/"06_redocking_validation.csv", index=False)
print("\ntersimpan: results/redocking_validation.csv")
