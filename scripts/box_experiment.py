#!/usr/bin/env python
"""Uji hipotesis: kotak 30 A terlalu longgar untuk ligan kecil.

Docking ulang HUP dan GNT pada beberapa ukuran kotak, pusat tetap di tengah
gorge (CAS-PAS). Kalau RMSD membaik saat kotak mengecil, diagnosisnya benar
dan solusinya kotak adaptif menurut ukuran ligan.
"""
import warnings, sys
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, Descriptors
from vina import Vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
sys.path.insert(0, str(Path(__file__).resolve().parent))
from protonation import protonate
from validate_docking import _ca, canon           # pakai ulang helper yang sudah teruji

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
STRUCT, DOCK, RES = ROOT/"structures", ROOT/"docking", ROOT/"results"
REF = "4ey7"

# ---------- bagian 1: sebaran ukuran ligan kampanye ----------
a = pd.read_csv(RES/"02_ligand_set_CHEMBL220.csv")
b = pd.read_csv(RES/"04_validation_ligands_CHEMBL220.csv")
lig = pd.concat([a[["molecule_chembl_id","smiles_md"]],
                 b[["molecule_chembl_id","smiles_md"]]]).drop_duplicates("molecule_chembl_id")
mw, dia = [], []
for s in lig.smiles_md:
    m0 = Chem.MolFromSmiles(s); mw.append(Descriptors.MolWt(m0))
    m = Chem.AddHs(m0); p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(m, p) == 0:
        AllChem.MMFFOptimizeMolecule(m, maxIters=300)
        c = m.GetConformer().GetPositions()
        dia.append(float(np.max(np.linalg.norm(c[:,None]-c[None],axis=2))))
    else:
        dia.append(np.nan)
mw, dia = np.array(mw), np.array(dia)
print(f"SEBARAN {len(mw)} LIGAN KAMPANYE")
for lo,hi,lab in [(0,300,"<300 Da  (seukuran HUP/GNT yang gagal)"),
                  (300,400,"300-400  (~donepezil, LULUS)"),
                  (400,550,"400-550  (~territrem B, LULUS)"),
                  (550,999,">550     (bis-ligan besar)")]:
    n=int(((mw>=lo)&(mw<hi)).sum()); print(f"  {lab:40} {n:>4}  ({100*n/len(mw):>4.1f}%)")
print(f"  median BM {np.median(mw):.0f} Da | diameter median {np.nanmedian(dia):.1f} A")
print(f"  ligan berdiameter <=20 A (muat kotak 22 A): "
      f"{100*np.nanmean(dia<=20):.0f}%\n")

# ---------- bagian 2: percobaan ukuran kotak ----------
TESTS = [("HUP","4ey5","huperzine A","C/C=C/1\\[C@@H]2CC3=C([C@]1(CC(=C2)C)N)C=CC(=O)N3"),
         ("GNT","4ey6","galantamin","C[N@@]1CC[C@@]23C=C[C@@H](C[C@@H]2Oc4c3c(ccc4OC)C1)O"),
         ("E20","4ey7","donepezil","COc1cc2c(cc1OC)C(=O)[C@@H](C2)CC3CCN(CC3)Cc4ccccc4")]

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
    tmp = DOCK/f"ref_{comp}.pdb"; l.write(str(tmp))
    m = Chem.MolFromPDBFile(str(tmp), removeHs=True, sanitize=False)
    m = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smi), m)
    cf = m.GetConformer()
    for k in range(m.GetNumAtoms()): cf.SetAtomPosition(k, xyz[k].tolist())
    return m

center = np.load(DOCK/"box_center.npy")
print("PERCOBAAN UKURAN KOTAK (pusat tetap di tengah gorge, exhaustiveness 32)")
print(f"{'ligan':>5} {'kotak':>7} {'skor':>8} {'RMSD-1':>8} {'RMSD-min':>9}")
rows=[]
for comp, pdb, label, smi in TESTS:
    rm = canon(ref_mol(comp, pdb, smi))
    p = protonate(smi)
    mol = Chem.AddHs(Chem.MolFromSmiles(p["smiles"]))
    par = AllChem.ETKDGv3(); par.randomSeed = 42
    AllChem.EmbedMolecule(mol, par); AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    setups = MoleculePreparation().prepare(mol)
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    f = DOCK/f"box_{comp}.pdbqt"; f.write_text(s)
    for size in (18, 22, 26, 30):
        v = Vina(sf_name="vina", cpu=0, verbosity=0)
        v.set_receptor(str(DOCK/f"{REF}_rec.pdbqt"))
        v.compute_vina_maps(center=[float(c) for c in center], box_size=[size]*3)
        v.set_ligand_from_file(str(f))
        try:
            v.dock(exhaustiveness=32, n_poses=10)
        except Exception as e:
            print(f"{comp:>5} {size:>6}A  gagal: {str(e)[:40]}"); continue
        o = DOCK/f"box_{comp}_{size}.pdbqt"; v.write_poses(str(o), n_poses=10, overwrite=True)
        en = v.energies(n_poses=10)
        poses = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(str(o), skip_typing=True))[0]
        rs=[]
        for ci in range(poses.GetNumConformers()):
            one = Chem.Mol(poses); one.RemoveAllConformers()
            one.AddConformer(poses.GetConformer(ci), assignId=True)
            try: rs.append(rdMolAlign.CalcRMS(canon(one), rm))
            except Exception: rs.append(np.nan)
        rs=np.array(rs,float)
        print(f"{comp:>5} {size:>6}A {float(en[0][0]):>8.2f} {rs[0]:>8.2f} {np.nanmin(rs):>9.2f}")
        rows.append({"ligan":comp,"kotak":size,"skor":float(en[0][0]),
                     "rmsd_1":float(rs[0]),"rmsd_min":float(np.nanmin(rs))})
    print()
pd.DataFrame(rows).to_csv(RES/"06_box_experiment.csv", index=False)
print("tersimpan: results/box_experiment.csv")
