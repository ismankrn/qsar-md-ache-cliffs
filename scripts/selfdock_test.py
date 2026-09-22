#!/usr/bin/env python
"""Uji penentu: self-docking vs cross-docking.

Kalau HUP dan GNT berhasil didocking ke reseptornya SENDIRI tetapi gagal ke
4EY7, penyebabnya kekakuan reseptor - rantai samping 4EY7 tertata untuk
donepezil - bukan cacat protokol.
"""
import warnings, sys
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from vina import Vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
sys.path.insert(0, str(Path(__file__).resolve().parent))
from protonation import protonate
from validate_docking import canon, _ca

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
DOCK = ROOT / "docking"

import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
STRUCT = ROOT / "structures"
REF = "4ey7"


def ref_mol(comp, pdb, smi):
    """Ligan kristal dengan koordinatnya, disuperposisi ke kerangka 4EY7."""
    mob = mda.Universe(str(STRUCT / f"{pdb}.pdb"))
    l = mob.select_atoms(f"resname {comp} and chainID A and not name H*")
    xyz = l.positions.copy()
    if pdb != REF:
        ref = mda.Universe(str(STRUCT / f"{REF}.pdb"))
        r_ca, m_ca = _ca(ref), _ca(mob)
        common = np.intersect1d(r_ca.residues.resids, m_ca.residues.resids)
        rp = r_ca.positions[np.isin(r_ca.resids, common)]
        mp = m_ca.positions[np.isin(m_ca.resids, common)]
        n = min(len(rp), len(mp)); rp, mp = rp[:n], mp[:n]
        rc, mc = rp.mean(0), mp.mean(0)
        R, _ = mda_align.rotation_matrix(mp - mc, rp - rc)
        xyz = (xyz - mc) @ R.T + rc
    tmp = DOCK / f"ref_{comp}.pdb"; l.write(str(tmp))
    m = Chem.MolFromPDBFile(str(tmp), removeHs=True, sanitize=False)
    m = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smi), m)
    cf = m.GetConformer()
    for k in range(m.GetNumAtoms()):
        cf.SetAtomPosition(k, xyz[k].tolist())
    return m


CASES = [
    ("HUP", "4ey5", "huperzine A", "C/C=C/1\\[C@@H]2CC3=C([C@]1(CC(=C2)C)N)C=CC(=O)N3"),
    ("GNT", "4ey6", "galantamin",  "C[N@@]1CC[C@@]23C=C[C@@H](C[C@@H]2Oc4c3c(ccc4OC)C1)O"),
]

print(f"{'ligan':>5} {'reseptor':>10} {'jenis':>7} {'skor':>8} {'RMSD-1':>8} {'RMSD-min':>9}")
print("-" * 56)
for comp, own, label, smi in CASES:
    rm = canon(ref_mol(comp, own, smi))          # referensi selalu di kerangka 4EY7
    p = protonate(smi)
    mol = Chem.AddHs(Chem.MolFromSmiles(p["smiles"]))
    par = AllChem.ETKDGv3(); par.randomSeed = 42
    AllChem.EmbedMolecule(mol, par); AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    s, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol)[0])
    f = DOCK / f"self_{comp}.pdbqt"; f.write_text(s)

    for rec, kind in ((own, "self"), ("4ey7", "cross")):
        # referensi kristal ada di kerangka 4EY7; struktur lain sudah
        # tersuperposisi ke sana (RMSD CA 0,15-0,21 A), jadi pose dari
        # reseptor mana pun bisa dibandingkan langsung dengan referensi itu
        center = np.load(DOCK / f"{rec}_box_center.npy")
        v = Vina(sf_name="vina", cpu=0, verbosity=0)
        v.set_receptor(str(DOCK / f"{rec}_rec.pdbqt"))
        v.compute_vina_maps(center=[float(c) for c in center], box_size=[30]*3)
        v.set_ligand_from_file(str(f))
        v.dock(exhaustiveness=32, n_poses=10)
        o = DOCK / f"self_{comp}_{rec}.pdbqt"
        v.write_poses(str(o), n_poses=10, overwrite=True)
        en = v.energies(n_poses=10)
        poses = RDKitMolCreate.from_pdbqt_mol(
            PDBQTMolecule.from_file(str(o), skip_typing=True))[0]
        rs = []
        for ci in range(poses.GetNumConformers()):
            one = Chem.Mol(poses); one.RemoveAllConformers()
            one.AddConformer(poses.GetConformer(ci), assignId=True)
            try: rs.append(rdMolAlign.CalcRMS(canon(one), rm))
            except Exception: rs.append(np.nan)
        rs = np.array(rs, float)
        print(f"{comp:>5} {rec:>10} {kind:>7} {float(en[0][0]):>8.2f} "
              f"{rs[0]:>8.2f} {np.nanmin(rs):>9.2f}")
    print()
