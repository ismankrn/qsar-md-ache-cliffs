#!/usr/bin/env python
"""Validasi protokol docking terhadap ligan ko-kristal AChE.

Semua 354 ligan akan didocking ke SATU reseptor (4EY7). Uji yang benar karena
itu bukan self-docking, melainkan CROSS-docking: ligan dari struktur lain
disuperposisi ke kerangka 4EY7, lalu didocking ke reseptor 4EY7 dan pose
hasilnya dibandingkan dengan pose kristalnya.

Ambang keberhasilan yang lazim: RMSD < 2,0 A terhadap pose kristal.
Ligan uji sengaja dipilih beragam mode pengikatan:
  E20 donepezil   - dual-site (CAS+PAS), 4EY7  -> self-docking
  HUP huperzine A - CAS saja, 4EY5             -> cross
  GNT galantamin  - CAS saja, 4EY6             -> cross
  1YK territrem B - besar (527 Da), 4M0F       -> cross
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from vina import Vina
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate

sys.path.insert(0, str(Path(__file__).resolve().parent))
from protonation import protonate

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
STRUCT, DOCK = ROOT / "structures", ROOT / "docking"
REF_PDB, EXHAUST, N_POSES = "4ey7", 32, 10

TESTS = [
    ("E20", "4ey7", "donepezil (dual-site)",
     "COc1cc2c(cc1OC)C(=O)[C@@H](C2)CC3CCN(CC3)Cc4ccccc4"),
    ("HUP", "4ey5", "huperzine A (CAS)",
     "C/C=C/1\\[C@@H]2CC3=C([C@]1(CC(=C2)C)N)C=CC(=O)N3"),
    ("GNT", "4ey6", "galantamin (CAS)",
     "C[N@@]1CC[C@@]23C=C[C@@H](C[C@@H]2Oc4c3c(ccc4OC)C1)O"),
    ("1YK", "4m0f", "territrem B (besar)",
     "C[C@@]12CC[C@@]3([C@@]([C@]1(CC4=C(O2)C=C(OC4=O)c5cc(c(c(c5)OC)OC)OC)O)"
     "(C(=O)C=CC3(C)C)C)O"),
]


def _ca(u):
    """CA rantai A, satu per residu (buang altloc ganda)."""
    sel = u.select_atoms("protein and chainID A and name CA")
    seen, idx = set(), []
    for i, at in enumerate(sel):
        if at.resid not in seen:
            seen.add(at.resid); idx.append(i)
    return sel[idx]


def ref_mol(comp: str, pdb: str, smi: str):
    """Molekul ligan kristal DENGAN koordinatnya, disuperposisi ke 4EY7.

    Dibangun dari file PDB (bukan dari SMILES) supaya koordinat menempel pada
    atom yang benar. Ikatan diambil dari template SMILES lewat
    AssignBondOrdersFromTemplate, yang mempertahankan urutan atom PDB.
    """
    mob = mda.Universe(str(STRUCT / f"{pdb}.pdb"))
    lig = mob.select_atoms(f"resname {comp} and chainID A and not name H*")
    xyz = lig.positions.copy()

    if pdb != REF_PDB:
        ref = mda.Universe(str(STRUCT / f"{REF_PDB}.pdb"))
        r_ca, m_ca = _ca(ref), _ca(mob)
        common = np.intersect1d(r_ca.residues.resids, m_ca.residues.resids)
        rmask = np.isin(r_ca.resids, common)
        mmask = np.isin(m_ca.resids, common)
        rp, mp = r_ca.positions[rmask], m_ca.positions[mmask]
        n = min(len(rp), len(mp)); rp, mp = rp[:n], mp[:n]
        rc, mc = rp.mean(axis=0), mp.mean(axis=0)
        R, rmsd = mda_align.rotation_matrix(mp - mc, rp - rc)
        print(f"     superposisi {pdb}->{REF_PDB}: {n} CA, RMSD {rmsd:.2f} A")
        xyz = (xyz - mc) @ R.T + rc

    tmp = DOCK / f"ref_{comp}.pdb"
    lig.write(str(tmp))
    m = Chem.MolFromPDBFile(str(tmp), removeHs=True, sanitize=False)
    if m is None:
        return None
    tpl = Chem.MolFromSmiles(smi)
    try:
        m = AllChem.AssignBondOrdersFromTemplate(tpl, m)
    except Exception as e:
        print(f"     peringatan: gagal terapkan template ({str(e)[:50]})")
        return None
    conf = m.GetConformer()
    for k in range(m.GetNumAtoms()):
        conf.SetAtomPosition(k, xyz[k].tolist())
    return m


def canon(mol):
    """Samakan graf untuk RMSD: buang H, netralkan muatan formal.

    Pose hasil docking berbentuk terprotonasi (+1), sedangkan referensi dibangun
    dari template SMILES netral. Tanpa penyeragaman ini CalcRMS menolak
    memetakan keduanya.
    """
    m = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
    for a in m.GetAtoms():
        a.SetFormalCharge(0); a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    return m


def main():
    v = Vina(sf_name="vina", cpu=0, verbosity=0)
    v.set_receptor(str(DOCK / f"{REF_PDB}_rec.pdbqt"))
    center = np.load(DOCK / "box_center.npy")
    v.compute_vina_maps(center=[float(c) for c in center], box_size=[30, 30, 30])
    print(f"reseptor {REF_PDB}, kotak 30 A di {np.round(center,2)}, exhaustiveness {EXHAUST}\n")

    rows = []
    for comp, pdb, label, smi in TESTS:
        print(f"[{comp}] {label}")
        rmol = ref_mol(comp, pdb, smi)
        if rmol is None:
            print("     referensi tidak terbangun, dilewati\n"); continue

        p = protonate(smi)
        mol = Chem.AddHs(Chem.MolFromSmiles(p["smiles"]))
        par = AllChem.ETKDGv3(); par.randomSeed = 42
        AllChem.EmbedMolecule(mol, par); AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        print(f"     muatan {p['charge']:+d}, aturan: {','.join(p['rules']) or '-'}")

        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            print(f"     GAGAL siapkan ligan: {err}"); continue
        lig_f = DOCK / f"val_{comp}.pdbqt"
        lig_f.write_text(pdbqt)

        v.set_ligand_from_file(str(lig_f))
        v.dock(exhaustiveness=EXHAUST, n_poses=N_POSES)
        out = DOCK / f"val_{comp}_out.pdbqt"
        v.write_poses(str(out), n_poses=N_POSES, overwrite=True)
        energies = v.energies(n_poses=N_POSES)

        pm = PDBQTMolecule.from_file(str(out), skip_typing=True)
        poses = RDKitMolCreate.from_pdbqt_mol(pm)[0]
        ref_c = canon(rmol)

        rmsds = []
        for ci in range(poses.GetNumConformers()):
            one = Chem.Mol(poses)
            one.RemoveAllConformers()
            one.AddConformer(poses.GetConformer(ci), assignId=True)
            try:
                rmsds.append(rdMolAlign.CalcRMS(canon(one), ref_c))
            except Exception:
                rmsds.append(float("nan"))
        rmsds = np.array(rmsds, dtype=float)
        best_i = int(np.nanargmin(rmsds)) if np.any(~np.isnan(rmsds)) else 0
        best_rmsd = float(rmsds[best_i])
        rmsd_top = float(rmsds[0])
        top = float(energies[0][0]) if len(energies) else float("nan")
        print(f"     skor {top:.2f} kcal/mol | RMSD pose-1 {rmsd_top:.2f} A"
              f" | RMSD terbaik {best_rmsd:.2f} A (pose #{best_i+1}/{len(rmsds)})")
        rows.append((comp, label, top, rmsd_top, best_rmsd))
        print()

    print("=" * 66)
    print(f"{'ligan':>5} {'skor':>8} {'RMSD-1':>8} {'RMSD-min':>9}   status")
    for comp, label, e, r1, rb in rows:
        st = "LULUS" if r1 < 2.0 else ("pose benar ada" if rb < 2.0 else "GAGAL")
        print(f"{comp:>5} {e:>8.2f} {r1:>8.2f} {rb:>9.2f}   {st}  {label}")
    print("\nRMSD-1 = pose peringkat teratas; RMSD-min = pose terbaik dari 10.")


if __name__ == "__main__":
    main()
