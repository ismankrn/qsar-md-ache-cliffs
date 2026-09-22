#!/usr/bin/env python
"""Prepare coordinate files for the huperzine A self- vs cross-docking panels.

Also writes Tyr337 from 4EY5 transformed into the 4EY7 frame, so the rotamer
shift that blocks the true huperzine A pose can be drawn directly.
"""
import warnings, sys
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate
sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_docking import canon, _ca
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
S, D = ROOT/"structures", ROOT/"docking"
SMI = "C/C=C/1\\[C@@H]2CC3=C([C@]1(CC(=C2)C)N)C=CC(=O)N3"

ref7 = mda.Universe(str(S/"4ey7.pdb"))
mob5 = mda.Universe(str(S/"4ey5.pdb"))

# transformasi 4EY5 -> 4EY7
r_ca, m_ca = _ca(ref7), _ca(mob5)
common = np.intersect1d(r_ca.residues.resids, m_ca.residues.resids)
rp = r_ca.positions[np.isin(r_ca.resids, common)]
mp = m_ca.positions[np.isin(m_ca.resids, common)]
n = min(len(rp), len(mp)); rp, mp = rp[:n], mp[:n]
rc, mc = rp.mean(0), mp.mean(0)
R, ca_rmsd = mda_align.rotation_matrix(mp-mc, rp-rc)
print(f"4EY5 -> 4EY7 superposition: {n} CA, RMSD {ca_rmsd:.3f} A")

def write_xyz(sel, xyz, path):
    tmp = sel.universe.copy().atoms[sel.indices]
    tmp.positions = xyz
    tmp.write(str(path))

# --- huperzine A kristal: kerangka asli 4EY5, dan hasil transformasi ke 4EY7
hup5 = mob5.select_atoms("resname HUP and chainID A and not name H*")
hup5.write(str(D/"panel_HUP_xtal_4ey5.pdb"))
write_xyz(hup5, (hup5.positions-mc)@R.T+rc, D/"panel_HUP_xtal_in4ey7.pdb")

# --- Tyr337: milik 4EY7, dan milik 4EY5 yang ditransformasi ke kerangka 4EY7
y7 = ref7.select_atoms("chainID A and resid 337 and not name H*")
y5 = mob5.select_atoms("chainID A and resid 337 and not name H*")
y7.write(str(D/"panel_Y337_4ey7.pdb"))
write_xyz(y5, (y5.positions-mc)@R.T+rc, D/"panel_Y337_4ey5_in4ey7.pdb")
shift = float(np.sqrt((((y7.positions) - ((y5.positions-mc)@R.T+rc))**2).sum(1).mean()))
print(f"Tyr337 side-chain shift 4EY5 vs 4EY7: {shift:.2f} A")

# --- pose hasil docking -> SDF
def pose_to_sdf(pdbqt, out):
    poses = RDKitMolCreate.from_pdbqt_mol(
        PDBQTMolecule.from_file(str(pdbqt), skip_typing=True))[0]
    one = Chem.Mol(poses); one.RemoveAllConformers()
    one.AddConformer(poses.GetConformer(0), assignId=True)
    one = Chem.RemoveHs(one)
    w = Chem.SDWriter(str(out)); w.write(one); w.close()
    return one

self_pose  = pose_to_sdf(D/"self_HUP_4ey5.pdbqt", D/"panel_HUP_self.sdf")
cross_pose = pose_to_sdf(D/"self_HUP_4ey7.pdbqt", D/"panel_HUP_cross.sdf")

# --- RMSD untuk pose yang benar-benar digambar
def ref_from_pdb(path):
    m = Chem.MolFromPDBFile(str(path), removeHs=True, sanitize=False)
    return canon(AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(SMI), m))

r_in7 = ref_from_pdb(D/"panel_HUP_xtal_in4ey7.pdb")
r_in5 = ref_from_pdb(D/"panel_HUP_xtal_4ey5.pdb")
print(f"\nRMSD of the poses drawn in the panels:")
print(f"  self-docking  (into 4EY5): {rdMolAlign.CalcRMS(canon(self_pose),  r_in5):.2f} A")
print(f"  cross-docking (into 4EY7): {rdMolAlign.CalcRMS(canon(cross_pose), r_in7):.2f} A")
