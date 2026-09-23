#!/usr/bin/env python
"""Rung 2 — geometric descriptors from the static docked pose.

The raw docking score already failed: it gave cliff pairs (|dp| = 2.66) and
control pairs (|dp| = 0.13) the same score difference, 0.52 vs 0.55 kcal/mol.
But a score is one number summarising a whole pose. This module goes back to
the pose itself and extracts what that number threw away.

The descriptors are chosen to test a specific mechanism, not to cast a wide
net. AChE's active site sits at the base of a narrow gorge with the catalytic
anionic site (Trp86) at the bottom and the peripheral anionic site (Trp286) at
the mouth, 15.2 A apart. Half of these cliffs are homologation cliffs, where
the only difference is linker length. The hypothesis is that linker length
decides whether a ligand can engage both sites at once.

So the headline descriptors are bridging geometry and cation-pi contacts --
every ligand here carries +1 or +2, and AChE binds cations through cation-pi
with the two tryptophans.

Descriptors are computed on the top pose of each of the three seeds and
reduced to the median, with the spread kept so unstable ones can be filtered.
"""
from __future__ import annotations
import warnings, sys
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
import MDAnalysis as mda
from rdkit import Chem, RDLogger
from meeko import PDBQTMolecule, RDKitMolCreate

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES, DOCK = ROOT / "results", ROOT / "docking"
POSES = DOCK / "poses"
TARGET, N_SEEDS = "CHEMBL220", 3
CONTACT = 4.5          # A, heavy-atom contact cutoff
CATPI = 6.0            # A, cation-pi is a long-range interaction

# gorge-lining residues, by their role
SITES = {
    "cas":   86,    # Trp86  - catalytic anionic site, gorge floor
    "pas":   286,   # Trp286 - peripheral anionic site, gorge mouth
    "gate":  337,   # Tyr337 - the swinging gate, shifts 2.53 A between structures
}
LINING = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]

# ---- receptor, loaded once -------------------------------------------
u = mda.Universe(str(DOCK / "4ey7_rec_clean.pdb"))
prot = u.select_atoms("protein and not name H*")
site_sel = {k: u.select_atoms(f"resid {r} and not backbone and not name H*")
            for k, r in SITES.items()}
# tryptophan ring centroids, for cation-pi
ring = {k: u.select_atoms(
            f"resid {SITES[k]} and name CD2 CE2 CE3 CZ2 CZ3 CH2").center_of_geometry()
        for k in ("cas", "pas")}
lining_sel = {r: u.select_atoms(f"resid {r} and not name H*") for r in LINING}


def pose_descriptors(xyz: np.ndarray, charges: np.ndarray) -> dict:
    """Geometry of one pose. xyz: ligand heavy atoms. charges: formal charges."""
    d = {}
    # distance from the ligand to each key site
    for k, sel in site_sel.items():
        dm = np.linalg.norm(xyz[:, None] - sel.positions[None], axis=2)
        d[f"d_{k}"] = float(dm.min())

    # --- the hypothesis: does the ligand engage BOTH anionic sites?
    d["bridges"] = int(d["d_cas"] <= CONTACT and d["d_pas"] <= CONTACT)
    # continuous version: the worse of the two contacts. Small means the
    # ligand reaches both; large means it only manages one.
    d["bridge_gap"] = float(max(d["d_cas"], d["d_pas"]))

    # --- cation-pi: these ligands are all cationic, and AChE binds cations
    # through the face of Trp86 and Trp286
    cat = xyz[charges > 0]
    for k in ("cas", "pas"):
        d[f"catpi_{k}"] = (float(np.linalg.norm(cat - ring[k], axis=1).min())
                           if len(cat) else np.nan)
    d["catpi_both"] = int(d["catpi_cas"] <= CATPI and d["catpi_pas"] <= CATPI)

    # --- burial and contact density
    dm = np.linalg.norm(xyz[:, None] - prot.positions[None], axis=2)
    d["n_contacts_4A"] = int((dm.min(axis=1) <= 4.0).sum())
    d["n_prot_5A"] = int((dm <= 5.0).any(axis=0).sum())
    d["n_prot_8A"] = int((dm <= 8.0).any(axis=0).sum())
    d["frac_lig_buried"] = float((dm.min(axis=1) <= 4.5).mean())

    # --- shape of the bound conformation. end_to_end is the descriptor that
    # should separate homologation cliffs: a longer linker can span further.
    d["rgyr"] = float(np.sqrt(((xyz - xyz.mean(0))**2).sum(1).mean()))
    intra = np.linalg.norm(xyz[:, None] - xyz[None], axis=2)
    d["end_to_end"] = float(intra.max())
    d["n_heavy"] = len(xyz)

    # --- per-residue contact counts along the gorge
    for r, sel in lining_sel.items():
        dm_r = np.linalg.norm(xyz[:, None] - sel.positions[None], axis=2)
        d[f"c_{r}"] = int((dm_r <= CONTACT).sum())
    return d


# ---- run over every ligand -------------------------------------------
score = pd.read_csv(RES / "07_docking_scores.csv")
ligands = sorted(score.molecule_chembl_id.unique())
print(f"ligands: {len(ligands)}   seeds: {N_SEEDS}\n", flush=True)

rows, failed = [], []
for i, cid in enumerate(ligands, 1):
    per_seed = []
    for seed in range(1, N_SEEDS + 1):
        f = POSES / f"{cid}_seed{seed}.pdbqt"
        if not f.exists():
            continue
        try:
            pm = PDBQTMolecule.from_file(str(f), skip_typing=True)
            mol = RDKitMolCreate.from_pdbqt_mol(pm)[0]
            mol = Chem.RemoveHs(mol)
            conf = mol.GetConformer(0)          # top-ranked pose of this seed
            xyz = np.array(conf.GetPositions())
            q = np.array([a.GetFormalCharge() for a in mol.GetAtoms()])
            per_seed.append(pose_descriptors(xyz, q))
        except Exception as e:
            failed.append((cid, seed, f"{type(e).__name__}: {str(e)[:50]}"))
    if not per_seed:
        continue
    df = pd.DataFrame(per_seed)
    rec = {"molecule_chembl_id": cid, "n_seeds_used": len(per_seed)}
    rec.update({c: df[c].median() for c in df.columns})
    # spread of the two headline descriptors, to flag unstable poses
    rec["spread_bridge_gap"] = float(df.bridge_gap.max() - df.bridge_gap.min())
    rec["spread_d_cas"] = float(df.d_cas.max() - df.d_cas.min())
    rows.append(rec)
    if i % 50 == 0:
        print(f"  {i}/{len(ligands)}", flush=True)

out = pd.DataFrame(rows)
out.to_csv(RES / "08_rung2_descriptors.csv", index=False)
print(f"\ndescriptors: {out.shape[0]} ligands x {out.shape[1]-2} features")
if failed:
    print(f"failed poses: {len(failed)}")
    for f in failed[:5]:
        print("   ", f)
print(f"\nbridging (both sites within {CONTACT} A): "
      f"{int(out.bridges.sum())} of {len(out)} ligands")
print(f"cation-pi to both tryptophans within {CATPI} A: "
      f"{int(out.catpi_both.sum())} of {len(out)}")
print(f"pose stability: median bridge_gap spread across seeds "
      f"{out.spread_bridge_gap.median():.2f} A")
print(f"\nwritten: 08_rung2_descriptors.csv")
