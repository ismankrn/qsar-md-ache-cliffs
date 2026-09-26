"""Stage 15 — rung-2 pose geometry on the ENSEMBLE poses.

Stage 11 showed that scoring against eight gorge conformations adds no predictive power:
correlation with pChEMBL was -0.327 for 4EY7 alone and -0.329 for the ensemble best, and
framing B fell from R2 = 0.027 to 0.008. But that compared *scores*. Rung-2 descriptors are
geometric, and two poses with the same score can sit differently in the gorge. This closes
that gap with data instead of recording it as a limitation.

Each receptor gets its own selections rebuilt -- CAS/PAS/gate sidechains, the Trp86 and
Trp286 ring centroids, every protein heavy atom, and the fourteen lining residues -- because
the whole point of an ensemble is that this geometry differs between structures. All eight
were verified beforehand to carry the right residues with the right identities and six-atom
rings; an empty selection would produce wrong descriptors without raising anything.

Descriptors are computed exactly as scripts/rung2_descriptors.py computes them, so the
numbers are comparable to 08_rung2_descriptors.csv. Aggregation across receptors is done
two ways, mirroring stage 11's ens_best and ens_mean:
  best  -- descriptors from whichever receptor scored the ligand lowest
  med   -- median of each descriptor across all eight receptors
"""
from __future__ import annotations
import re, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
import MDAnalysis as mda
from rdkit import Chem, RDLogger
from meeko import PDBQTMolecule, RDKitMolCreate
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RES, DOCK = ROOT/"results", ROOT/"docking"
ENS = DOCK/"ens_poses"
CONTACT, CATPI = 4.5, 6.0
SITES = {"cas": 86, "pas": 286, "gate": 337}
LINING = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]
RINGSEL = "name CD2 CE2 CE3 CZ2 CZ3 CH2"


def build_receptor(pdb: Path):
    """Selections for one receptor. Raises if anything the descriptors need is missing."""
    u = mda.Universe(str(pdb))
    prot = u.select_atoms("protein and not name H*")
    if prot.n_atoms == 0:
        raise RuntimeError(f"{pdb.name}: tidak ada atom protein")
    site = {}
    for k, r in SITES.items():
        s = u.select_atoms(f"resid {r} and not backbone and not name H*")
        if s.n_atoms == 0:
            raise RuntimeError(f"{pdb.name}: sidechain {k} (resid {r}) kosong")
        site[k] = s.positions.copy()
    ring = {}
    for k in ("cas", "pas"):
        s = u.select_atoms(f"resid {SITES[k]} and {RINGSEL}")
        if s.n_atoms != 6:
            raise RuntimeError(f"{pdb.name}: cincin {k} punya {s.n_atoms} atom, harus 6")
        ring[k] = s.center_of_geometry()
    lin = {}
    for r in LINING:
        s = u.select_atoms(f"resid {r} and not name H*")
        if s.n_atoms == 0:
            raise RuntimeError(f"{pdb.name}: residu pelapis {r} kosong")
        lin[r] = s.positions.copy()
    return {"prot": prot.positions.copy(), "site": site, "ring": ring, "lin": lin}


def pose_descriptors(xyz, charges, R):
    d = {}
    for k, pos in R["site"].items():
        d[f"d_{k}"] = float(np.linalg.norm(xyz[:, None] - pos[None], axis=2).min())
    d["bridges"] = int(d["d_cas"] <= CONTACT and d["d_pas"] <= CONTACT)
    d["bridge_gap"] = float(max(d["d_cas"], d["d_pas"]))
    cat = xyz[charges > 0]
    for k in ("cas", "pas"):
        d[f"catpi_{k}"] = (float(np.linalg.norm(cat - R["ring"][k], axis=1).min())
                           if len(cat) else np.nan)
    d["catpi_both"] = int(d["catpi_cas"] <= CATPI and d["catpi_pas"] <= CATPI)
    dm = np.linalg.norm(xyz[:, None] - R["prot"][None], axis=2)
    d["n_contacts_4A"] = int((dm.min(axis=1) <= 4.0).sum())
    d["n_prot_5A"] = int((dm <= 5.0).any(axis=0).sum())
    d["n_prot_8A"] = int((dm <= 8.0).any(axis=0).sum())
    d["frac_lig_buried"] = float((dm.min(axis=1) <= 4.5).mean())
    d["rgyr"] = float(np.sqrt(((xyz - xyz.mean(0))**2).sum(1).mean()))
    d["end_to_end"] = float(np.linalg.norm(xyz[:, None] - xyz[None], axis=2).max())
    d["n_heavy"] = len(xyz)
    for r, pos in R["lin"].items():
        d[f"c_{r}"] = int((np.linalg.norm(xyz[:, None] - pos[None], axis=2) <= CONTACT).sum())
    return d


recs = sorted(re.sub(r".*ens_|_clean\.pdb", "", str(p))
              for p in DOCK.glob("ens_*_clean.pdb"))
print(f"reseptor: {recs}", flush=True)
REC = {}
for r in recs:
    REC[r] = build_receptor(DOCK/f"ens_{r}_clean.pdb")
print(f"  seleksi terbangun untuk {len(REC)} reseptor, semua lolos pemeriksaan\n", flush=True)

score = pd.read_csv(RES/"11_ensemble_docking.csv")
ligands = sorted(score.molecule_chembl_id.unique())
print(f"ligan: {len(ligands)}  -> {len(ligands)*len(recs)} pose\n", flush=True)

rows, failed, t0 = [], [], time.time()
for i, cid in enumerate(ligands, 1):
    for r in recs:
        f = ENS/f"{cid}_{r}.pdbqt"
        if not f.exists():
            failed.append((cid, r, "berkas pose tidak ada")); continue
        try:
            pm = PDBQTMolecule.from_file(str(f), skip_typing=True)
            mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
            xyz = np.array(mol.GetConformer(0).GetPositions())   # pose peringkat teratas
            q = np.array([a.GetFormalCharge() for a in mol.GetAtoms()])
            d = pose_descriptors(xyz, q, REC[r])
            d.update({"molecule_chembl_id": cid, "receptor": r})
            rows.append(d)
        except Exception as e:
            failed.append((cid, r, f"{type(e).__name__}: {str(e)[:60]}"))
    if i % 40 == 0:
        el = time.time()-t0
        print(f"  {i}/{len(ligands)} ligan, {el:.0f} s, ETA {el/i*(len(ligands)-i):.0f} s",
              flush=True)

per_pose = pd.DataFrame(rows)
per_pose.to_csv(RES/"15_rung2_ensemble_per_pose.csv", index=False)
print(f"\nper-pose: {len(per_pose)} baris, gagal {len(failed)}", flush=True)
for x in failed[:5]:
    print("  ", x)

feat = [c for c in per_pose.columns if c not in ("molecule_chembl_id", "receptor")]
s1 = score.set_index(["molecule_chembl_id", "receptor"]).score_1
pp = per_pose.set_index(["molecule_chembl_id", "receptor"])
pp["score_1"] = s1
pp = pp.reset_index()

best = (pp.sort_values("score_1").groupby("molecule_chembl_id").first()
          .reset_index()[["molecule_chembl_id", "receptor"] + feat])
best = best.rename(columns={"receptor": "best_receptor"})
med = pp.groupby("molecule_chembl_id")[feat].median().reset_index()

best.to_csv(RES/"15_rung2_ensemble_best.csv", index=False)
med.to_csv(RES/"15_rung2_ensemble_median.csv", index=False)
print(f"agregat: best {best.shape}, median {med.shape}")
print(f"reseptor terpilih sebagai 'best': "
      f"{best.best_receptor.value_counts().to_dict()}")
