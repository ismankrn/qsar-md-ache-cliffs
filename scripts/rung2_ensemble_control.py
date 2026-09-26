"""Stage 15b — control: is the ensemble-median improvement denoising, or receptor geometry?

Ensemble-median rung-2 descriptors gave framing B R2 = 0.274 against 0.160 for 4EY7 alone,
while ensemble-BEST gave 0.018. That pattern is what noise reduction predicts, not what
"the ensemble finds better poses" predicts: if some receptor genuinely bound a ligand
better, taking that receptor's descriptors would win, not averaging all eight. Vina is
stochastic, so one pose is one noisy sample and a median over eight is simply a more
precise estimate of the same geometry.

The control separates the two. Compute descriptors for all eight ensemble poses but measure
every one of them against the 4EY7 receptor only, then take the median. Receptor
conformation is now held constant and only the eight-fold averaging remains:

  ~0.27  -> the gain was denoising; ensemble conformations add nothing
  ~0.16  -> the gain came from the receptor conformations themselves
"""
from __future__ import annotations
import time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda
from rdkit import Chem, RDLogger
from meeko import PDBQTMolecule, RDKitMolCreate
RDLogger.DisableLog("rdApp.*")
import importlib.util

ROOT = Path(__file__).resolve().parent.parent
RES, DOCK = ROOT/"results", ROOT/"docking"
ENS = DOCK/"ens_poses"

spec = importlib.util.spec_from_file_location("re2", ROOT/"scripts"/"rung2_ensemble.py")
# reuse the descriptor code without executing the script body
src = (ROOT/"scripts"/"rung2_ensemble.py").read_text()
head = src.split("recs = sorted(")[0]
ns: dict = {"__file__": str(ROOT/"scripts"/"rung2_ensemble.py"),
            "__name__": "rung2_ensemble_head"}
exec(compile(head, "rung2_ensemble_head", "exec"), ns)
build_receptor, pose_descriptors = ns["build_receptor"], ns["pose_descriptors"]

R4 = build_receptor(DOCK/"ens_4ey7_clean.pdb")     # the ONLY receptor used here
score = pd.read_csv(RES/"11_ensemble_docking.csv")
recs = sorted(score.receptor.unique())
ligands = sorted(score.molecule_chembl_id.unique())
print(f"reseptor pose: {recs}")
print(f"reseptor pengukur: 4ey7 saja (dikunci)")
print(f"ligan: {len(ligands)}\n", flush=True)

rows, failed, t0 = [], [], time.time()
for i, cid in enumerate(ligands, 1):
    for r in recs:
        f = ENS/f"{cid}_{r}.pdbqt"
        if not f.exists():
            failed.append((cid, r)); continue
        try:
            pm = PDBQTMolecule.from_file(str(f), skip_typing=True)
            mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm)[0])
            xyz = np.array(mol.GetConformer(0).GetPositions())
            q = np.array([a.GetFormalCharge() for a in mol.GetAtoms()])
            d = pose_descriptors(xyz, q, R4)        # <- always 4ey7
            d.update({"molecule_chembl_id": cid, "receptor": r})
            rows.append(d)
        except Exception as e:
            failed.append((cid, r))
    if i % 80 == 0:
        print(f"  {i}/{len(ligands)}  {time.time()-t0:.0f} s", flush=True)

pp = pd.DataFrame(rows)
feat = [c for c in pp.columns if c not in ("molecule_chembl_id", "receptor")]
med = pp.groupby("molecule_chembl_id")[feat].median().reset_index()
med.to_csv(RES/"15b_rung2_ensemble_median_fixed4ey7.csv", index=False)
print(f"\nper-pose {len(pp)} baris, gagal {len(failed)}; median {med.shape} tersimpan")
