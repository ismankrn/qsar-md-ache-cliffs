#!/usr/bin/env python
"""Rebuild the control pairs with molecular size matched to the cliff pairs.

The original controls were matched on Tanimoto similarity and potency, but
not on size. That left a confound: cliff partners differ by 4.4 heavy atoms
on average, controls by only 1.8. Any size-correlated descriptor therefore
appeared to discriminate cliffs from controls for the wrong reason -- buried
surface area looked highly significant (p < 0.0001) purely because it tracks
heavy-atom difference at rho = +0.70.

That flaw affects every rung, including rung 1, so it is fixed before any
further comparison is made.

Cliff pairs are untouched; only the controls are reselected. Candidates whose
ligands have already been docked are preferred, so the extra compute is kept
small. The two largest size bins are thin in the candidate pool (35 available
against 37 needed, and 3 against 9), so matching there is as close as the
data allows and the shortfall is reported.
"""
from __future__ import annotations
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import BulkTanimotoSimilarity

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED = "CHEMBL220", 20260902
rng = np.random.default_rng(SEED)
SIM_CUTOFF, CTRL_MAX_DELTA = 0.90, 0.30

pairs = pd.read_csv(RES/f"02_pairs_{TARGET}.csv")
old_lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
pool = pd.read_csv(RES/f"01_compounds_{TARGET}.csv")
already_docked = set(pd.read_csv(RES/"07_docking_scores.csv").molecule_chembl_id)

cliff = pairs[pairs.group == "cliff"].reset_index(drop=True)
cliff_ids = set(cliff.ligand_a) | set(cliff.ligand_b)
print(f"cliff pairs kept unchanged : {len(cliff)}  ({len(cliff_ids)} ligands)")

mols = [Chem.MolFromSmiles(s) for s in pool.smiles_std]
keep = [i for i, m in enumerate(mols) if m is not None]
pool, mols = pool.iloc[keep].reset_index(drop=True), [mols[i] for i in keep]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [gen.GetFingerprint(m) for m in mols]
nh = np.array([m.GetNumHeavyAtoms() for m in mols])
flat = [Chem.MolToSmiles(m, isomericSmiles=False) for m in mols]
pot = pool.pchembl.to_numpy()
cid = pool.molecule_chembl_id.tolist()

nh_map = dict(zip(cid, nh))
d_cliff = np.array([abs(nh_map.get(a, 0) - nh_map.get(b, 0))
                    for a, b in zip(cliff.ligand_a, cliff.ligand_b)])
mean_p_cliff = ((cliff.pchembl_a + cliff.pchembl_b) / 2).to_numpy()

print("\nenumerating control candidates ...")
cand = []
for i in range(len(fps) - 1):
    if cid[i] in cliff_ids:
        continue
    sims = np.array(BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
    dp = np.abs(pot[i + 1:] - pot[i])
    hit = np.where((sims >= SIM_CUTOFF) & (dp <= CTRL_MAX_DELTA))[0]
    for k in hit:
        j = i + 1 + k
        if cid[j] in cliff_ids or flat[i] == flat[j]:
            continue
        cand.append({"ligand_a": cid[i], "ligand_b": cid[j],
                     "smiles_a": pool.smiles_std[i], "smiles_b": pool.smiles_std[j],
                     "pchembl_a": pot[i], "pchembl_b": pot[j],
                     "delta": float(dp[k]), "tanimoto": float(sims[k]),
                     "d_heavy": int(abs(nh[i] - nh[j])),
                     "mean_p": float((pot[i] + pot[j]) / 2),
                     "n_docked": int((cid[i] in already_docked)
                                     + (cid[j] in already_docked))})
cand = pd.DataFrame(cand)
print(f"candidates: {len(cand):,}")

# match the cliff joint distribution of (size difference, potency)
SZ = [0, 2, 4, 6, 9, 100]
PB = [0, 5, 6, 7, 8, 20]
cand["_s"] = pd.cut(cand.d_heavy, SZ, labels=False, include_lowest=True)
cand["_p"] = pd.cut(cand.mean_p, PB, labels=False)
want = pd.DataFrame({"_s": pd.cut(d_cliff, SZ, labels=False, include_lowest=True),
                     "_p": pd.cut(mean_p_cliff, PB, labels=False)}
                    ).groupby(["_s", "_p"], observed=True).size()

picked, used, short = [], set(), 0
for (s, p), n in want.items():
    block = cand[(cand._s == s) & (cand._p == p)]
    block = block[~block.ligand_a.isin(used) & ~block.ligand_b.isin(used)]
    # prefer pairs already docked, so the extra compute stays small
    block = block.sort_values("n_docked", ascending=False, kind="stable")
    took = 0
    for _, r in block.iterrows():
        if took >= n:
            break
        if r.ligand_a in used or r.ligand_b in used:
            continue
        picked.append(r); used.update([r.ligand_a, r.ligand_b]); took += 1
    if took < n:
        short += n - took
ctrl = pd.DataFrame(picked).drop(columns=["_s", "_p"], errors="ignore")
print(f"controls selected: {len(ctrl)}  (short of target by {short})")

d_ctrl = ctrl.d_heavy.to_numpy()
print(f"\nSIZE MATCHING")
print(f"  cliff pairs   : mean {d_cliff.mean():.1f}, median {np.median(d_cliff):.0f} heavy atoms")
print(f"  controls, old : mean 1.8")
print(f"  controls, new : mean {d_ctrl.mean():.1f}, median {np.median(d_ctrl):.0f}")
new_lig = (set(ctrl.ligand_a) | set(ctrl.ligand_b)) - already_docked
print(f"\nligands needing new docking + MM/GBSA: {len(new_lig)}")

# ---- write the rebuilt pair and ligand tables ------------------------
rows = []
for k, (_, r) in enumerate(ctrl.iterrows()):
    rows.append({"pair_id": f"control_{k:03d}", "group": "control",
                 "ligand_a": r.ligand_a, "ligand_b": r.ligand_b,
                 "pchembl_a": r.pchembl_a, "pchembl_b": r.pchembl_b,
                 "delta": r.delta, "tanimoto": r.tanimoto})
new_pairs = pd.concat([cliff[["pair_id","group","ligand_a","ligand_b",
                              "pchembl_a","pchembl_b","delta","tanimoto"]],
                       pd.DataFrame(rows)], ignore_index=True)
new_pairs.to_csv(RES/f"02_pairs_{TARGET}.csv", index=False)

keep_ids = set(new_pairs.ligand_a) | set(new_pairs.ligand_b)
old = old_lig.set_index("molecule_chembl_id")
lig_rows = []
for c in sorted(keep_ids):
    grp = "cliff" if c in cliff_ids else "control"
    if c in old.index:
        r = old.loc[c]
        lig_rows.append({"molecule_chembl_id": c, "smiles_std": r.smiles_std,
                         "pchembl": r.pchembl, "group": grp,
                         "smiles_md": r.smiles_md, "formal_charge": r.formal_charge,
                         "prot_ambiguous": r.prot_ambiguous})
    else:
        row = pool[pool.molecule_chembl_id == c].iloc[0]
        lig_rows.append({"molecule_chembl_id": c, "smiles_std": row.smiles_std,
                         "pchembl": row.pchembl, "group": grp,
                         "smiles_md": None, "formal_charge": None,
                         "prot_ambiguous": None})
pd.DataFrame(lig_rows).to_csv(RES/f"02_ligand_set_{TARGET}.csv", index=False)
print(f"\nwritten: 02_pairs_{TARGET}.csv  ({len(new_pairs)} pairs)")
print(f"         02_ligand_set_{TARGET}.csv  ({len(lig_rows)} ligands)")
print(f"\nNEXT: protonate the new ligands, re-run make_splits.py, then dock them")
