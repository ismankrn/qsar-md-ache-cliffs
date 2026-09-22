#!/usr/bin/env python
"""Rung 1 of the descriptor ladder: the complete 2D-QSAR baseline.

Supersedes baseline_2d.py, which was provisional: it computed its own folds
instead of loading the frozen split, and filtered zero-variance descriptors
across the whole dataset before splitting (a small leak).

Implements all three framings on IDENTICAL frozen folds:
  A  per-ligand regression        y = pChEMBL
  B  per-pair regression          y = dpChEMBL (signed, antisymmetric)
  C  per-pair classification      y = cliff / non-cliff (symmetric)

Every preprocessing step is fitted INSIDE the training fold. The sealed test
set is never touched here; all numbers are cross-validation on the dev set.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import (r2_score, mean_squared_error, mean_absolute_error,
                             roc_auc_score, average_precision_score,
                             matthews_corrcoef, balanced_accuracy_score)
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
TARGET, SEED, NF = "CHEMBL220", 20260902, 5

lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
lig = lig.merge(sl[["molecule_chembl_id","split","cv_fold"]], on="molecule_chembl_id")

dev = lig[lig.split == "dev"].reset_index(drop=True)
dpair = pairs[pairs.split == "dev"].reset_index(drop=True)
print(f"dev: {len(dev)} ligands, {len(dpair)} pairs "
      f"(cliff {(dpair.group=='cliff').sum()}, control {(dpair.group=='control').sum()})")
print(f"test set sealed: {(lig.split=='test').sum()} ligands — not touched here\n")

# ---------------- descriptors ----------------
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
names = [d[0] for d in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)

def featurise(smiles):
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    bit = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], float)
    cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], float)
    rdk = np.nan_to_num(np.array([calc.CalcDescriptors(m) for m in mols], float),
                        nan=0.0, posinf=0.0, neginf=0.0)
    return {"binary ECFP4": bit, "count ECFP4 + RDKit2D": np.hstack([cnt, rdk])}

FEAT = featurise(dev.smiles_std)
idx = {c: i for i, c in enumerate(dev.molecule_chembl_id)}
fold = dev.cv_fold.to_numpy()
y_lig = dev.pchembl.to_numpy()

def clean_fit(Xtr, Xte):
    """Variance filter fitted on TRAIN ONLY — never on the full dataset."""
    keep = Xtr.std(axis=0) > 0
    return Xtr[:, keep], Xte[:, keep]

def rf_reg():
    return RandomForestRegressor(n_estimators=500, random_state=SEED, n_jobs=-1)
def rf_clf():
    return RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1)

# ---------------- Framing A ----------------
print("=" * 74)
print("FRAMING A — per-ligand regression, y = pChEMBL")
resA = {}
for name, X in FEAT.items():
    pred = np.zeros(len(y_lig))
    for k in range(NF):
        tr, te = fold != k, fold == k
        Xtr, Xte = clean_fit(X[tr], X[te])
        m = rf_reg(); m.fit(Xtr, y_lig[tr]); pred[te] = m.predict(Xte)
    resA[name] = pred
    print(f"  {name:24} R2={r2_score(y_lig,pred):6.3f}  "
          f"RMSE={mean_squared_error(y_lig,pred)**0.5:5.3f}  "
          f"rho={spearmanr(y_lig,pred).statistic:6.3f}")
    pm = dict(zip(dev.molecule_chembl_id, pred))
    sub = dpair[dpair.ligand_a.isin(pm) & dpair.ligand_b.isin(pm)]
    for g in ("cliff", "control"):
        s = sub[sub.group == g]
        pd_ = np.abs([pm[a]-pm[b] for a, b in zip(s.ligand_a, s.ligand_b)])
        print(f"      {g:8} n={len(s):3}  true |dp|={s.delta.mean():.2f}"
              f"  predicted |dp|={pd_.mean():.2f}")

# ---------------- pair matrices ----------------
def pair_xy(tab, X, mode):
    """mode 'anti' -> signed difference (framing B); 'sym' -> symmetric (C)."""
    xs, ys, fs = [], [], []
    for _, r in tab.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx:
            continue
        da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        pa = dev.pchembl[idx[r.ligand_a]]; pb = dev.pchembl[idx[r.ligand_b]]
        if mode == "anti":
            # label antisymmetric -> features antisymmetric, both orderings
            xs += [da-db, db-da]; ys += [pa-pb, pb-pa]; fs += [r.cv_fold]*2
        else:
            # label symmetric -> features must be symmetric too
            xs.append(np.concatenate([da+db, np.abs(da-db)]))
            ys.append(1 if r.group == "cliff" else 0); fs.append(r.cv_fold)
    return np.array(xs), np.array(ys), np.array(fs)

# ---------------- Framing B ----------------
print("\n" + "=" * 74)
print("FRAMING B — per-pair regression, y = signed dpChEMBL  [PRIMARY ENDPOINT]")
for name, X in FEAT.items():
    Xp, yp, fp = pair_xy(dpair, X, "anti")
    pred = np.zeros(len(yp))
    for k in range(NF):
        tr, te = fp != k, fp == k
        if te.sum() == 0: continue
        Xtr, Xte = clean_fit(Xp[tr], Xp[te])
        m = rf_reg(); m.fit(Xtr, yp[tr]); pred[te] = m.predict(Xte)
    grp = np.repeat((dpair.group == "cliff").to_numpy(), 2)[:len(yp)]
    print(f"  {name:24} R2={r2_score(yp,pred):6.3f}  "
          f"MAE={mean_absolute_error(yp,pred):5.3f}  "
          f"rho={spearmanr(yp,pred).statistic:6.3f}   (n={len(yp)} rows)")
    for g, mask in (("cliff", grp), ("control", ~grp)):
        print(f"      {g:8} true |dp|={np.abs(yp[mask]).mean():.2f}"
              f"  predicted |dp|={np.abs(pred[mask]).mean():.2f}"
              f"  MAE={mean_absolute_error(yp[mask],pred[mask]):.2f}")

# ---------------- Framing C ----------------
print("\n" + "=" * 74)
print("FRAMING C — per-pair classification, y = cliff / non-cliff")
for name, X in FEAT.items():
    Xp, yp, fp = pair_xy(dpair, X, "sym")
    prob = np.zeros(len(yp))
    for k in range(NF):
        tr, te = fp != k, fp == k
        if te.sum() == 0: continue
        Xtr, Xte = clean_fit(Xp[tr], Xp[te])
        m = rf_clf(); m.fit(Xtr, yp[tr]); prob[te] = m.predict_proba(Xte)[:, 1]
    pred = (prob >= 0.5).astype(int)
    print(f"  {name:24} AUC-ROC={roc_auc_score(yp,prob):.3f}  "
          f"AUC-PR={average_precision_score(yp,prob):.3f}  "
          f"MCC={matthews_corrcoef(yp,pred):.3f}  "
          f"BalAcc={balanced_accuracy_score(yp,pred):.3f}")

# ---------------- proof: identical fingerprints ----------------
print("\n" + "=" * 74)
fp_map = {c: gen.GetFingerprint(Chem.MolFromSmiles(s))
          for c, s in zip(dev.molecule_chembl_id, dev.smiles_std)}
from rdkit.DataStructs import TanimotoSimilarity
ident = [(r.pair_id, r.delta) for _, r in dpair.iterrows()
         if r.group == "cliff" and r.ligand_a in fp_map and r.ligand_b in fp_map
         and TanimotoSimilarity(fp_map[r.ligand_a], fp_map[r.ligand_b]) >= 0.9999]
print(f"Dev cliff pairs whose binary ECFP4 vectors are IDENTICAL: "
      f"{len(ident)}/{(dpair.group=='cliff').sum()}")
if ident:
    d = np.array([x[1] for x in ident])
    print(f"  their true |dp| ranges {d.min():.2f}-{d.max():.2f} "
          f"({10**d.min():,.0f}x to {10**d.max():,.0f}x in potency)")
    print("  a fingerprint-only model must predict dp = 0 for every one of them")
