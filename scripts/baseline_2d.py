#!/usr/bin/env python
"""Rung 1: the 2D-QSAR baseline, and the go/no-go test for the study premise.

If a 2D model already predicts the cliff pairs, there is nothing for MD to
explain and the project should stop here. Costs seconds, so it runs first.

Two baselines are fitted, because a strawman baseline proves nothing:
  binary ECFP4          - the literature-standard fingerprint
  count ECFP4 + RDKit2D - a deliberately STRONGER, fairer baseline

Evaluation is pair-aware: ligands share pairs, so folds are built from
connected components of the ligand-pair graph. Splitting naively would put
one partner of a cliff in train and the other in test, and the model would
simply memorise the scaffold.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
TARGET = "CHEMBL220"
SEED = 20260902

lig = pd.read_csv(RESULTS / f"02_ligand_set_{TARGET}.csv")
pairs = pd.read_csv(RESULTS / f"02_pairs_{TARGET}.csv")
mols = [Chem.MolFromSmiles(s) for s in lig.smiles_std]
y = lig.pchembl.to_numpy()
print(f"ligands {len(lig)} | pairs {len(pairs)} "
      f"(cliff {(pairs.group=='cliff').sum()}, control {(pairs.group=='control').sum()})")

# ---- features -----------------------------------------------------
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
X_bit = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype=float)
X_cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], dtype=float)

names = [d[0] for d in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
X_rdk = np.array([calc.CalcDescriptors(m) for m in mols], dtype=float)
X_rdk = np.nan_to_num(X_rdk, nan=0.0, posinf=0.0, neginf=0.0)
keep = X_rdk.std(axis=0) > 0
X_rdk = X_rdk[:, keep]
X_strong = np.hstack([X_cnt, X_rdk])
print(f"features: binary ECFP {X_bit.shape[1]} | count ECFP+2D {X_strong.shape[1]}")

# ---- pair-aware folds: connected components of the ligand-pair graph
idx = {c: i for i, c in enumerate(lig.molecule_chembl_id)}
parent = list(range(len(lig)))
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb
for _, r in pairs.iterrows():
    if r.ligand_a in idx and r.ligand_b in idx:
        union(idx[r.ligand_a], idx[r.ligand_b])
groups = np.array([find(i) for i in range(len(lig))])
print(f"connected components (fold groups): {len(set(groups))}")

# ---- out-of-fold predictions --------------------------------------
def oof(X):
    pred = np.zeros(len(y))
    gkf = GroupKFold(n_splits=5)
    for tr, te in gkf.split(X, y, groups):
        m = RandomForestRegressor(n_estimators=500, random_state=SEED,
                                  n_jobs=-1, min_samples_leaf=1)
        m.fit(X[tr], y[tr])
        pred[te] = m.predict(X[te])
    return pred

results = {}
for label, X in (("binary ECFP4", X_bit), ("count ECFP4 + RDKit2D", X_strong)):
    p = oof(X)
    results[label] = p
    print(f"\n--- {label} ---")
    print(f"  global  R2={r2_score(y,p):.3f}  RMSE={mean_squared_error(y,p)**0.5:.3f}"
          f"  Spearman={spearmanr(y,p).statistic:.3f}")

    pmap = dict(zip(lig.molecule_chembl_id, p))
    rows = []
    for _, r in pairs.iterrows():
        if r.ligand_a in pmap and r.ligand_b in pmap:
            rows.append({"group": r.group, "actual": r.delta,
                         "pred": abs(pmap[r.ligand_a] - pmap[r.ligand_b])})
    d = pd.DataFrame(rows)
    for g in ("cliff", "control"):
        sub = d[d.group == g]
        print(f"  {g:8} pairs n={len(sub):3}  actual |dp| {sub.actual.mean():.2f}"
              f"   PREDICTED |dp| {sub.pred.mean():.2f}"
              f"  (recovers {100*sub.pred.mean()/max(sub.actual.mean(),1e-9):.0f}%)")
    cl = d[d.group == "cliff"]
    rho = spearmanr(cl.actual, cl.pred).statistic
    print(f"  correlation of predicted vs actual |dp| on cliffs: rho={rho:.3f}")

# ---- the guaranteed-failure subset --------------------------------
bitgen = gen
fp = {c: bitgen.GetFingerprint(Chem.MolFromSmiles(s))
      for c, s in zip(lig.molecule_chembl_id, lig.smiles_std)}
from rdkit.DataStructs import TanimotoSimilarity
ident = [r.pair_id for _, r in pairs.iterrows()
         if r.group == "cliff" and r.ligand_a in fp and r.ligand_b in fp
         and TanimotoSimilarity(fp[r.ligand_a], fp[r.ligand_b]) >= 0.9999]
print(f"\ncliff pairs where binary ECFP4 is IDENTICAL: {len(ident)}")
print("  a fingerprint-only model must predict |dp| = 0 for these, by construction")

out = lig[["molecule_chembl_id", "pchembl", "group"]].copy()
for k, v in results.items():
    out[f"pred_{k.split()[0]}"] = v
out.to_csv(RESULTS / "05_baseline_2d_predictions.csv", index=False)
print(f"\nwritten: baseline_2d_predictions.csv")
