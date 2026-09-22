#!/usr/bin/env python
"""Robustness test: does the 2D baseline fail because it was starved of data?

The rung-1 models saw 217 ligands. A reviewer can object that this is simply
too small for QSAR. This script retrains the same 2D models on the FULL
filtered AChE pool -- roughly 25x more compounds -- with every curated ligand
excluded, then asks the same question of the same cliff pairs.

Two numbers come out:
  1. Pool-internal performance on a random split. This is the reference against
     published AChE QSAR studies: it shows whether the baseline is competitive
     or weak.
  2. Cliff Δp recovery on the curated development pairs, predicted by the
     pool-trained model. This is the actual test.

A representational failure is immune to more data; a statistical one is not.
"""
from __future__ import annotations
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from scipy.stats import spearmanr
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_ligand_set import md_ok

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED = "CHEMBL220", 20260902
rng = np.random.default_rng(SEED)

pool = pd.read_csv(RES/f"01_compounds_{TARGET}.csv")
cur  = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
val  = pd.read_csv(RES/f"04_validation_ligands_{TARGET}.csv")
sl   = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
dpair = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
dpair = dpair[dpair.split == "dev"].reset_index(drop=True)

# MD-feasible pool, then remove every curated / validation ligand (no leakage)
pool = pool[[c[0] for c in pool.smiles_std.map(md_ok)]].reset_index(drop=True)
excl = set(cur.molecule_chembl_id) | set(val.molecule_chembl_id)
train_pool = pool[~pool.molecule_chembl_id.isin(excl)].reset_index(drop=True)
print(f"MD-feasible pool          : {len(pool):,}")
print(f"excluded (curated + val)  : {len(excl):,}")
print(f"pool available for training: {len(train_pool):,}   "
      f"({len(train_pool)/217:.0f}x the rung-1 training set)\n")

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
names = [d[0] for d in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)

def feats(smiles):
    ms = [Chem.MolFromSmiles(s) for s in smiles]
    bit = np.array([gen.GetFingerprintAsNumPy(m) for m in ms], float)
    cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in ms], float)
    rdk = np.nan_to_num(np.array([calc.CalcDescriptors(m) for m in ms], float),
                        nan=0.0, posinf=0.0, neginf=0.0)
    return {"binary ECFP4": bit, "count ECFP4 + RDKit2D": np.hstack([cnt, rdk])}

print("featurising ...", flush=True)
Ftr = feats(train_pool.smiles_std)
dev = cur.merge(sl[["molecule_chembl_id","split"]], on="molecule_chembl_id")
dev = dev[dev.split == "dev"].reset_index(drop=True)
Fdev = feats(dev.smiles_std)
y_tr, y_dev = train_pool.pchembl.to_numpy(), dev.pchembl.to_numpy()

def pipe(est, scale=True):
    st = [("var", VarianceThreshold(0.0))]
    if scale: st.append(("sc", StandardScaler()))
    st.append(("est", est)); return Pipeline(st)

ALGOS = {
 "PLS":               lambda: pipe(PLSRegression(n_components=5)),
 "Random Forest":     lambda: pipe(RandomForestRegressor(
                          n_estimators=400, random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting": lambda: pipe(HistGradientBoostingRegressor(
                          random_state=SEED), False),
}

rows = []
for fname, Xtr in Ftr.items():
    Xdev = Fdev[fname]
    for aname, make in ALGOS.items():
        # ---- 1. pool-internal 5-fold CV: the literature reference point
        oof = np.zeros(len(y_tr))
        for tr, te in KFold(5, shuffle=True, random_state=SEED).split(Xtr):
            m = make(); m.fit(Xtr[tr], y_tr[tr])
            oof[te] = np.asarray(m.predict(Xtr[te])).ravel()
        pool_r2 = r2_score(y_tr, oof)
        pool_rho = spearmanr(y_tr, oof).statistic
        pool_rmse = mean_squared_error(y_tr, oof) ** 0.5

        # ---- 2. train on the whole pool, predict the curated dev ligands
        m = make(); m.fit(Xtr, y_tr)
        pred = np.asarray(m.predict(Xdev)).ravel()
        pm = dict(zip(dev.molecule_chembl_id, pred))
        sub = dpair[dpair.ligand_a.isin(pm) & dpair.ligand_b.isin(pm)]
        cl = sub[sub.group == "cliff"]; ct = sub[sub.group == "control"]
        dcl = np.abs([pm[a]-pm[b] for a, b in zip(cl.ligand_a, cl.ligand_b)])
        dct = np.abs([pm[a]-pm[b] for a, b in zip(ct.ligand_a, ct.ligand_b)])
        # signed delta correlation - can the model rank cliff magnitude at all?
        signed_true = np.array([r.delta for _, r in cl.iterrows()])
        rho_cliff = spearmanr(signed_true, dcl).statistic

        rows.append(dict(features=fname, algorithm=aname,
            pool_r2=pool_r2, pool_rmse=pool_rmse, pool_rho=pool_rho,
            dev_r2=r2_score(y_dev, pred),
            cliff_true=cl.delta.mean(), cliff_pred=dcl.mean(),
            ctrl_true=ct.delta.mean(), ctrl_pred=dct.mean(),
            rho_cliff_rank=rho_cliff))
        print(f"{fname:22} {aname:18} pool R2={pool_r2:6.3f} | "
              f"cliff pred |dp|={dcl.mean():.2f} (true {cl.delta.mean():.2f}) | "
              f"rank rho={rho_cliff:6.3f}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(RES/"05_rung1_robustness.csv", index=False)

print("\n" + "="*78)
print("1. POOL-INTERNAL PERFORMANCE  (random split, the literature reference)")
print(f"   best R2 = {df.pool_r2.max():.3f}  "
      f"({df.loc[df.pool_r2.idxmax(),'algorithm']}, {df.loc[df.pool_r2.idxmax(),'features']})")
print(f"   range   = {df.pool_r2.min():.3f} to {df.pool_r2.max():.3f} over {len(df)} models")
print("\n2. CLIFF RECOVERY WITH ~25x MORE TRAINING DATA")
print(f"   true |dp| on cliff pairs      : {df.cliff_true.iloc[0]:.2f}")
print(f"   predicted |dp|, best model    : {df.cliff_pred.max():.2f} "
      f"({100*df.cliff_pred.max()/df.cliff_true.iloc[0]:.0f}% recovered)")
print(f"   predicted |dp|, worst model   : {df.cliff_pred.min():.2f}")
print(f"   rank correlation with true dp : {df.rho_cliff_rank.min():.3f} to {df.rho_cliff_rank.max():.3f}")
print("\n   rung-1 comparison (217 ligands): best recovery was 30%, best rank rho ~0.05")
print(f"\nwritten: rung1_robustness.csv")
