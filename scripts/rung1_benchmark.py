#!/usr/bin/env python
"""Rung 1 benchmark: many algorithms, so the baseline cannot be called weak.

The study's central claim is that 2D descriptors cannot predict activity
cliffs. A single untuned Random Forest is not enough to support that -- a
reviewer would object that the failure is the model's, not the descriptors'.
Six algorithm families are therefore tested per framing, each inside a
Pipeline so that scaling and variance filtering are fitted on training folds
only. k-NN is included deliberately: it IS the similarity principle made
explicit, so its behaviour is conceptually informative.
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
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.cross_decomposition import PLSRegression
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.svm import SVR, SVC
from sklearn.ensemble import (RandomForestRegressor, RandomForestClassifier,
                              HistGradientBoostingRegressor,
                              HistGradientBoostingClassifier)
from sklearn.metrics import (r2_score, mean_squared_error, mean_absolute_error,
                             roc_auc_score, average_precision_score,
                             matthews_corrcoef, balanced_accuracy_score)
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED, NF = "CHEMBL220", 20260902, 5

lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
lig = lig.merge(sl[["molecule_chembl_id","split","cv_fold"]], on="molecule_chembl_id")
dev = lig[lig.split=="dev"].reset_index(drop=True)
dpair = pairs[pairs.split=="dev"].reset_index(drop=True)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
names = [d[0] for d in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
mols = [Chem.MolFromSmiles(s) for s in dev.smiles_std]
bit = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], float)
cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], float)
rdk = np.nan_to_num(np.array([calc.CalcDescriptors(m) for m in mols], float),
                    nan=0.0, posinf=0.0, neginf=0.0)
FEAT = {"binary ECFP4": bit, "count ECFP4 + RDKit2D": np.hstack([cnt, rdk])}
idx = {c:i for i,c in enumerate(dev.molecule_chembl_id)}
fold = dev.cv_fold.to_numpy(); y_lig = dev.pchembl.to_numpy()

def pipe(est, scale=True):
    steps = [("var", VarianceThreshold(0.0))]
    if scale: steps.append(("sc", StandardScaler()))
    steps.append(("est", est))
    return Pipeline(steps)

REG = {
 "Ridge (linear)":        lambda: pipe(Ridge(alpha=10.0)),
 "PLS":                   lambda: pipe(PLSRegression(n_components=5)),
 "k-NN (k=5)":            lambda: pipe(KNeighborsRegressor(n_neighbors=5)),
 "SVM (RBF)":             lambda: pipe(SVR(C=10.0, gamma="scale")),
 "Random Forest":         lambda: pipe(RandomForestRegressor(
                              n_estimators=500, random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting":     lambda: pipe(HistGradientBoostingRegressor(
                              random_state=SEED), False),
}
CLF = {
 "Logistic regression":   lambda: pipe(LogisticRegression(max_iter=3000, C=1.0)),
 "k-NN (k=5)":            lambda: pipe(KNeighborsClassifier(n_neighbors=5)),
 "SVM (RBF)":             lambda: pipe(SVC(C=10.0, probability=True, random_state=SEED)),
 "Random Forest":         lambda: pipe(RandomForestClassifier(
                              n_estimators=500, random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting":     lambda: pipe(HistGradientBoostingClassifier(
                              random_state=SEED), False),
}

def cv_reg(X, y, f, make):
    pred = np.zeros(len(y))
    for k in range(NF):
        tr, te = f != k, f == k
        if te.sum() == 0: continue
        m = make(); m.fit(X[tr], y[tr])
        p = m.predict(X[te])
        pred[te] = np.asarray(p).ravel()
    return pred

def pair_xy(tab, X, mode):
    xs, ys, fs, gs = [], [], [], []
    for _, r in tab.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx: continue
        da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        pa, pb = dev.pchembl[idx[r.ligand_a]], dev.pchembl[idx[r.ligand_b]]
        if mode == "anti":
            xs += [da-db, db-da]; ys += [pa-pb, pb-pa]
            fs += [r.cv_fold]*2;  gs += [r.group]*2
        else:
            xs.append(np.concatenate([da+db, np.abs(da-db)]))
            ys.append(1 if r.group=="cliff" else 0)
            fs.append(r.cv_fold); gs.append(r.group)
    return np.array(xs), np.array(ys), np.array(fs), np.array(gs)

rows = []
for fname, X in FEAT.items():
    # ---- Framing A
    for aname, make in REG.items():
        p = cv_reg(X, y_lig, fold, make)
        pm = dict(zip(dev.molecule_chembl_id, p))
        sub = dpair[dpair.ligand_a.isin(pm) & dpair.ligand_b.isin(pm)]
        cl = sub[sub.group=="cliff"]
        dp = np.abs([pm[a]-pm[b] for a,b in zip(cl.ligand_a, cl.ligand_b)])
        rows.append(dict(framing="A", features=fname, algorithm=aname,
            r2=r2_score(y_lig,p), rmse=mean_squared_error(y_lig,p)**0.5,
            rho=spearmanr(y_lig,p).statistic,
            cliff_true=cl.delta.mean(), cliff_pred=dp.mean()))
        print(f"A {fname:22} {aname:20} R2={rows[-1]['r2']:6.3f}", flush=True)
    # ---- Framing B
    Xp, yp, fp, gp = pair_xy(dpair, X, "anti")
    for aname, make in REG.items():
        p = cv_reg(Xp, yp, fp, make)
        m = gp == "cliff"
        rows.append(dict(framing="B", features=fname, algorithm=aname,
            r2=r2_score(yp,p), mae=mean_absolute_error(yp,p),
            rho=spearmanr(yp,p).statistic,
            cliff_true=np.abs(yp[m]).mean(), cliff_pred=np.abs(p[m]).mean(),
            cliff_mae=mean_absolute_error(yp[m],p[m]),
            ctrl_true=np.abs(yp[~m]).mean(), ctrl_pred=np.abs(p[~m]).mean()))
        print(f"B {fname:22} {aname:20} R2={rows[-1]['r2']:6.3f}", flush=True)
    # ---- Framing C
    Xc, yc, fc, _ = pair_xy(dpair, X, "sym")
    for aname, make in CLF.items():
        prob = np.zeros(len(yc))
        for k in range(NF):
            tr, te = fc != k, fc == k
            if te.sum() == 0: continue
            mdl = make(); mdl.fit(Xc[tr], yc[tr])
            prob[te] = mdl.predict_proba(Xc[te])[:,1]
        pr = (prob >= 0.5).astype(int)
        rows.append(dict(framing="C", features=fname, algorithm=aname,
            auc_roc=roc_auc_score(yc,prob), auc_pr=average_precision_score(yc,prob),
            mcc=matthews_corrcoef(yc,pr), bal_acc=balanced_accuracy_score(yc,pr)))
        print(f"C {fname:22} {aname:20} AUC={rows[-1]['auc_roc']:.3f}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(RES/"05_rung1_benchmark.csv", index=False)
print(f"\nwritten: rung1_benchmark.csv  ({len(df)} model runs)")
