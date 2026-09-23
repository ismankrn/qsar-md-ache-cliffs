#!/usr/bin/env python
"""Rung 2 evaluation — the same comparison rung 1 received, on the same folds.

Mirrors rung1_benchmark.py exactly: identical frozen folds, identical three
framings, identical algorithm families, preprocessing fitted inside each
training fold. Anything that differs between the two rungs is therefore the
descriptors, which is the whole point of the ladder.
"""
from __future__ import annotations
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
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

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED, NF = "CHEMBL220", 20260902, 5

d2 = pd.read_csv(RES/"08_rung2_descriptors.csv")
br = pd.read_csv(RES/"08_rung2_bridging.csv")
desc = d2.merge(br, on="molecule_chembl_id", how="left").set_index("molecule_chembl_id")
desc = desc.drop(columns=[c for c in desc.columns if c.startswith("spread_")
                          or c == "n_seeds_used"], errors="ignore")

sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id","pchembl"]]
lig = lig.merge(sl[["molecule_chembl_id","split","cv_fold"]], on="molecule_chembl_id")
dev = lig[lig.split=="dev"].reset_index(drop=True)
dev = dev[dev.molecule_chembl_id.isin(desc.index)].reset_index(drop=True)
dpair = pairs[pairs.split=="dev"].reset_index(drop=True)

X = desc.loc[dev.molecule_chembl_id].to_numpy(float)
y = dev.pchembl.to_numpy(); fold = dev.cv_fold.to_numpy()
idx = {c:i for i,c in enumerate(dev.molecule_chembl_id)}
print(f"rung 2: {len(dev)} dev ligands x {X.shape[1]} descriptors")
print(f"        {len(dpair)} dev pairs\n")

def pipe(est, scale=True):
    st = [("imp", SimpleImputer(strategy="median")), ("var", VarianceThreshold(0.0))]
    if scale: st.append(("sc", StandardScaler()))
    st.append(("est", est)); return Pipeline(st)

REG = {
 "Ridge (linear)":    lambda: pipe(Ridge(alpha=10.0)),
 "PLS":               lambda: pipe(PLSRegression(n_components=5)),
 "k-NN (k=5)":        lambda: pipe(KNeighborsRegressor(n_neighbors=5)),
 "SVM (RBF)":         lambda: pipe(SVR(C=10.0)),
 "Random Forest":     lambda: pipe(RandomForestRegressor(n_estimators=500,
                          random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting": lambda: pipe(HistGradientBoostingRegressor(random_state=SEED), False),
}
CLF = {
 "Logistic regression": lambda: pipe(LogisticRegression(max_iter=3000)),
 "k-NN (k=5)":          lambda: pipe(KNeighborsClassifier(n_neighbors=5)),
 "SVM (RBF)":           lambda: pipe(SVC(C=10.0, probability=True, random_state=SEED)),
 "Random Forest":       lambda: pipe(RandomForestClassifier(n_estimators=500,
                            random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting":   lambda: pipe(HistGradientBoostingClassifier(random_state=SEED), False),
}

def cv(Xa, ya, fa, make, clf=False):
    out = np.zeros(len(ya))
    for k in range(NF):
        tr, te = fa != k, fa == k
        if te.sum() == 0: continue
        m = make(); m.fit(Xa[tr], ya[tr])
        out[te] = (m.predict_proba(Xa[te])[:,1] if clf
                   else np.asarray(m.predict(Xa[te])).ravel())
    return out

def pair_xy(mode):
    xs, ys, fs, gs = [], [], [], []
    for _, r in dpair.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx: continue
        da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        pa, pb = y[idx[r.ligand_a]], y[idx[r.ligand_b]]
        if mode == "anti":
            xs += [da-db, db-da]; ys += [pa-pb, pb-pa]
            fs += [r.cv_fold]*2;  gs += [r.group]*2
        else:
            xs.append(np.concatenate([da+db, np.abs(da-db)]))
            ys.append(1 if r.group=="cliff" else 0)
            fs.append(r.cv_fold); gs.append(r.group)
    return np.array(xs), np.array(ys), np.array(fs), np.array(gs)

rows = []
print("="*72); print("FRAMING A — per-ligand regression")
for a, mk in REG.items():
    p = cv(X, y, fold, mk)
    pm = dict(zip(dev.molecule_chembl_id, p))
    cl = dpair[(dpair.group=="cliff") & dpair.ligand_a.isin(pm) & dpair.ligand_b.isin(pm)]
    dp = np.abs([pm[x]-pm[z] for x,z in zip(cl.ligand_a, cl.ligand_b)])
    rows.append(dict(framing="A", algorithm=a, r2=r2_score(y,p),
        rho=spearmanr(y,p).statistic, cliff_true=cl.delta.mean(), cliff_pred=dp.mean()))
    print(f"  {a:20} R2={rows[-1]['r2']:6.3f}  cliff |dp| pred {dp.mean():.2f} "
          f"(true {cl.delta.mean():.2f})")

print("\n"+"="*72); print("FRAMING B — per-pair regression  [PRIMARY]")
Xp, yp, fp, gp = pair_xy("anti")
for a, mk in REG.items():
    p = cv(Xp, yp, fp, mk); m = gp=="cliff"
    rows.append(dict(framing="B", algorithm=a, r2=r2_score(yp,p),
        rho=spearmanr(yp,p).statistic, mae=mean_absolute_error(yp,p),
        cliff_true=np.abs(yp[m]).mean(), cliff_pred=np.abs(p[m]).mean(),
        ctrl_pred=np.abs(p[~m]).mean()))
    print(f"  {a:20} R2={rows[-1]['r2']:6.3f}  rho={rows[-1]['rho']:+.3f}  "
          f"cliff |dp| pred {np.abs(p[m]).mean():.2f} (true {np.abs(yp[m]).mean():.2f})")

print("\n"+"="*72); print("FRAMING C — per-pair classification")
Xc, yc, fc, _ = pair_xy("sym")
for a, mk in CLF.items():
    pr = cv(Xc, yc, fc, mk, clf=True)
    rows.append(dict(framing="C", algorithm=a, auc_roc=roc_auc_score(yc,pr),
        auc_pr=average_precision_score(yc,pr),
        mcc=matthews_corrcoef(yc,(pr>=.5).astype(int)),
        bal_acc=balanced_accuracy_score(yc,(pr>=.5).astype(int))))
    print(f"  {a:20} AUC={rows[-1]['auc_roc']:.3f}  MCC={rows[-1]['mcc']:+.3f}")

df = pd.DataFrame(rows); df.to_csv(RES/"08_rung2_benchmark.csv", index=False)

print("\n"+"="*72); print("RUNG 1 vs RUNG 2")
r1 = pd.read_csv(RES/"05_rung1_benchmark.csv")
for fr, metric, lab in [("A","r2","R2"), ("B","r2","R2"), ("C","auc_roc","AUC")]:
    a = r1[r1.framing==fr][metric].max(); b = df[df.framing==fr][metric].max()
    print(f"  framing {fr}  best {lab}:  rung 1 = {a:6.3f}   rung 2 = {b:6.3f}"
          f"   {'improved' if b>a else 'no gain'}")
print(f"\nwritten: 08_rung2_benchmark.csv")
