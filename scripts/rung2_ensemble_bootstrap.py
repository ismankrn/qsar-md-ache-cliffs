"""Stage 15c — is the ensemble gain real, or within noise on 145 pairs?

Framing B R2 came out 0.160 for 4EY7, 0.207 for eight poses measured against 4EY7, and
0.274 for eight poses each measured against its own receptor. Differences of 0.05-0.07 on
145 pairs with 5-fold CV are not obviously outside sampling error, and taking the best of
six algorithms per set biases every figure upward. This puts a confidence interval on them.

The algorithm is fixed to SVM (RBF) a priori -- it wins for three of the four sets, so
holding it constant removes the selection effect without cherry-picking. The frozen CV folds
are untouched; the bootstrap resamples the 145 evaluation PAIRS, keeping each pair's two
antisymmetric rows together so the resampling unit matches the unit of independence.
"""
from __future__ import annotations
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.svm import SVR
from sklearn.metrics import r2_score

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED, NF, NBOOT = "CHEMBL220", 20260902, 5, 4000

sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id", "pchembl"]]
lig = lig.merge(sl[["molecule_chembl_id", "split", "cv_fold"]], on="molecule_chembl_id")
dpair = pairs[pairs.split == "dev"].reset_index(drop=True)

SETS = {
 "4ey7":        (RES/"08_rung2_descriptors.csv", ["n_seeds_used"]),
 "best":        (RES/"15_rung2_ensemble_best.csv", ["best_receptor"]),
 "med_fix4ey7": (RES/"15b_rung2_ensemble_median_fixed4ey7.csv", []),
 "median":      (RES/"15_rung2_ensemble_median.csv", []),
}

def mk():
    return Pipeline([("imp", SimpleImputer(strategy="median")),
                     ("vt", VarianceThreshold()), ("sc", StandardScaler()),
                     ("est", SVR(C=10.0, gamma="scale"))])

preds, truth, pair_of = {}, None, None
for name, (path, drop) in SETS.items():
    d = pd.read_csv(path).set_index("molecule_chembl_id")
    d = d.drop(columns=[c for c in d.columns if c.startswith("spread_") or c in drop],
               errors="ignore").select_dtypes(include=[np.number])
    dev = lig[lig.split == "dev"]
    dev = dev[dev.molecule_chembl_id.isin(d.index)].reset_index(drop=True)
    X = d.loc[dev.molecule_chembl_id].to_numpy(float)
    y = dev.pchembl.to_numpy(); fold = dev.cv_fold.to_numpy()
    idx = {c: i for i, c in enumerate(dev.molecule_chembl_id)}
    xs, ys, fs, po = [], [], [], []
    for pi, r in dpair.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx: continue
        da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        pa, pb = y[idx[r.ligand_a]], y[idx[r.ligand_b]]
        xs += [da-db, db-da]; ys += [pa-pb, pb-pa]; fs += [r.cv_fold]*2; po += [pi]*2
    Xp, yp, fp = np.array(xs), np.array(ys), np.array(fs)
    out = np.zeros(len(yp))
    for k in range(NF):
        tr, te = fp != k, fp == k
        if te.sum() == 0: continue
        m = mk(); m.fit(Xp[tr], yp[tr]); out[te] = np.asarray(m.predict(Xp[te])).ravel()
    preds[name] = out
    if truth is None:
        truth, pair_of = yp, np.array(po)

upairs = np.unique(pair_of)
rows_by_pair = {p: np.where(pair_of == p)[0] for p in upairs}
rng = np.random.default_rng(SEED)

print(f"pasangan evaluasi: {len(upairs)}   baris: {len(truth)}   bootstrap: {NBOOT}\n")
boot = {n: np.empty(NBOOT) for n in SETS}
for b in range(NBOOT):
    samp = rng.choice(upairs, size=len(upairs), replace=True)
    ii = np.concatenate([rows_by_pair[p] for p in samp])
    for n in SETS:
        boot[n][b] = r2_score(truth[ii], preds[n][ii])

print(f"{'set':14s}{'R2':>8s}{'95% CI':>20s}")
for n in SETS:
    lo, hi = np.percentile(boot[n], [2.5, 97.5])
    print(f"  {n:12s}{r2_score(truth, preds[n]):8.3f}   [{lo:+.3f}, {hi:+.3f}]")

print(f"\nselisih berpasangan (bootstrap yang sama, jadi berkorelasi):")
def diff(a, b):
    d = boot[a]-boot[b]
    lo, hi = np.percentile(d, [2.5, 97.5])
    p = 2*min((d <= 0).mean(), (d >= 0).mean())
    return r2_score(truth, preds[a])-r2_score(truth, preds[b]), lo, hi, p
for a, b, label in (("med_fix4ey7", "4ey7", "denoising saja (8 sampel vs 3)"),
                    ("median", "med_fix4ey7", "geometri reseptor saja"),
                    ("median", "4ey7", "keduanya bersama"),
                    ("best", "4ey7", "pose berskor terbaik vs 4EY7")):
    dd, lo, hi, p = diff(a, b)
    sig = "SIGNIFIKAN" if (lo > 0 or hi < 0) else "dalam derau"
    print(f"  {label:34s} {dd:+.3f}  [{lo:+.3f}, {hi:+.3f}]  p={p:.3f}  {sig}")
