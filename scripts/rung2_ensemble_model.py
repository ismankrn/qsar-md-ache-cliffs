"""Stage 15 — does rung-2 geometry from ENSEMBLE poses beat geometry from 4EY7 alone?

Runs the identical model suite, framings, CV folds and preprocessing as
scripts/rung2_model.py over three descriptor sets, so the numbers are comparable to each
other by construction:

  4ey7   -- 08_rung2_descriptors.csv, the single rigid receptor
  best   -- 15_rung2_ensemble_best.csv, descriptors from whichever of eight receptors
            scored each ligand lowest
  median -- 15_rung2_ensemble_median.csv, median of each descriptor across all eight

The published rung-2 figure of R2 = 0.167 used these pose descriptors PLUS the
axis-projection bridging features in 08_rung2_bridging.csv, which have no ensemble
counterpart. The 4ey7 baseline here therefore excludes bridging, so all three sets carry
the same features and the comparison is matched rather than inherited.
"""
from __future__ import annotations
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.cross_decomposition import PLSRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED, NF = "CHEMBL220", 20260902, 5

sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id", "pchembl"]]
lig = lig.merge(sl[["molecule_chembl_id", "split", "cv_fold"]], on="molecule_chembl_id")
dpair = pairs[pairs.split == "dev"].reset_index(drop=True)


def pipe(est, scale=True):
    steps = [("imp", SimpleImputer(strategy="median")), ("vt", VarianceThreshold())]
    if scale: steps.append(("sc", StandardScaler()))
    steps.append(("est", est))
    return Pipeline(steps)

REG = {
 "Ridge (linear)":    lambda: pipe(Ridge(alpha=10.0)),
 "PLS":               lambda: pipe(PLSRegression(n_components=5)),
 "k-NN (k=5)":        lambda: pipe(KNeighborsRegressor(5)),
 "SVM (RBF)":         lambda: pipe(SVR(C=10.0, gamma="scale")),
 "Random Forest":     lambda: pipe(RandomForestRegressor(n_estimators=500,
                                    random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting": lambda: pipe(HistGradientBoostingRegressor(random_state=SEED), False),
}

SETS = {
 "4ey7":   (RES/"08_rung2_descriptors.csv",       ["n_seeds_used"]),
 "best":   (RES/"15_rung2_ensemble_best.csv",     ["best_receptor"]),
 "median": (RES/"15_rung2_ensemble_median.csv",   []),
 # control: same eight poses, but every one measured against 4EY7 only. Isolates the
 # eight-fold averaging from the receptor conformations.
 "med_fix4ey7": (RES/"15b_rung2_ensemble_median_fixed4ey7.csv", []),
}

def cv(Xa, ya, fa, make):
    out = np.zeros(len(ya))
    for k in range(NF):
        tr, te = fa != k, fa == k
        if te.sum() == 0: continue
        m = make(); m.fit(Xa[tr], ya[tr])
        out[te] = np.asarray(m.predict(Xa[te])).ravel()
    return out

rows = []
for name, (path, drop) in SETS.items():
    d = pd.read_csv(path).set_index("molecule_chembl_id")
    d = d.drop(columns=[c for c in d.columns
                        if c.startswith("spread_") or c in drop], errors="ignore")
    d = d.select_dtypes(include=[np.number])
    dev = lig[lig.split == "dev"]
    dev = dev[dev.molecule_chembl_id.isin(d.index)].reset_index(drop=True)
    X = d.loc[dev.molecule_chembl_id].to_numpy(float)
    y = dev.pchembl.to_numpy(); fold = dev.cv_fold.to_numpy()
    idx = {c: i for i, c in enumerate(dev.molecule_chembl_id)}

    xs, ys, fs, gs = [], [], [], []
    for _, r in dpair.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx: continue
        da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        pa, pb = y[idx[r.ligand_a]], y[idx[r.ligand_b]]
        xs += [da-db, db-da]; ys += [pa-pb, pb-pa]
        fs += [r.cv_fold]*2;  gs += [r.group]*2
    Xp, yp, fp, gp = np.array(xs), np.array(ys), np.array(fs), np.array(gs)

    print(f"\n{'='*70}\n{name}: {len(dev)} ligan dev x {X.shape[1]} deskriptor, "
          f"{len(dpair)} pasangan -> {len(yp)} baris antisimetris")
    print(f"{'  algoritma':22s}{'A: R2':>9s}{'B: R2':>9s}{'B: rho':>9s}{'B: MAE':>9s}")
    for a, mk in REG.items():
        pa_ = cv(X, y, fold, mk)
        r2a = r2_score(y, pa_)
        pb_ = cv(Xp, yp, fp, mk)
        r2b = r2_score(yp, pb_)
        rho = spearmanr(yp, pb_).statistic
        mae = mean_absolute_error(yp, pb_)
        print(f"  {a:20s}{r2a:9.3f}{r2b:9.3f}{rho:+9.3f}{mae:9.3f}")
        rows.append(dict(descriptors=name, algorithm=a, r2_A=r2a, r2_B=r2b,
                         rho_B=rho, mae_B=mae))

out = pd.DataFrame(rows)
out.to_csv(RES/"15_rung2_ensemble_benchmark.csv", index=False)
print(f"\n{'='*70}\nRINGKASAN — framing B (endpoint utama), R2 terbaik per set deskriptor")
for name in SETS:
    s = out[out.descriptors == name]
    b = s.loc[s.r2_B.idxmax()]
    print(f"  {name:8s} R2 = {b.r2_B:+.3f}  ({b.algorithm})")
