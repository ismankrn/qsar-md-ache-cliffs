#!/usr/bin/env python
"""Does the 2D baseline change if it sees the protonated ligand?

Rung 1 is computed on the standardised neutral structures, while docking, MM/GBSA and
the simulations all use the dominant protonation state at pH 7.4. The ladder therefore
does not describe one chemical species from bottom to top, and the two forms are not a
cosmetic difference: 232 of the 271 ligands differ between the two SMILES columns, and
ECFP4 scores the neutral and protonated form of the same molecule at a median Tanimoto
of 0.706.

Which way that cuts cannot be guessed. A protonated baseline that scores higher means
the published rung 1 is understated and every structure-based rung looks better than it
should; one that scores lower means the opposite. Either way it bears on the headline
claim, so it is measured rather than argued.

Everything except the SMILES column is held identical to rung1_benchmark.py, i.e. the
same frozen folds, the same two feature sets, the same six and five algorithm families,
the same pipelines and the same three framings.

Writes results/33_rung1_protonation.csv.

Usage:
    python scripts/rung1_protonation.py
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
from sklearn.metrics import r2_score, roc_auc_score

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED, NF = "CHEMBL220", 20260902, 5


def pipe(est, scale=True):
    steps = [("var", VarianceThreshold(0.0))]
    if scale:
        steps.append(("sc", StandardScaler()))
    steps.append(("est", est))
    return Pipeline(steps)


REG = {
 "Ridge (linear)":      lambda: pipe(Ridge(alpha=10.0)),
 "PLS":                 lambda: pipe(PLSRegression(n_components=5)),
 "k-NN (k=5)":          lambda: pipe(KNeighborsRegressor(n_neighbors=5)),
 "SVM (RBF)":           lambda: pipe(SVR(C=10.0, gamma="scale")),
 "Random Forest":       lambda: pipe(RandomForestRegressor(
                            n_estimators=500, random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting":   lambda: pipe(HistGradientBoostingRegressor(
                            random_state=SEED), False),
}
CLF = {
 "Logistic regression": lambda: pipe(LogisticRegression(max_iter=3000, C=1.0)),
 "k-NN (k=5)":          lambda: pipe(KNeighborsClassifier(n_neighbors=5)),
 "SVM (RBF)":           lambda: pipe(SVC(C=10.0, probability=True, random_state=SEED)),
 "Random Forest":       lambda: pipe(RandomForestClassifier(
                            n_estimators=500, random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting":   lambda: pipe(HistGradientBoostingClassifier(
                            random_state=SEED), False),
}


def cv(X, y, f, make, clf=False):
    out = np.zeros(len(y))
    for k in range(NF):
        tr, te = f != k, f == k
        if te.sum() == 0 or tr.sum() == 0:
            continue
        m = make(); m.fit(X[tr], y[tr])
        out[te] = (m.predict_proba(X[te])[:, 1] if clf
                   else np.asarray(m.predict(X[te])).ravel())
    return out


def main() -> None:
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig = lig.merge(sl[["molecule_chembl_id", "split", "cv_fold"]],
                    on="molecule_chembl_id")
    dev = lig[lig.split == "dev"].reset_index(drop=True)
    dpair = pairs[pairs.split == "dev"].reset_index(drop=True)
    idx = {c: i for i, c in enumerate(dev.molecule_chembl_id)}
    fold, y_lig = dev.cv_fold.to_numpy(), dev.pchembl.to_numpy()

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(
        [d[0] for d in Descriptors._descList])

    def features(col):
        mols = [Chem.MolFromSmiles(s) for s in dev[col]]
        bit = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], float)
        cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], float)
        rdk = np.nan_to_num(np.array([calc.CalcDescriptors(m) for m in mols], float),
                            nan=0.0, posinf=0.0, neginf=0.0)
        return {"binary ECFP4": bit, "count ECFP4 + RDKit2D": np.hstack([cnt, rdk])}

    def pair_xy(X, mode):
        xs, ys, fs = [], [], []
        for _, r in dpair.iterrows():
            if r.ligand_a not in idx or r.ligand_b not in idx:
                continue
            da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
            pa, pb = y_lig[idx[r.ligand_a]], y_lig[idx[r.ligand_b]]
            if mode == "anti":
                xs += [da-db, db-da]; ys += [pa-pb, pb-pa]; fs += [r.cv_fold]*2
            else:
                xs.append(np.concatenate([da+db, np.abs(da-db)]))
                ys.append(1 if r.group == "cliff" else 0); fs.append(r.cv_fold)
        return np.array(xs), np.array(ys), np.array(fs)

    rows = []
    for form, col in (("neutral", "smiles_std"), ("protonated", "smiles_md")):
        FEAT = features(col)
        for fname, X in FEAT.items():
            for aname, make in REG.items():
                rows.append(dict(form=form, framing="A", features=fname,
                                 algorithm=aname,
                                 score=r2_score(y_lig, cv(X, y_lig, fold, make))))
            Xp, yp, fp = pair_xy(X, "anti")
            for aname, make in REG.items():
                rows.append(dict(form=form, framing="B", features=fname,
                                 algorithm=aname,
                                 score=r2_score(yp, cv(Xp, yp, fp, make))))
            Xc, yc, fc = pair_xy(X, "sym")
            for aname, make in CLF.items():
                rows.append(dict(form=form, framing="C", features=fname,
                                 algorithm=aname,
                                 score=roc_auc_score(yc, cv(Xc, yc, fc, make, clf=True))))
            print(f"  {form:11s} {fname:22s} done", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(RES/"33_rung1_protonation.csv", index=False)
    print("\n  -> results/33_rung1_protonation.csv\n")

    # rung 1 takes the best of the two feature sets, so the comparison does too
    best = d.loc[d.groupby(["form", "framing"]).score.idxmax()]
    print(f"  {'framing':9s}{'neutral':>10s}{'protonated':>12s}{'change':>9s}   winning arm")
    for fr in "ABC":
        n = best[(best.form == "neutral") & (best.framing == fr)].iloc[0]
        p = best[(best.form == "protonated") & (best.framing == fr)].iloc[0]
        print(f"  {fr:9s}{n.score:10.3f}{p.score:12.3f}{p.score-n.score:+9.3f}   "
              f"{n.features} / {n.algorithm}  ->  {p.features} / {p.algorithm}")


if __name__ == "__main__":
    main()
