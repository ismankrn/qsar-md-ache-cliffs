#!/usr/bin/env python
"""Rung 3 analysis — does energetics separate cliffs where geometry did not?

Three questions, in order of how much they matter:

  1. Do the MM/GBSA components differ more between cliff partners than
     between matched control partners?
  2. Does the per-residue decomposition support bridging? A ligand engaging
     both anionic sites should draw energy from Trp86 AND Trp286 at once,
     while its short-linker partner draws from only one. This tests the same
     hypothesis rung 2 tested geometrically, but through energy instead.
  3. Full model on the frozen folds, directly comparable to rungs 1 and 2.

Absolute dG values are not interpretable here -- single-snapshot MM/GBSA on
a charged ligand in an anionic gorge produced a total of about +8 kcal/mol,
which is unphysical. Differences within a matched pair are used throughout,
because the two partners share a protocol and a receptor, so systematic error
largely cancels.
"""
from __future__ import annotations
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import r2_score

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED, NF = "CHEMBL220", 20260902, 5

d = pd.read_csv(RES/"09_rung3_mmgbsa.csv").drop_duplicates(
        "molecule_chembl_id").set_index("molecule_chembl_id")
pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
print(f"ligands with MM/GBSA: {len(d)}\n")

COMP = ["dE_vdw", "dE_elec", "dG_polar_gb", "dG_nonpolar", "dG_mm_gb",
        "dG_total", "dSASA"]

# ---- Q1: component differences, cliff vs control ---------------------
rows = []
for _, r in pairs.iterrows():
    if r.ligand_a in d.index and r.ligand_b in d.index:
        rec = {"group": r.group, "split": r.split, "dp": r.delta}
        for c in COMP:
            rec["d_"+c] = abs(d.loc[r.ligand_a, c] - d.loc[r.ligand_b, c])
        # signed bridging contrast: does one partner draw from BOTH
        # tryptophans while the other draws from only one?
        for res in (86, 286, 337):
            rec[f"d_res_{res}"] = abs(d.loc[r.ligand_a, f"res_{res}"]
                                      - d.loc[r.ligand_b, f"res_{res}"])
        a86, b86 = d.loc[r.ligand_a,"res_86"], d.loc[r.ligand_b,"res_86"]
        a286, b286 = d.loc[r.ligand_a,"res_286"], d.loc[r.ligand_b,"res_286"]
        # "bridge energy" = how much the weaker of the two anchors contributes
        rec["d_bridge_energy"] = abs(max(a86, a286) - max(b86, b286))
        rows.append(rec)
x = pd.DataFrame(rows)
cl, ct = x[x.group=="cliff"], x[x.group=="control"]

print("Q1 — DIFFERENCE WITHIN A PAIR")
print(f"  {'component':18} {'cliff':>9} {'control':>9} {'ratio':>7} {'p':>9}")
print("  " + "-"*56)
for c in COMP + ["d_res_86", "d_res_286", "d_res_337", "d_bridge_energy"]:
    col = "d_"+c if c in COMP else c
    if col not in x.columns: continue
    a, b = cl[col].dropna(), ct[col].dropna()
    if len(a) < 3 or len(b) < 3: continue
    pv = mannwhitneyu(a, b).pvalue
    ratio = a.mean()/b.mean() if b.mean() else np.nan
    print(f"  {col.replace('d_',''):18} {a.mean():>9.2f} {b.mean():>9.2f} "
          f"{ratio:>7.2f} {pv:>9.4f}{'  *' if pv < 0.05 else ''}")
print(f"\n  n: {len(cl)} cliff, {len(ct)} control")

# ---- Q2: does any component RANK cliff magnitude? -------------------
print("\nQ2 — RANK CORRELATION WITH |dp|, cliff pairs only")
for c in COMP + ["d_res_86", "d_res_286", "d_bridge_energy"]:
    col = "d_"+c if c in COMP else c
    if col not in x.columns: continue
    rho = spearmanr(cl.dp, cl[col]).statistic
    print(f"  {col.replace('d_',''):18} rho = {rho:+.3f}")

# ---- Q3: model on the frozen folds ----------------------------------
sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id","pchembl"]]
lig = lig.merge(sl[["molecule_chembl_id","split","cv_fold"]], on="molecule_chembl_id")
dev = lig[(lig.split=="dev") & lig.molecule_chembl_id.isin(d.index)].reset_index(drop=True)
feat = [c for c in d.columns if c.startswith(("dE_","dG_","res_","dSASA"))]
X = d.loc[dev.molecule_chembl_id, feat].to_numpy(float)
y = dev.pchembl.to_numpy(); fold = dev.cv_fold.to_numpy()
idx = {c:i for i,c in enumerate(dev.molecule_chembl_id)}
dpair = pairs[pairs.split=="dev"]

def pipe(est, sc=True):
    st=[("imp",SimpleImputer(strategy="median")),("var",VarianceThreshold(0.0))]
    if sc: st.append(("s",StandardScaler()))
    st.append(("e",est)); return Pipeline(st)
REG = {"Ridge": lambda: pipe(Ridge(alpha=10.)),
       "PLS": lambda: pipe(PLSRegression(n_components=5)),
       "k-NN": lambda: pipe(KNeighborsRegressor(5)),
       "SVM": lambda: pipe(SVR(C=10.)),
       "RF": lambda: pipe(RandomForestRegressor(n_estimators=500,
                          random_state=SEED, n_jobs=-1), False),
       "GB": lambda: pipe(HistGradientBoostingRegressor(random_state=SEED), False)}

xs, ys, fs = [], [], []
for _, r in dpair.iterrows():
    if r.ligand_a in idx and r.ligand_b in idx:
        a, b = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        pa, pb = y[idx[r.ligand_a]], y[idx[r.ligand_b]]
        xs += [a-b, b-a]; ys += [pa-pb, pb-pa]; fs += [r.cv_fold]*2
Xp, yp, fp = np.array(xs), np.array(ys), np.array(fs)

print(f"\nQ3 — FRAMING B, frozen folds  ({len(dev)} ligands, {len(yp)} pair rows)")
best = -9
for name, mk in REG.items():
    pred = np.zeros(len(yp))
    for k in range(NF):
        tr, te = fp != k, fp == k
        if te.sum()==0: continue
        m = mk(); m.fit(Xp[tr], yp[tr])
        pred[te] = np.asarray(m.predict(Xp[te])).ravel()
    r2 = r2_score(yp, pred); best = max(best, r2)
    print(f"  {name:6} R2 = {r2:+.3f}   rho = {spearmanr(yp,pred).statistic:+.3f}")

print("\nLADDER COMPARISON — framing B, best R2")
r1 = pd.read_csv(RES/"05_rung1_benchmark.csv")
r2f = pd.read_csv(RES/"08_rung2_benchmark.csv")
print(f"  rung 1  (2D)            {r1[r1.framing=='B'].r2.max():+.3f}")
print(f"  rung 2  (pose geometry) {r2f[r2f.framing=='B'].r2.max():+.3f}")
print(f"  rung 3  (MM/GBSA)       {best:+.3f}")
x.to_csv(RES/"09_rung3_pair_analysis.csv", index=False)
print(f"\nwritten: 09_rung3_pair_analysis.csv")
