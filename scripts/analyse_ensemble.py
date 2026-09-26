"""Stage 11 analysis — does ensemble docking rescue what rigid 4EY7 could not?

Rungs 1-3 all failed on the primary endpoint (framing B): 2D fingerprints R2 = +0.200,
static pose geometry +0.167, single-snapshot MM/GBSA -0.055, and the raw 4EY7 docking score
separated cliff from control pairs by -0.03 kcal/mol, i.e. not at all. Three explanations
remained: the poses are wrong, dynamics matter, or the phenomenon is not structural.

This stage tests the first. If cross-docking failures were an artefact of holding one
receptor rigid, then scoring against eight gorge conformations should recover signal.
If it does not, MD becomes the only untested explanation left.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"

ens = pd.read_csv(RES/"11_ensemble_docking.csv")
lig = pd.read_csv(RES/"02_ligand_set_CHEMBL220.csv")
prs = pd.read_csv(RES/"03_splits_pairs_CHEMBL220.csv")

print(f"data: {len(ens)} baris, {ens.molecule_chembl_id.nunique()} ligan, "
      f"{ens.receptor.nunique()} reseptor")
print(f"reseptor: {sorted(ens.receptor.unique())}\n")

# ---- per-ligand scores: single receptor vs ensemble best
best = (ens.groupby("molecule_chembl_id").score_1.min().rename("ens_best"))
mean_ = (ens.groupby("molecule_chembl_id").score_1.mean().rename("ens_mean"))
single = (ens[ens.receptor == "4ey7"].set_index("molecule_chembl_id").score_1.rename("s_4ey7"))
per = pd.concat([single, best, mean_], axis=1).reset_index()
per = per.merge(lig[["molecule_chembl_id", "pchembl"]], on="molecule_chembl_id", how="inner")
print(f"ligan dengan pChEMBL: {len(per)}")
print(f"perbaikan ensemble atas 4EY7: median {(per.s_4ey7-per.ens_best).median():.2f}, "
      f"maks {(per.s_4ey7-per.ens_best).max():.2f} kcal/mol\n")

print("=== 1. korelasi skor dengan potensi, per ligan ===")
for col, name in (("s_4ey7", "4EY7 saja"), ("ens_best", "ensemble terbaik"), ("ens_mean", "ensemble rata-rata")):
    r, p = stats.pearsonr(per[col], per.pchembl)
    rho, _ = stats.spearmanr(per[col], per.pchembl)
    print(f"  {name:20s} r = {r:+.3f}  (p = {p:.2g})   rho = {rho:+.3f}")

print("\n=== 2. apakah selisih skor memisahkan cliff dari kontrol? ===")
sc = {c: per.set_index("molecule_chembl_id")[c].to_dict()
      for c in ("s_4ey7", "ens_best", "ens_mean")}
rows = []
for _, r in prs.iterrows():
    e = {"pair_id": r.pair_id, "group": r.group, "split": r.split, "delta": r.delta}
    ok = True
    for c in sc:
        a, b = sc[c].get(r.ligand_a), sc[c].get(r.ligand_b)
        if a is None or b is None:
            ok = False; break
        e[f"d_{c}"] = a - b          # signed, antisymmetric (framing B)
        e[f"absd_{c}"] = abs(a - b)  # symmetric (framing C)
    if ok:
        rows.append(e)
pa = pd.DataFrame(rows)
print(f"  pasangan terpakai: {len(pa)} (cliff {int((pa.group=='cliff').sum())}, "
      f"kontrol {int((pa.group=='control').sum())})")
for c in ("s_4ey7", "ens_best", "ens_mean"):
    cl = pa.loc[pa.group == "cliff", f"absd_{c}"]
    ko = pa.loc[pa.group == "control", f"absd_{c}"]
    u, pv = stats.mannwhitneyu(cl, ko)
    d = (cl.mean() - ko.mean()) / np.sqrt((cl.var() + ko.var())/2)
    print(f"  {c:10s} |Δskor| cliff {cl.mean():.3f} vs kontrol {ko.mean():.3f} "
          f"-> selisih {cl.mean()-ko.mean():+.3f} kcal/mol, Cohen d {d:+.2f}, p {pv:.3f}")

print("\n=== 3. framing B: apakah Δskor memprediksi ΔpChEMBL bertanda? ===")
dev = pa[pa.split == "dev"]
print(f"  hanya split dev: {len(dev)} pasangan (test tetap tersegel)")
for c in ("s_4ey7", "ens_best", "ens_mean"):
    x, y = dev[f"d_{c}"], dev.delta
    r, p = stats.pearsonr(x, y)
    # R2 dari regresi linier satu prediktor, sebanding dengan angka rung
    r2 = r**2
    print(f"  {c:10s} r = {r:+.3f}  R2 = {r2:.3f}  (p = {p:.2g})")

print("\n=== 4. apakah ada satu reseptor yang lebih baik dari 4EY7? ===")
for rec in sorted(ens.receptor.unique()):
    s = ens[ens.receptor == rec].set_index("molecule_chembl_id").score_1
    m = per.molecule_chembl_id.map(s)
    ok = m.notna()
    r, _ = stats.pearsonr(m[ok], per.pchembl[ok])
    print(f"  {rec:6s} r dengan pChEMBL = {r:+.3f}   n = {int(ok.sum())}")

pa.to_csv(RES/"11_ensemble_pair_analysis.csv", index=False)
per.to_csv(RES/"11_ensemble_per_ligand.csv", index=False)
print(f"\ntersimpan: 11_ensemble_pair_analysis.csv, 11_ensemble_per_ligand.csv")
