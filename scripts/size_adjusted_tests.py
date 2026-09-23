#!/usr/bin/env python
"""Cliff-vs-control comparisons, adjusted for molecular size difference.

Cliff partners differ by 4.4 heavy atoms on average, controls by 1.8, so any
size-correlated descriptor separates the two groups for the wrong reason.
Buried surface area looked highly significant that way (p < 0.0001) purely
because it tracks heavy-atom difference at rho = +0.70.

Matching harder does not fix this. Across all 1,911 similar pairs in the AChE
pool, potency difference grows with size difference (rho = +0.26): pairs
differing by 0-2 heavy atoms are iso-potent 48% of the time, pairs differing
by 6-9 only 18% of the time. Size-matched iso-potent controls are therefore
intrinsically rare, and forcing the match cost 26 of 86 control pairs while
still leaving a gap (2.9 vs 4.4).

So the confound is removed statistically instead, two ways:

  stratified   compare cliff vs control only within the same size bin
  residual     regress the descriptor on size difference and test the part
               size cannot explain

A descriptor that survives both is not a size artefact.
"""
from __future__ import annotations
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LinearRegression

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET = "CHEMBL220"

pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
lig = pd.concat([
    pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id","smiles_std"]],
    pd.read_csv(RES/f"04_validation_ligands_{TARGET}.csv")[["molecule_chembl_id","smiles_std"]],
]).drop_duplicates("molecule_chembl_id").set_index("molecule_chembl_id")
nh = {c: Chem.MolFromSmiles(s).GetNumHeavyAtoms()
      for c, s in lig.smiles_std.items() if Chem.MolFromSmiles(s)}

r2 = pd.read_csv(RES/"08_rung2_descriptors.csv").set_index("molecule_chembl_id")
br = pd.read_csv(RES/"08_rung2_bridging.csv").set_index("molecule_chembl_id")
r3 = pd.read_csv(RES/"09_rung3_mmgbsa.csv").drop_duplicates(
        "molecule_chembl_id").set_index("molecule_chembl_id")
D = r2.join(br, how="outer").join(r3, how="outer", rsuffix="_e")

TEST = ["bridge_gap", "d_cas", "d_pas", "catpi_cas", "catpi_pas",
        "axis_span_frac", "axis_reach_pas", "end_to_end", "n_contacts_4A",
        "frac_lig_buried", "dE_vdw", "dE_elec", "dG_polar_gb", "dG_nonpolar",
        "dG_total", "dSASA", "res_86", "res_286", "res_337"]
TEST = [c for c in TEST if c in D.columns]

rows = []
for _, r in pairs.iterrows():
    a, b = r.ligand_a, r.ligand_b
    if a not in D.index or b not in D.index or a not in nh or b not in nh:
        continue
    rec = {"group": r.group, "dp": r.delta, "d_heavy": abs(nh[a]-nh[b])}
    for c in TEST:
        va, vb = D.loc[a, c], D.loc[b, c]
        rec[c] = abs(va - vb) if pd.notna(va) and pd.notna(vb) else np.nan
    rows.append(rec)
x = pd.DataFrame(rows)
cl, ct = x[x.group=="cliff"], x[x.group=="control"]
print(f"pairs: {len(cl)} cliff, {len(ct)} control")
print(f"size difference: cliff {cl.d_heavy.mean():.1f}, control {ct.d_heavy.mean():.1f} "
      f"heavy atoms\n")

BINS = [0, 2, 4, 6, 100]
x["_sz"] = pd.cut(x.d_heavy, BINS, labels=False, include_lowest=True)

print(f"{'descriptor':18} {'raw p':>9} {'stratified p':>14} {'residual p':>12}  verdict")
print("-"*74)
out = []
for c in TEST:
    s = x[["group", c, "d_heavy", "_sz"]].dropna()
    a, b = s[s.group=="cliff"][c], s[s.group=="control"][c]
    if len(a) < 5 or len(b) < 5:
        continue
    p_raw = mannwhitneyu(a, b).pvalue

    # stratified: combine within-bin tests, weighted by bin size
    ps, ws = [], []
    for z in sorted(s._sz.dropna().unique()):
        g = s[s._sz == z]
        ga, gb = g[g.group=="cliff"][c], g[g.group=="control"][c]
        if len(ga) >= 4 and len(gb) >= 4:
            ps.append(mannwhitneyu(ga, gb).pvalue); ws.append(len(g))
    p_str = float(np.average(ps, weights=ws)) if ps else np.nan

    # residual: strip the part size explains, then test what is left
    m = LinearRegression().fit(s[["d_heavy"]], s[c])
    res = s[c] - m.predict(s[["d_heavy"]])
    p_res = mannwhitneyu(res[s.group=="cliff"], res[s.group=="control"]).pvalue

    ok = (p_str < 0.05) and (p_res < 0.05)
    verdict = "SURVIVES" if ok else ("size artefact" if p_raw < 0.05 else "-")
    print(f"{c:18} {p_raw:>9.4f} {p_str:>14.4f} {p_res:>12.4f}  {verdict}")
    out.append({"descriptor": c, "p_raw": p_raw, "p_stratified": p_str,
                "p_residual": p_res, "survives": ok,
                "cliff_mean": a.mean(), "control_mean": b.mean()})

pd.DataFrame(out).to_csv(RES/"10_size_adjusted_tests.csv", index=False)
n_ok = sum(o["survives"] for o in out)
n_art = sum(1 for o in out if o["p_raw"] < 0.05 and not o["survives"])
print(f"\n{n_ok} of {len(out)} descriptors survive size adjustment")
print(f"{n_art} were significant raw but are size artefacts")
print(f"\nwritten: 10_size_adjusted_tests.csv")
