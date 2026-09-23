#!/usr/bin/env python
"""Stage 07 analysis — what does the raw docking score actually explain?

Two separate questions, and they have different answers in the literature:

  1. Does the score rank potency across the whole set?
  2. Does the DIFFERENCE in score predict the potency difference within a
     matched pair?

The second is what matters here. It is the crudest possible structure-based
descriptor, so it sets the floor that rung 2 and rung 3 must beat.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
TARGET = "CHEMBL220"

d = pd.read_csv(RES / "07_docking_scores.csv")
pairs = pd.read_csv(RES / f"03_splits_pairs_{TARGET}.csv")

print(f"docking runs : {len(d):,}")
print(f"ligands      : {d.molecule_chembl_id.nunique()}")
print(f"failures     : {d.score_1.isna().sum()}\n")

# ---- seed reproducibility -------------------------------------------
g = d.groupby("molecule_chembl_id").score_1
rep = pd.DataFrame({"score": g.median(), "spread": g.max() - g.min(),
                    "sd": g.std(), "n": g.size()})
print("SEED REPRODUCIBILITY")
print(f"  spread across 3 seeds : median {rep.spread.median():.2f}, "
      f"p95 {rep.spread.quantile(.95):.2f}, max {rep.spread.max():.2f} kcal/mol")
unstable = rep[rep.spread > 1.0]
print(f"  ligands with spread > 1.0 kcal/mol : {len(unstable)} "
      f"({100*len(unstable)/len(rep):.1f}%)")
print("  -> the median pose is used below; unstable ligands are flagged\n")

# ---- Q1: does the score rank potency at all? ------------------------
lig = pd.read_csv(RES / f"02_ligand_set_{TARGET}.csv")[
    ["molecule_chembl_id", "pchembl"]]
val = pd.read_csv(RES / f"04_validation_ligands_{TARGET}.csv")[
    ["molecule_chembl_id", "pchembl"]]
pot = pd.concat([lig, val]).drop_duplicates("molecule_chembl_id")
m = rep.reset_index().merge(pot, on="molecule_chembl_id")

print("Q1 — SCORE vs POTENCY, across all ligands")
print(f"  n = {len(m)}")
print(f"  Pearson  r  = {pearsonr(m.score, m.pchembl)[0]:+.3f}")
print(f"  Spearman rho= {spearmanr(m.score, m.pchembl).statistic:+.3f}")
print("  (a useful score would be NEGATIVE: lower energy, higher potency)\n")

# ---- Q2: does the score difference predict the potency difference? ---
sm = dict(zip(rep.index, rep.score))
rows = []
for _, r in pairs.iterrows():
    if r.ligand_a in sm and r.ligand_b in sm:
        rows.append({"pair_id": r.pair_id, "group": r.group, "split": r.split,
                     "d_potency": r.delta,
                     "d_score": abs(sm[r.ligand_a] - sm[r.ligand_b])})
p = pd.DataFrame(rows)
print("Q2 — SCORE DIFFERENCE vs POTENCY DIFFERENCE, within matched pairs")
print(f"{'group':>9} {'n':>4} {'true |dp|':>10} {'|d score|':>10} {'rho':>8}")
print("-" * 46)
for grp in ("cliff", "control"):
    s = p[p.group == grp]
    rho = spearmanr(s.d_potency, s.d_score).statistic if len(s) > 2 else np.nan
    print(f"{grp:>9} {len(s):>4} {s.d_potency.mean():>10.2f} "
          f"{s.d_score.mean():>10.2f} {rho:>8.3f}")
cl, ct = p[p.group == "cliff"], p[p.group == "control"]
print(f"\n  cliff vs control separation in |d score|: "
      f"{cl.d_score.mean() - ct.d_score.mean():+.2f} kcal/mol")
rho_all = spearmanr(p.d_potency, p.d_score).statistic
print(f"  rho over all pairs: {rho_all:+.3f}")

p.to_csv(RES / "07_docking_pair_analysis.csv", index=False)
rep.reset_index().to_csv(RES / "07_docking_per_ligand.csv", index=False)
print(f"\nwritten: 07_docking_pair_analysis.csv, 07_docking_per_ligand.csv")
