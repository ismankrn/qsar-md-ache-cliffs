#!/usr/bin/env python
"""Freeze the train/test split ONCE, for every rung of the comparison.

Why this is a separate script producing a file, rather than a line inside
each model script:

  1. All four rungs (2D / static-3D / MM-GBSA / dynamic) must be scored on
     IDENTICAL folds. Recomputing folds per script means the comparison
     measures fold luck as much as descriptor quality.
  2. The split must be fixed BEFORE anyone sees MD results, so it cannot be
     tuned - knowingly or not - to flatter the conclusion.
  3. Ligands share pairs, so folds are built from connected components of
     the ligand-pair graph. Splitting by ligand would put one partner of a
     cliff in train and the other in test, letting the model memorise the
     scaffold and inflating every score.

Produces an external test set that stays SEALED until the final model is
chosen, plus 5 CV folds over the development set. With four rungs and
several algorithms each, many comparisons happen; a sealed test set is the
only defence against selecting a winner on noise.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
TARGET = "CHEMBL220"
SEED = 20260902
N_FOLDS = 5
TEST_FRACTION = 0.20

rng = np.random.default_rng(SEED)
lig = pd.read_csv(RESULTS / f"02_ligand_set_{TARGET}.csv")
pairs = pd.read_csv(RESULTS / f"02_pairs_{TARGET}.csv")

# ---- connected components over the ligand-pair graph ---------------
idx = {c: i for i, c in enumerate(lig.molecule_chembl_id)}
parent = list(range(len(lig)))
def find(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]; a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb
for _, r in pairs.iterrows():
    if r.ligand_a in idx and r.ligand_b in idx:
        union(idx[r.ligand_a], idx[r.ligand_b])
lig["component"] = [find(i) for i in range(len(lig))]

comp = lig.groupby("component").agg(
    n_ligands=("molecule_chembl_id", "size"),
    mean_p=("pchembl", "mean"),
).reset_index()

# does each component contain a cliff pair?
pair_comp = {}
for _, r in pairs.iterrows():
    if r.ligand_a in idx:
        pair_comp.setdefault(find(idx[r.ligand_a]), set()).add(r.group)
comp["has_cliff"] = comp.component.map(lambda c: "cliff" in pair_comp.get(c, set()))

# cliff pairs per component - the quantity the folds must balance
cliff_count = {}
for _, r in pairs.iterrows():
    if r.group == "cliff" and r.ligand_a in idx:
        c = find(idx[r.ligand_a])
        cliff_count[c] = cliff_count.get(c, 0) + 1
comp["n_cliff_pairs"] = comp.component.map(lambda c: cliff_count.get(c, 0))
print(f"ligands {len(lig)} | components {len(comp)} "
      f"({comp.has_cliff.sum()} contain cliff pairs)")

# ---- stratified assignment: balance cliff-containing components and
# ---- potency across the test set and every CV fold
comp["p_stratum"] = pd.qcut(comp.mean_p, 3, labels=False, duplicates="drop")
comp["stratum"] = comp.has_cliff.astype(str) + "_" + comp.p_stratum.astype(str)

comp["split"] = "dev"
comp["cv_fold"] = -1
for _, block in comp.groupby("stratum"):
    ids = block.component.to_numpy()
    rng.shuffle(ids)
    n_test = max(1, int(round(TEST_FRACTION * len(ids))))
    comp.loc[comp.component.isin(ids[:n_test]), "split"] = "test"

# Greedy bin-packing over the dev components, balancing CLIFF PAIRS first.
# Round-robin over components fails badly here: components differ hugely in
# size, so some folds ended up with 25 cliff pairs and others with 4, making
# per-fold cliff metrics uselessly noisy.
dev = comp[comp.split == "dev"].sort_values(
    ["n_cliff_pairs", "n_ligands"], ascending=False)
# Balance cliff pairs AND ligand count jointly: optimising cliffs alone
# gave folds of 13 vs 63 ligands, because dense cliff clusters are small
# while the large components are cliff-free.
tgt_cliff = max(dev.n_cliff_pairs.sum() / N_FOLDS, 1)
tgt_lig = max(dev.n_ligands.sum() / N_FOLDS, 1)
load = {k: {"cliff": 0, "lig": 0} for k in range(N_FOLDS)}
for _, r in dev.iterrows():
    k = min(load, key=lambda f: (load[f]["cliff"] / tgt_cliff
                                 + load[f]["lig"] / tgt_lig))
    comp.loc[comp.component == r.component, "cv_fold"] = k
    load[k]["cliff"] += r.n_cliff_pairs
    load[k]["lig"] += r.n_ligands

lig = lig.merge(comp[["component", "split", "cv_fold", "has_cliff"]],
                on="component", how="left")

# ---- pair-level split (a pair inherits its component's split) ------
lig_split = dict(zip(lig.molecule_chembl_id, lig.split))
lig_fold = dict(zip(lig.molecule_chembl_id, lig.cv_fold))
pairs["split"] = pairs.ligand_a.map(lig_split)
pairs["cv_fold"] = pairs.ligand_a.map(lig_fold)
bad = sum(lig_split.get(a) != lig_split.get(b)
          for a, b in zip(pairs.ligand_a, pairs.ligand_b))
print(f"pairs split across train/test boundary: {bad}  (must be 0)")

# ---- report --------------------------------------------------------
print(f"\n{'':10} {'ligands':>8} {'cliff pairs':>12} {'control pairs':>14}")
for s in ("dev", "test"):
    ls = (lig.split == s).sum()
    cp = ((pairs.split == s) & (pairs.group == "cliff")).sum()
    op = ((pairs.split == s) & (pairs.group == "control")).sum()
    print(f"{s:10} {ls:>8} {cp:>12} {op:>14}")

print(f"\nCV folds within dev set:")
for k in range(N_FOLDS):
    ls = ((lig.split == "dev") & (lig.cv_fold == k)).sum()
    cp = ((pairs.split == "dev") & (pairs.cv_fold == k) & (pairs.group == "cliff")).sum()
    print(f"  fold {k}: {ls:>3} ligands, {cp:>2} cliff pairs")

out_l = RESULTS / f"03_splits_ligands_{TARGET}.csv"
out_p = RESULTS / f"03_splits_pairs_{TARGET}.csv"
lig[["molecule_chembl_id", "component", "split", "cv_fold", "group"]].to_csv(out_l, index=False)
pairs[["pair_id", "group", "ligand_a", "ligand_b", "delta", "split", "cv_fold"]].to_csv(out_p, index=False)
print(f"\nFROZEN. seed={SEED}. Every rung must load these files.")
print(f"written: {out_l.name}, {out_p.name}")
