#!/usr/bin/env python
"""Build the continuous-delta external validation set.

The curated study set is case-control: it uses only |dp| <= 0.30 and >= 2.00,
discarding the 56% of similar pairs that fall in between. A model scored only
on that set cannot be claimed to predict potency differences - only to
separate two hand-picked extremes.

This set samples pairs across the FULL delta range with no gap, so the model
can be tested on the continuum that actually exists in nature.

Leakage rule: no ligand here may appear in the DEV split of the study set.
Reusing dev ligands would let the model exploit activities it trained on.
Ligands from the sealed TEST split are permitted (never trained on) and are
free - their MD is already budgeted.

Usage: python build_validation_set.py [--per-bin 18]
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import BulkTanimotoSimilarity

from build_ligand_set import md_ok, embed_ok, protonate, SIM_CUTOFF
from protonation import protonate as _p

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
TARGET = "CHEMBL220"
SEED = 20260902

# contiguous bins covering the whole range - deliberately NO gap
BINS = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 2.0, 2.5, 3.0, 6.0]

ap = argparse.ArgumentParser()
ap.add_argument("--per-bin", type=int, default=18)
args = ap.parse_args()
rng = np.random.default_rng(SEED)

pool = pd.read_csv(RESULTS / f"01_compounds_{TARGET}.csv")
study = pd.read_csv(RESULTS / f"02_ligand_set_{TARGET}.csv")
splits = pd.read_csv(RESULTS / f"03_splits_ligands_{TARGET}.csv")
old_pairs = pd.read_csv(RESULTS / f"02_pairs_{TARGET}.csv")

dev_ids = set(splits[splits.split == "dev"].molecule_chembl_id)
test_ids = set(splits[splits.split == "test"].molecule_chembl_id)
print(f"study set: {len(dev_ids)} dev ligands (forbidden), "
      f"{len(test_ids)} test ligands (allowed, already simulated)")

# ---- MD-feasible pool
chk = pool["smiles_std"].map(md_ok)
pool = pool[[c[0] for c in chk]].reset_index(drop=True)
mols = [Chem.MolFromSmiles(s) for s in pool.smiles_std]
keep = [i for i, m in enumerate(mols) if m is not None]
pool, mols = pool.iloc[keep].reset_index(drop=True), [mols[i] for i in keep]
print(f"MD-feasible pool: {len(pool):,}")

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [gen.GetFingerprint(m) for m in mols]
flat = [Chem.MolToSmiles(m, isomericSmiles=False) for m in mols]
pot = pool.pchembl.to_numpy()
cid = pool.molecule_chembl_id.tolist()

seen = {frozenset((r.ligand_a, r.ligand_b)) for _, r in old_pairs.iterrows()}

print("enumerating candidate pairs ...")
cand = []
for i in range(len(fps) - 1):
    if cid[i] in dev_ids:
        continue
    sims = np.array(BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
    hit = np.where(sims >= SIM_CUTOFF)[0]
    for h in hit:
        j = i + 1 + h
        if cid[j] in dev_ids:                 # leakage guard
            continue
        if flat[i] == flat[j]:                # stereochemical pair
            continue
        if frozenset((cid[i], cid[j])) in seen:
            continue
        cand.append({"ligand_a": cid[i], "ligand_b": cid[j],
                     "smiles_a": pool.smiles_std[i], "smiles_b": pool.smiles_std[j],
                     "pchembl_a": pot[i], "pchembl_b": pot[j],
                     "delta": float(abs(pot[j] - pot[i])),
                     "tanimoto": float(sims[h])})
cand = pd.DataFrame(cand)
print(f"candidate pairs (no dev ligand, not already used): {len(cand):,}")

# prefer pairs that reuse already-simulated test ligands: those cost nothing
cand["n_new"] = [(a not in test_ids) + (b not in test_ids)
                 for a, b in zip(cand.ligand_a, cand.ligand_b)]
cand["bin"] = pd.cut(cand.delta, BINS, labels=False, include_lowest=True)

print(f"\n{'bin':>14} {'tersedia':>9} {'diambil':>8}")
picked, used = [], set()
for b in range(len(BINS) - 1):
    block = cand[cand.bin == b].sort_values(
        ["n_new"], kind="stable")                      # cheapest first
    block = block.sample(frac=1.0, random_state=SEED).sort_values("n_new", kind="stable")
    take = []
    for _, r in block.iterrows():
        if len(take) >= args.per_bin:
            break
        take.append(r)
        used.update([r.ligand_a, r.ligand_b])
    picked.extend(take)
    print(f"  {BINS[b]:.1f}-{BINS[b+1]:.1f} {len(block):>9,} {len(take):>8}")

val = pd.DataFrame(picked).reset_index(drop=True)
val["pair_id"] = [f"val_{i:03d}" for i in range(len(val))]

# ---- ligand table for the new set
rows = {}
for _, r in val.iterrows():
    for s in ("a", "b"):
        c = r[f"ligand_{s}"]
        rows.setdefault(c, {"molecule_chembl_id": c,
                            "smiles_std": r[f"smiles_{s}"],
                            "pchembl": r[f"pchembl_{s}"]})
vlig = pd.DataFrame(rows.values()).reset_index(drop=True)
vlig["already_simulated"] = vlig.molecule_chembl_id.isin(test_ids)

vlig["embed_ok"] = vlig.smiles_std.map(embed_ok)
dead = set(vlig[~vlig.embed_ok].molecule_chembl_id)
if dead:
    val = val[~val.ligand_a.isin(dead) & ~val.ligand_b.isin(dead)]
    vlig = vlig[vlig.embed_ok]
vlig = vlig.drop(columns=["embed_ok"]).reset_index(drop=True)

prot = vlig.smiles_std.map(protonate)
vlig["smiles_md"] = [p[0] for p in prot]
vlig["formal_charge"] = [p[1] for p in prot]
vlig["prot_ambiguous"] = [bool(p[2]) for p in prot]

used_ids = set(val.ligand_a) | set(val.ligand_b)
vlig = vlig[vlig.molecule_chembl_id.isin(used_ids)].reset_index(drop=True)
n_new = int((~vlig.already_simulated).sum())

val.drop(columns=["bin"]).to_csv(RESULTS / f"04_validation_pairs_{TARGET}.csv", index=False)
vlig.to_csv(RESULTS / f"04_validation_ligands_{TARGET}.csv", index=False)

print("\n" + "=" * 62)
print(f"PASANGAN VALIDASI : {len(val)}")
print(f"LIGAN TERLIBAT    : {len(vlig)}"
      f"  ({len(vlig)-n_new} sudah disimulasikan, {n_new} BARU)")
print(f"rentang delta p   : {val.delta.min():.2f} - {val.delta.max():.2f}  (tanpa celah)")
print(f"ligan dev yang bocor: {len(used_ids & dev_ids)}  (harus 0)")
total = 271 + n_new
ns = total * 2 * 25 * 1.15 + 300
print(f"\ntotal ligan MD    : 271 + {n_new} = {total}")
print(f"biaya MD          : ${ns/350*24*0.34:,.0f}"
      f"  (naik dari $370, tambahan ${ns/350*24*0.34-370:,.0f})")
