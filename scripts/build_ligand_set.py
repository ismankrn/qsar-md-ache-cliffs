#!/usr/bin/env python
"""Curate the AChE simulation set: constitutional cliff pairs + matched controls.

The control set is what makes the study interpretable. If dynamic descriptors
separate cliff pairs from structurally-similar NON-cliff pairs, that is a
result. Without matched controls, any signal could just be "similar molecules
behave similarly".

Controls are matched to cliffs on both similarity and potency, so the two
groups differ in the property under test (activity difference) and as little
else as possible.

Usage:
    python build_ligand_set.py --target CHEMBL220 [--check-params]
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, AllChem
from rdkit.DataStructs import BulkTanimotoSimilarity

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

# --- MD feasibility limits -----------------------------------------
MW_MIN, MW_MAX = 150, 700
MAX_HEAVY = 55
MAX_ROTB = 18               # beyond this, 25 ns won't sample conformers
ALLOWED = set("H C N O F S Cl Br I P".split())

# --- control matching ----------------------------------------------
SIM_CUTOFF = 0.90           # same similarity floor as the cliff definition
CONTROL_MAX_DELTA = 0.30    # "same potency": under 2-fold
PH = 7.4

_fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def flat(smiles: str) -> str | None:
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m, isomericSmiles=False) if m else None


# ------------------------------------------------------- MD feasibility
def md_ok(smiles: str) -> tuple[bool, str]:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False, "unparseable"
    els = {a.GetSymbol() for a in m.GetAtoms()}
    if not els <= ALLOWED:
        return False, f"element {sorted(els - ALLOWED)}"
    if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
        return False, "radical"
    mw = Descriptors.MolWt(m)
    if not (MW_MIN <= mw <= MW_MAX):
        return False, f"MW {mw:.0f}"
    if m.GetNumHeavyAtoms() > MAX_HEAVY:
        return False, f"{m.GetNumHeavyAtoms()} heavy atoms"
    if Descriptors.NumRotatableBonds(m) > MAX_ROTB:
        return False, f"{Descriptors.NumRotatableBonds(m)} rot bonds"
    return True, "ok"


def embed_ok(smiles: str) -> bool:
    """A ligand that will not embed in 3D cannot be simulated."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    m = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = 0xC0FFEE
    return AllChem.EmbedMolecule(m, p) == 0


# --------------------------------------------------------- protonation
from protonation import protonate as _protonate_rules


def protonate(smiles: str) -> tuple[str, int, int]:
    """Dominant protonation state at pH 7.4 via validated rules.

    AChE's catalytic anionic site binds cationic ligands through cation-pi
    with Trp86, so simulating a neutralised amine would be plainly wrong.
    Third element is 1 when the assignment is ambiguous (e.g. imidazole).
    """
    r = _protonate_rules(smiles, ph=PH)
    return r["smiles"], r["charge"], int(r["ambiguous"])


# ------------------------------------------------------------ controls
def find_controls(pool: pd.DataFrame, cliffs: pd.DataFrame,
                  exclude: set[str], rng: np.random.Generator) -> pd.DataFrame:
    """Structurally similar pairs with near-identical potency, matched to the
    cliff pairs' joint distribution of similarity and potency."""
    sub = pool[~pool["molecule_chembl_id"].isin(exclude)].reset_index(drop=True)
    mols = [Chem.MolFromSmiles(s) for s in sub["smiles_std"]]
    keep = [i for i, m in enumerate(mols) if m is not None]
    sub = sub.iloc[keep].reset_index(drop=True)
    fps = [_fpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in sub["smiles_std"]]
    pot = sub["pchembl"].to_numpy()

    cand = []
    for i in range(len(fps) - 1):
        sims = np.array(BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
        d = np.abs(pot[i + 1:] - pot[i])
        hit = np.where((sims >= SIM_CUTOFF) & (d <= CONTROL_MAX_DELTA))[0]
        for h in hit:
            j = i + 1 + h
            cand.append({
                "cid_a": sub.molecule_chembl_id[i], "cid_b": sub.molecule_chembl_id[j],
                "smiles_a": sub.smiles_std[i], "smiles_b": sub.smiles_std[j],
                "pchembl_a": pot[i], "pchembl_b": pot[j],
                "delta": float(d[h]), "tanimoto": float(sims[h]),
            })
    cand = pd.DataFrame(cand)
    if cand.empty:
        return cand

    # drop stereochemical pairs: same trap as before
    cand = cand[[flat(a) != flat(b)
                 for a, b in zip(cand.smiles_a, cand.smiles_b)]].reset_index(drop=True)

    # match the cliff set's (similarity, potency) distribution
    cand["mean_p"] = (cand.pchembl_a + cand.pchembl_b) / 2
    cliffs = cliffs.copy()
    cliffs["mean_p"] = (cliffs.pchembl_a + cliffs.pchembl_b) / 2
    sim_bins = [0.90, 0.93, 0.96, 0.99, 1.01]
    pot_bins = [0, 5, 6, 7, 8, 20]
    cand["_sb"] = pd.cut(cand.tanimoto, sim_bins, labels=False)
    cand["_pb"] = pd.cut(cand.mean_p, pot_bins, labels=False)
    cliffs["_sb"] = pd.cut(cliffs.tanimoto, sim_bins, labels=False)
    cliffs["_pb"] = pd.cut(cliffs.mean_p, pot_bins, labels=False)

    want = cliffs.groupby(["_sb", "_pb"], observed=True).size()
    picked, used = [], set()
    for (sb, pb), n in want.items():
        block = cand[(cand._sb == sb) & (cand._pb == pb)]
        block = block[~block.cid_a.isin(used) & ~block.cid_b.isin(used)]
        if block.empty:
            continue
        take = block.sample(min(n, len(block)), random_state=int(rng.integers(1e6)))
        for _, r in take.iterrows():
            if r.cid_a in used or r.cid_b in used:
                continue
            picked.append(r)
            used.update([r.cid_a, r.cid_b])
    return pd.DataFrame(picked).drop(columns=["_sb", "_pb"], errors="ignore")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="CHEMBL220")
    ap.add_argument("--check-params", action="store_true",
                    help="test OpenFF parameterisation (slow, ~2s/ligand)")
    ap.add_argument("--seed", type=int, default=20260902)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    cliffs_all = pd.read_csv(RESULTS / f"01_cliffs_{a.target}.csv")
    cliffs_all = cliffs_all[cliffs_all.cliff_class == "constitutional"].reset_index(drop=True)
    pool = pd.read_csv(RESULTS / f"01_compounds_{a.target}.csv")
    print(f"constitutional cliff pairs: {len(cliffs_all)}")
    print(f"compound pool: {len(pool):,}")

    # ---- filter the WHOLE pool for MD feasibility FIRST.
    # Filtering ligands after pairing silently destroys pairs (a pair dies if
    # either partner fails); filtering the pool up front keeps every surviving
    # pair intact by construction.
    print("\nfiltering pool for MD feasibility ...")
    checks = pool["smiles_std"].map(md_ok)
    pool["md_ok"] = [c[0] for c in checks]
    pool["md_reason"] = [c[1] for c in checks]
    rej = pool[~pool.md_ok]
    print(f"  rejected {len(rej):,}/{len(pool):,}; top reasons:")
    print("   ", rej.md_reason.str.replace(r"[0-9]+", "N", regex=True)
          .value_counts().head(4).to_dict())
    pool = pool[pool.md_ok].drop(columns=["md_ok", "md_reason"]).reset_index(drop=True)
    ok_ids = set(pool.molecule_chembl_id)
    print(f"  pool after filtering: {len(pool):,}")

    # keep only cliff pairs whose BOTH partners survived
    cliffs = cliffs_all[cliffs_all.cid_a.isin(ok_ids)
                        & cliffs_all.cid_b.isin(ok_ids)].reset_index(drop=True)
    print(f"\ncliff pairs surviving intact: {len(cliffs)}/{len(cliffs_all)}")
    cliff_ids = set(cliffs.cid_a) | set(cliffs.cid_b)
    print(f"cliff compounds: {len(cliff_ids)}")

    print("\nbuilding matched control pairs ...")
    controls = find_controls(pool, cliffs, cliff_ids, rng)
    n_ctrl_cmpd = len(set(controls.cid_a) | set(controls.cid_b)) if len(controls) else 0
    print(f"control pairs: {len(controls)}  ({n_ctrl_cmpd} compounds)")

    # ---- assemble. A ligand is simulated ONCE; pairs REFERENCE ligands.
    # Collapsing the two into one table (and de-duplicating) destroys pair
    # membership, because a compound can belong to several pairs.
    pair_rows, lig_rows = [], {}
    for grp, tab in (("cliff", cliffs), ("control", controls)):
        if tab.empty:
            continue
        for pid, (_, r) in enumerate(tab.iterrows()):
            pair_rows.append({
                "pair_id": f"{grp}_{pid:03d}", "group": grp,
                "ligand_a": r["cid_a"], "ligand_b": r["cid_b"],
                "pchembl_a": r["pchembl_a"], "pchembl_b": r["pchembl_b"],
                "delta": r["delta"], "tanimoto": r["tanimoto"],
            })
            for side in ("a", "b"):
                cid = r[f"cid_{side}"]
                lig_rows.setdefault(cid, {
                    "molecule_chembl_id": cid,
                    "smiles_std": r[f"smiles_{side}"],
                    "pchembl": r[f"pchembl_{side}"],
                    "group": grp,
                })
    pairs = pd.DataFrame(pair_rows)
    lig = pd.DataFrame(lig_rows.values()).reset_index(drop=True)
    print(f"\nunique ligands to simulate: {len(lig)}")
    print(f"pairs to analyse: {len(pairs)}")

    # ---- 3D embedding
    lig["embed_ok"] = lig["smiles_std"].map(embed_ok)
    n_bad = int((~lig.embed_ok).sum())
    if n_bad:
        dead_ids = set(lig[~lig.embed_ok].molecule_chembl_id)
        lig = lig[lig.embed_ok]
        before = len(pairs)
        pairs = pairs[~pairs.ligand_a.isin(dead_ids) & ~pairs.ligand_b.isin(dead_ids)]
        print(f"3D embedding failures: {n_bad} ligands, breaking "
              f"{before - len(pairs)} pairs")
    lig = lig.drop(columns=["embed_ok"]).reset_index(drop=True)

    # ---- protonation at pH 7.4
    print(f"\nassigning protonation states at pH {PH} ...")
    prot = lig["smiles_std"].map(protonate)
    lig["smiles_md"] = [p[0] for p in prot]
    lig["formal_charge"] = [p[1] for p in prot]
    lig["prot_ambiguous"] = [bool(p[2]) for p in prot]
    print(lig.formal_charge.value_counts().sort_index().to_string())
    print(f"ambiguous protonation: {int(lig.prot_ambiguous.sum())}")

    # ---- optional OpenFF check
    if a.check_params:
        print("\ntesting OpenFF parameterisation ...")
        from openff.toolkit import Molecule
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        from openmm import app
        ok = []
        for i, s in enumerate(lig.smiles_md, 1):
            try:
                m = Molecule.from_smiles(s, allow_undefined_stereo=True)
                m.generate_conformers(n_conformers=1)
                gen = SMIRNOFFTemplateGenerator(molecules=m, forcefield="openff-2.3.0")
                ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
                ff.registerTemplateGenerator(gen.generator)
                ff.createSystem(m.to_topology().to_openmm())
                ok.append(True)
            except Exception:
                ok.append(False)
            if i % 25 == 0:
                print(f"   {i}/{len(lig)}", flush=True)
        lig["param_ok"] = ok
        print(f"parameterisation failures: {int((~lig.param_ok).sum())}")
        lig = lig[lig.param_ok].drop(columns=["param_ok"]).reset_index(drop=True)

    # drop any ligand no longer referenced by a surviving pair
    used = set(pairs.ligand_a) | set(pairs.ligand_b)
    lig = lig[lig.molecule_chembl_id.isin(used)].reset_index(drop=True)

    lig_out = RESULTS / f"02_ligand_set_{a.target}.csv"
    pair_out = RESULTS / f"02_pairs_{a.target}.csv"
    lig.to_csv(lig_out, index=False)
    pairs.to_csv(pair_out, index=False)

    n_cliff = int((pairs.group == "cliff").sum())
    n_ctrl = int((pairs.group == "control").sum())
    print("\n" + "=" * 62)
    print(f"LIGANDS TO SIMULATE : {len(lig)}"
          f"  (cliff {int((lig.group=='cliff').sum())},"
          f" control {int((lig.group=='control').sum())})")
    print(f"PAIRS TO ANALYSE    : {len(pairs)}  (cliff {n_cliff}, control {n_ctrl})")
    print(f"charged ligands     : {int((lig.formal_charge != 0).sum())}/{len(lig)}")
    print(f"ambiguous protonation: {int(lig.prot_ambiguous.sum())} (review before running)")
    print(f"pChEMBL range       : {lig.pchembl.min():.2f} - {lig.pchembl.max():.2f}")
    print(f"\ncliff   |delta pChEMBL|: "
          f"{pairs[pairs.group=='cliff'].delta.min():.2f}"
          f" - {pairs[pairs.group=='cliff'].delta.max():.2f}")
    print(f"control |delta pChEMBL|: "
          f"{pairs[pairs.group=='control'].delta.min():.2f}"
          f" - {pairs[pairs.group=='control'].delta.max():.2f}")
    print(f"\nwritten: {lig_out.name}, {pair_out.name}")

    # simulation cost at the design in the plan
    ns = len(lig) * 2 * 25 * 1.15 + 300
    print(f"\nat 2 replicas x 25 ns: {ns/1000:.1f} us"
          f"  ~= {ns/350:.0f} GPU-days (4 fs HMR)"
          f"  ~= ${ns/350*24*0.34:,.0f} on RunPod 4090")


if __name__ == "__main__":
    main()
