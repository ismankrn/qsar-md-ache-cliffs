#!/usr/bin/env python
"""Why the fold assignment balances two quantities and not one.

The fold assignment in make_splits.py looks more elaborate than it needs to, and the
manuscript owes the reader a reason. This reproduces the alternatives on the frozen
development components and shows what each one gives.

It exists because the reason was first quoted from comments in make_splits.py, which
record folds of "25 cliff pairs against 4" and "13 vs 63 ligands". Those numbers are
historical, from an earlier state of the pipeline, and no variant reproduces them
against the frozen splits. The joint criterion reproduces exactly, which is what shows
this harness is faithful, so the failure cases are quoted from here instead.

Writes results/34_fold_balance.csv.

Usage:
    python scripts/check_fold_balance.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, NF = "CHEMBL220", 5


def components() -> pd.DataFrame:
    """One row per development component, with its load on both quantities."""
    lig = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pair = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig, pair = lig[lig.split == "dev"], pair[pair.split == "dev"]
    rows = []
    for c, g in lig.groupby("component"):
        ids = set(g.molecule_chembl_id)
        rows.append(dict(component=c, n_lig=len(g),
                         n_cliff=int(((pair.group == "cliff")
                                      & pair.ligand_a.isin(ids)).sum())))
    # largest first, which is what greedy bin packing needs to work at all
    return pd.DataFrame(rows).sort_values(["n_cliff", "n_lig"], ascending=False)


def pack(d: pd.DataFrame, rule: str) -> dict:
    tgt_c = max(d.n_cliff.sum()/NF, 1)
    tgt_l = max(d.n_lig.sum()/NF, 1)
    load = {k: {"lig": 0, "cliff": 0} for k in range(NF)}
    for _, r in d.iterrows():
        if rule == "joint":
            k = min(load, key=lambda f: load[f]["cliff"]/tgt_c + load[f]["lig"]/tgt_l)
        else:
            k = min(load, key=lambda f: load[f][rule])
        load[k]["lig"] += r.n_lig
        load[k]["cliff"] += r.n_cliff
    return load


def main() -> None:
    d = components()
    print(f"  {len(d)} development components, {d.n_lig.sum()} ligands, "
          f"{d.n_cliff.sum()} cliff pairs\n")
    rows = []
    for rule, label in (("cliff", "cliff pairs only"), ("lig", "ligand count only"),
                        ("joint", "both, as used")):
        load = pack(d, rule)
        L = [load[k]["lig"] for k in range(NF)]
        C = [load[k]["cliff"] for k in range(NF)]
        rows.append(dict(rule=label, lig_min=min(L), lig_max=max(L),
                         cliff_min=min(C), cliff_max=max(C)))
        print(f"  {label:20s} ligands {min(L):3d}--{max(L):3d}   "
              f"cliff pairs {min(C):3d}--{max(C):3d}")
    pd.DataFrame(rows).to_csv(RES/"34_fold_balance.csv", index=False)
    print("\n  -> results/34_fold_balance.csv")


if __name__ == "__main__":
    main()
