#!/usr/bin/env python
"""Table — what the frozen splits contain, fold by fold.

The splitting is described over three paragraphs and involves four levels, i.e.
ligands, the connected components that are the unit of assignment, the development and
test split, and the five cross-validation folds. A reader following the prose has to
hold that hierarchy in their head while the balancing is explained. The table shows it
once.

The component column is the one that matters and the one a reader would not expect.
Assignment moves components, never ligands, so the component counts are the sample
sizes the balancing actually worked with, and they are small.

Every number is read from the written splits rather than recomputed, so the table
cannot drift from the files the models were scored on.

Writes results/35_splits.tex.

Usage:
    python scripts/table_splits.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, NF = "CHEMBL220", 5


def main() -> None:
    lig = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pair = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")

    def row(name, L, P):
        return (name, L.component.nunique(), len(L), int((L.group == "cliff").sum()),
                len(P), int((P.group == "cliff").sum()), int((P.group != "cliff").sum()))

    dev_l, dev_p = lig[lig.split == "dev"], pair[pair.split == "dev"]
    rows = [row(f"Fold {k}", dev_l[dev_l.cv_fold == k], dev_p[dev_p.cv_fold == k])
            for k in range(NF)]
    dev = row("Development, total", dev_l, dev_p)
    test = row("Test, sealed", lig[lig.split == "test"], pair[pair.split == "test"])
    tot = row("Whole set", lig, pair)

    def line(r):
        return " & ".join([r[0]] + [str(v) for v in r[1:]]) + r" \\"

    tex = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
           r"Split & Components & Ligands & of which cliff & Pairs & Cliff pairs & "
           r"Control pairs \\", r"\midrule"]
    tex += [line(r) for r in rows]
    tex += [r"\midrule", line(dev), line(test), r"\midrule", line(tot),
            r"\bottomrule", r"\end{tabular}"]
    (RES/"35_splits.tex").write_text("\n".join(tex) + "\n")
    print("  -> results/35_splits.tex")
    for r in rows + [dev, test, tot]:
        print(f"    {r[0]:22s} comp {r[1]:3d}  lig {r[2]:3d}  cliff lig {r[3]:3d}  "
              f"pairs {r[4]:3d}  cliff {r[5]:3d}  ctrl {r[6]:3d}")


if __name__ == "__main__":
    main()
