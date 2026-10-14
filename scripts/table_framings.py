#!/usr/bin/env python
"""Table — the three framings of the prediction problem, side by side.

Section 2.5 defines the framings in prose, where the differences between them are
spread over four paragraphs. This puts the unit, the predicted quantity, the row
construction, the metric and the number of algorithm families in one place, so a
reader comparing two numbers later can see in one glance which question each of them
answers.

The design of the table is fixed text, but the two unit counts are not, so they are
read from the frozen splits rather than typed. If the split is ever regenerated the
table follows it, and a table that quietly disagreed with Section 2.4 would be worse
than no table at all.

Writes results/31_framings.tex.

Usage:
    python scripts/table_framings.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET = "CHEMBL220"
# Framing C drops PLS, which has no classification counterpart; see rung4_model.CLF.
N_REG_FAMILIES, N_CLF_FAMILIES = 6, 5


def main() -> None:
    lig = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pair = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    n_lig = int((lig.split == "dev").sum())
    n_pair = int((pair.split == "dev").sum())
    print(f"  development split: {n_lig} ligands, {n_pair} pairs")

    # One math group per cell. Splitting these into several $...$ runs, which the
    # first version did by interpolating the symbols, renders with wrong spacing
    # around the operators.
    rows = [
        ("A", f"Ligand ({n_lig})", "pChEMBL",
         r"$\mathbf{d}_a$, the vector of the ligand itself",
         "$R^2$", N_REG_FAMILIES),
        ("B", f"Pair ({n_pair})", "Signed difference in pChEMBL",
         r"$\mathbf{d}_a-\mathbf{d}_b$ and $\mathbf{d}_b-\mathbf{d}_a$, "
         r"two rows per pair",
         "$R^2$", N_REG_FAMILIES),
        ("C", f"Pair ({n_pair})", "Cliff or matched control",
         r"$[\,\mathbf{d}_a+\mathbf{d}_b,\ |\mathbf{d}_a-\mathbf{d}_b|\,]$, "
         r"one row per pair",
         "ROC AUC", N_CLF_FAMILIES),
    ]
    # >{\raggedright\arraybackslash}p{} rather than a bare p{}, because a justified
    # 2.9 cm column stretches "Signed difference" across the cell. Spelled out rather
    # than hidden behind a \newcommand, since this fragment is \input and a macro
    # definition travelling with it would be defined in the caller's scope.
    col = r">{\raggedright\arraybackslash}p"
    tex = [r"\begin{tabular}{ll" + col + "{2.9cm}" + col + r"{4.3cm}ll}", r"\toprule",
           r"Framing & Unit & Quantity predicted & Features of one row "
           r"& Metric & Algorithms \\", r"\midrule"]
    tex += [" & ".join([f, u, q, feat, m, str(n)]) + r" \\"
            for f, u, q, feat, m, n in rows]
    tex += [r"\bottomrule", r"\end{tabular}"]
    (RES/"31_framings.tex").write_text("\n".join(tex) + "\n")
    print("  -> results/31_framings.tex")


if __name__ == "__main__":
    main()
