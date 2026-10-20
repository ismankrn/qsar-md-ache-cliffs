#!/usr/bin/env python
"""Supporting-information table of every molecular structure used.

The JCIM editorial that governs QSAR/QSPR submissions (Jorgensen, J. Chem. Inf. Model.
2006, 46, 937) requires that all data and molecular structures used to carry out the
study be reported in the paper or its supporting information, or be readily available
without restriction. The deposited archive satisfies the second route, and this table
satisfies the first, which is the one the editorial names first.

One row per ligand: the ChEMBL identifier, the standardised neutral SMILES that carries
the fingerprints and the pair definition, the protonated SMILES that enters the docking
and the simulations, the measured potency, the split and fold it was frozen into, and
whether it belongs to a cliff pair.

Usage:  python scripts/table_si_structures.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET = "CHEMBL220"


def esc(s: str) -> str:
    """SMILES into a breakable \\texttt cell.

    Order matters. The break points go in first, parked on a byte that cannot occur in
    a SMILES, because inserting them after the backslash has become \\textbackslash{}
    would split that macro. The backslash is parked for the same reason while the other
    LaTeX specials are escaped. seqsplit is not used: it re-reads its argument and the
    escaped # of a triple bond comes back as a macro parameter.
    """
    out = str(s)
    for ch in ")(.=/]":
        out = out.replace(ch, ch + "\x00")
    out = out.replace("\\", "\x01")
    for ch in "#$%&_~^{}":
        out = out.replace(ch, "\\" + ch)
    out = out.replace("\x01", r"\textbackslash{}")
    return out.replace("\x00", r"\allowbreak{}")


def main() -> None:
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
    spl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")[
        ["molecule_chembl_id", "split", "cv_fold"]]
    d = lig.merge(spl, on="molecule_chembl_id").sort_values("molecule_chembl_id")

    rows = []
    for _, r in d.iterrows():
        where = "test" if r.split == "test" else f"fold {int(r.cv_fold)}"
        rows.append(" & ".join([
            r.molecule_chembl_id,
            r"\texttt{" + esc(r.smiles_std) + "}",
            r"\texttt{" + esc(r.smiles_md) + "}",
            f"{r.pchembl:.2f}",
            f"{int(r.formal_charge):+d}",
            where,
            "cliff" if r.group == "cliff" else "control",
        ]) + r" \\")

    body = "\n".join(rows)
    out = (
        r"\begin{longtable}{l p{4.4cm} p{4.4cm} r r l l}" "\n"
        r"\caption{\textbf{Every ligand used in this study.} The standardised neutral"
        "\n"
        r"  structure carries the fingerprints and the cliff definition; the protonated"
        "\n"
        r"  structure enters the docking and the simulations. Potency is pChEMBL, the"
        "\n"
        r"  negative logarithm of the molar IC\textsubscript{50}. The split and fold were frozen"
        "\n"
        r"  before any descriptor was computed.}\label{tab:sistructures}\\" "\n"
        r"\toprule" "\n"
        r"ChEMBL ID & SMILES, neutral & SMILES, pH 7.4 & pChEMBL & $q$ & Split & Group \\"
        "\n" r"\midrule" "\n" r"\endfirsthead" "\n"
        r"\multicolumn{7}{l}{\itshape Table \thetable{} continued}\\" "\n"
        r"\toprule" "\n"
        r"ChEMBL ID & SMILES, neutral & SMILES, pH 7.4 & pChEMBL & $q$ & Split & Group \\"
        "\n" r"\midrule" "\n" r"\endhead" "\n"
        r"\midrule \multicolumn{7}{r}{\itshape continued on the next page}\\" "\n"
        r"\endfoot" "\n" r"\bottomrule" "\n" r"\endlastfoot" "\n"
        + body + "\n" r"\end{longtable}" "\n"
    )
    (RES/"38_si_structures.tex").write_text(out)
    print(f"  {len(d)} ligands, {(d.group=='cliff').sum()} in cliff pairs, "
          f"{(d.split=='test').sum()} sealed")
    print("  -> results/38_si_structures.tex")


if __name__ == "__main__":
    main()
