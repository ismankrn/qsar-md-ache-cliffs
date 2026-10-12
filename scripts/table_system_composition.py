#!/usr/bin/env python
"""Table — simulated system composition, from the build records of all 217 systems.

A Methods section that describes a molecular dynamics campaign without saying how large the
systems were, what was in the box, or which force fields were combined is incomplete, and
these numbers were only ever in per-system JSON on the cluster. This collects them.

It also carries two checks worth stating rather than assuming. Every system should be
electroneutral after counterions, and the build records the residual charge, so the worst
case across 217 systems is a number the manuscript can quote instead of a claim. And the
receptor-ligand minimum distance at build time should be a real contact rather than a clash,
which is what the ligand-placement step was supposed to guarantee.

Writes results/27_system_composition.csv and a LaTeX fragment.

Usage:
    python scripts/table_system_composition.py --builds <dir with */build.json>
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"


def esc(s) -> str:
    return str(s).replace("_", r"\_").replace("&", r"\&")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--builds", required=True, help="directory containing */build.json")
    a = ap.parse_args()

    rows = []
    for f in sorted(Path(a.builds).rglob("build.json")):
        d = json.loads(f.read_text())
        box = d.get("box_A") or [np.nan]*3
        rows.append(dict(
            ligand=d.get("ligand"), n_atoms=d.get("n_atoms"),
            n_receptor=d.get("n_receptor"), n_ligand=d.get("n_ligand"),
            box_A=float(box[0]),
            ligand_charge=d.get("ligand_formal_charge"),
            residual_charge=abs(float(d.get("system_charge", np.nan))),
            min_dist_A=d.get("min_dist_receptor_ligand_A"),
            ff_protein=(d.get("forcefield") or {}).get("protein"),
            ff_ligand=(d.get("forcefield") or {}).get("ligand"),
            ff_water=(d.get("forcefield") or {}).get("water")))
    d = pd.DataFrame(rows)
    if d.empty:
        raise SystemExit("no build.json found")
    # solvent is whatever is not receptor or ligand; TIP3P is 3 sites per molecule
    d["n_solvent_atoms"] = d.n_atoms - d.n_receptor - d.n_ligand
    d.to_csv(RES/"27_system_composition.csv", index=False)
    print(f"  {len(d)} systems -> results/27_system_composition.csv")

    def rng(col, fmt="{:.0f}"):
        s = d[col].dropna()
        if s.nunique() == 1:
            return fmt.format(s.iloc[0])
        return f"{fmt.format(s.min())}--{fmt.format(s.max())}"

    ff = {c: sorted(set(d[c].dropna())) for c in ("ff_protein", "ff_ligand", "ff_water")}
    for k, v in ff.items():
        if len(v) != 1:
            print(f"  WARNING: {k} is not uniform across systems: {v}")

    items = [
        ("Systems built", f"{len(d)}"),
        ("Atoms per system", rng("n_atoms")),
        ("Receptor atoms", rng("n_receptor")),
        ("Ligand atoms", rng("n_ligand")),
        ("Solvent atoms", rng("n_solvent_atoms")),
        ("Cubic box edge (\\AA)", rng("box_A", "{:.1f}")),
        ("Ligand formal charge", rng("ligand_charge", "{:+.0f}")),
        ("Receptor--ligand minimum distance (\\AA)", rng("min_dist_A", "{:.2f}")),
        ("Residual system charge ($e$), worst case",
         f"{d.residual_charge.max():.1e}"),
        ("Protein force field", esc(ff['ff_protein'][0]) if ff['ff_protein'] else "--"),
        ("Ligand force field", esc(ff['ff_ligand'][0]) if ff['ff_ligand'] else "--"),
        ("Water model", esc(ff['ff_water'][0]) if ff['ff_water'] else "--"),
    ]
    tex = [r"\begin{tabular}{ll}", r"\toprule", r"Quantity & Value \\", r"\midrule"]
    tex += [f"{k} & {v}" + r" \\" for k, v in items]
    tex += [r"\bottomrule", r"\end{tabular}"]
    (RES/"27_system_composition.tex").write_text("\n".join(tex) + "\n")
    print(f"  -> results/27_system_composition.tex")
    for k, v in items:
        print(f"    {k:44s} {v}")


if __name__ == "__main__":
    main()
