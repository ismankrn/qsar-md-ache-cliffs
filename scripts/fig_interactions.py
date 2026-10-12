#!/usr/bin/env python
"""Figure — the reference ligand's interactions, and whether docking reproduces them.

Donepezil in its crystal pose (4EY7) is the reference this study docks against, and the
conventional way to show a binding mode is an interaction diagram. This draws one, and then
makes it earn its place: the same fingerprint is computed for the redocked pose, so the
figure says not only what donepezil touches but which of those contacts the docking protocol
recovers. Figure 2 validates docking by RMSD; this validates it by interaction, which is the
level the rung-2 and rung-4 descriptors actually work at.

Interaction typing is ProLIF's, not ours. The rung-2 and rung-4 descriptors are distance-based
occupancies with no hydrogen-bond or pi-stacking distinction, so an independent tool is the
right thing here -- it shows what the in-house descriptors are summarising and what they throw
away. This is descriptive: nothing from this figure enters any model, and adding a typed
descriptor family after the ladder was scored would be post hoc.

The crystal ligand PDB carries no bond orders, so they are transferred from the docked SDF of
the same molecule with AssignBondOrdersFromTemplate. The two poses then share one chemistry
and differ only in coordinates, which is what makes the comparison meaningful.

Usage:
    python scripts/fig_interactions.py
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")
import MDAnalysis as mda
import prolif
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

ROOT = Path(__file__).resolve().parent.parent
DOCK, RES, FIGS = ROOT/"docking", ROOT/"results", ROOT/"figures"

REC = DOCK/"4ey7_rec_clean.pdb"
XTAL = DOCK/"4ey7_ligand_ref.pdb"
DOCKED = DOCK/"E20_docked_pose1.sdf"

# Residues the study names, so the figure and the descriptors speak of the same places.
NAMED = {86: "CAS", 286: "PAS", 337: "gate", 203: "triad", 447: "triad", 334: "triad"}
COLOUR = {
    "Hydrophobic": "#8c8c8c", "HBAcceptor": "#2166ac", "HBDonor": "#4393c3",
    "PiStacking": "#1b7837", "EdgeToFace": "#1b7837", "FaceToFace": "#1b7837",
    "CationPi": "#b2182b", "PiCation": "#b2182b", "Cationic": "#d6604d",
    "Anionic": "#d6604d", "VdWContact": "0.8", "XBAcceptor": "#762a83",
    "XBDonor": "#762a83",
}


def ligand_from(path: Path, template: Chem.Mol) -> Chem.Mol:
    """Crystal coordinates with the template's bond orders."""
    raw = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)
    return AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(template), raw)


def interactions(lig: Chem.Mol, prot) -> pd.DataFrame:
    fp = prolif.Fingerprint()
    fp.run_from_iterable([prolif.Molecule.from_rdkit(lig)], prot, progress=False)
    df = fp.to_dataframe()
    out = []
    for (_lig, resid, kind), col in df.items():
        if bool(col.iloc[0]):
            out.append(dict(residue=str(resid), kind=kind))
    return pd.DataFrame(out)


def resnum(label: str) -> int:
    digits = "".join(ch for ch in label if ch.isdigit())
    return int(digits) if digits else -1


def main() -> None:
    tmpl = Chem.MolFromMolFile(str(DOCKED), removeHs=False)
    if tmpl is None:
        sys.exit(f"cannot read {DOCKED}")
    lig_x = ligand_from(XTAL, tmpl)
    lig_d = Chem.RemoveHs(tmpl)

    # The docking receptor is protonated (4,008 H of 8,134 atoms), which H-bond detection
    # needs, but the PDB carries no CONECT records so MDAnalysis sees no bonds at all and
    # ProLIF reports every residue as fragmented. Guessing them from distances is the fix.
    u = mda.Universe(str(REC))
    prot_ag = u.select_atoms("protein")
    prot_ag.guess_bonds()
    prot = prolif.Molecule.from_mda(prot_ag)

    ix = interactions(lig_x, prot).assign(pose="crystal")
    idk = interactions(lig_d, prot).assign(pose="docked")
    both = pd.concat([ix, idk], ignore_index=True)
    both["resnum"] = both.residue.map(resnum)
    both["site"] = both.resnum.map(NAMED).fillna("")
    both.to_csv(RES/"24_reference_interactions.csv", index=False)
    print(f"  crystal: {len(ix)} interactions, docked: {len(idk)}")

    # presence matrix, residues x (pose, kind)
    key = both.assign(one=1).pivot_table(index=["resnum", "residue", "site"],
                                         columns=["kind", "pose"], values="one",
                                         aggfunc="max").fillna(0)
    key = key.sort_index(level="resnum")

    FIGS.mkdir(exist_ok=True)
    fig = plt.figure(figsize=(12.0, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.25], wspace=0.22)

    # ---- (a) the ligand and what it touches ----
    ax = fig.add_subplot(gs[0, 0])
    AllChem.Compute2DCoords(lig_x)
    img = Draw.MolToImage(lig_x, size=(620, 620), wedgeBonds=False)
    ax.imshow(np.asarray(img))
    ax.axis("off")
    xres = ix.copy()
    xres["resnum"] = xres.residue.map(resnum)
    uniq = xres.drop_duplicates("residue").sort_values("resnum")
    n = len(uniq)
    for i, (_, r) in enumerate(uniq.iterrows()):
        ang = 2*np.pi*i/max(n, 1) - np.pi/2
        x, y = 310 + 360*np.cos(ang), 310 + 360*np.sin(ang)
        kinds = sorted(set(xres[xres.residue == r.residue].kind))
        col = COLOUR.get(kinds[0], "0.4")
        site = NAMED.get(r.resnum, "")
        lab = r.residue + (f"\n{site}" if site else "")
        ax.annotate(lab, xy=(310, 310), xytext=(x, y), fontsize=7, ha="center",
                    va="center", color=col,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=col, lw=0.9),
                    arrowprops=dict(arrowstyle="-", lw=0.9, color=col, alpha=0.75,
                                    shrinkA=18, shrinkB=10))
    ax.set_xlim(-130, 750); ax.set_ylim(750, -130)
    ax.set_title("(a) Donepezil in the crystal pose and the residues it contacts",
                 fontsize=9.5)
    # The edge colours carry the interaction type and were unlabelled in the first draft.
    seen = []
    for k in sorted(set(xres.kind)):
        c = COLOUR.get(k, "0.4")
        if (k, c) not in seen:
            seen.append((k, c))
    ax.legend(handles=[Line2D([], [], color=c, lw=2, label=k) for k, c in seen],
              fontsize=7, frameon=False, loc="lower center", ncol=3,
              bbox_to_anchor=(0.5, -0.06))

    # ---- (b) crystal vs docked, by residue and interaction type ----
    ax = fig.add_subplot(gs[0, 1])
    rows = list(key.index)
    kinds = sorted({k for k, _p in key.columns})
    ypos = {r: i for i, r in enumerate(rows)}
    for j, kind in enumerate(kinds):
        for r in rows:
            c = key.loc[r]
            xo = j
            for pose, dx, mk in (("crystal", -0.16, "o"), ("docked", 0.16, "s")):
                v = c.get((kind, pose), 0)
                if v:
                    ax.scatter(xo + dx, ypos[r], marker=mk, s=46,
                               color=COLOUR.get(kind, "0.4"),
                               edgecolors="0.25", linewidths=0.5)
    ax.set_xticks(range(len(kinds)))
    ax.set_xticklabels(kinds, rotation=30, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{r[1]}" + (f"  [{r[2]}]" if r[2] else "") for r in rows],
                       fontsize=7)
    ax.set_ylim(-0.8, len(rows)-0.2)
    ax.grid(color="0.93", lw=0.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("(b) Does the docked pose recover the crystal interactions?",
                 fontsize=9.5)
    # Outside the axes: the first version sat on the TYR72 row and hid its only marker.
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color="0.35", label="crystal"),
                       Line2D([], [], marker="s", ls="", color="0.35", label="docked")],
              fontsize=7.5, frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.02),
              ncol=2)

    only_x = set(ix.residue) - set(idk.residue)
    only_d = set(idk.residue) - set(ix.residue)
    shared = set(ix.residue) & set(idk.residue)
    fig.suptitle("Reference binding mode, and docking checked at the level of contacts",
                 fontsize=11, y=0.98)
    fig.text(0.5, 0.012,
             f"Interaction typing by ProLIF {prolif.__version__}; bond orders transferred "
             f"from the docked SDF so both poses share one chemistry. "
             f"{len(shared)} residues contacted in both poses, {len(only_x)} in the crystal "
             f"only, {len(only_d)} in the docked pose only. Descriptive: nothing here enters "
             f"any model.", ha="center", fontsize=7, color="0.35")
    fig.tight_layout(rect=(0, 0.045, 1, 0.95))
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_interactions.{ext}", dpi=300)
    print(f"  -> figures/fig_interactions.pdf / .png")
    print(f"  -> results/24_reference_interactions.csv")
    print(f"  shared {sorted(shared)}")
    print(f"  crystal only {sorted(only_x)}   docked only {sorted(only_d)}")


if __name__ == "__main__":
    main()
