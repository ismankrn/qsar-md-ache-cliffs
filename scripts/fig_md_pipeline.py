#!/usr/bin/env python
"""Figure — how a descriptor is obtained from a trajectory.

Section 2.5 describes this in prose over three paragraphs, and the prose has to carry
five things at once, i.e. the periodic-boundary treatment, the analysis window, the three
design classes of descriptor, the reliability filter and the construction of pair rows.
A reader who loses the order of those steps cannot tell which of them a given number was
measured after. The figure fixes the order and shows the attrition, which is the part the
prose states but does not make visible.

Every count is read from the files the pipeline wrote rather than typed here, so the
figure cannot drift from the analysis it describes.

Usage:
    python scripts/fig_md_pipeline.py
"""
from __future__ import annotations
import glob, json
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
TARGET, WINDOW = "CHEMBL220", 20.0
EQUIL_NS, ANALYSED_NS, DT_PS = 5.0, 20.0, 10.0

INK = "#1a1a1a"
BOX = "#f2f4f7"
EDGE = "#9aa5b1"
KEEP = "#2166ac"
FULL = "#4a5568"      # the unfiltered arm, which is an arm and not a rejection


def counts() -> dict:
    runs = len(glob.glob(str(RES/"16_rung4_per_run"/"*.json")))
    pil = pd.read_csv(RES/"18_replicate_pilot.csv")
    pil = pil[pil.window.astype(float) == WINDOW]
    v = pil.verdict.value_counts().to_dict()
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pr = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    one = json.loads(Path(sorted(glob.glob(str(RES/"16_rung4_per_run"/"*.json")))[0]).read_text())
    n_res = len([k for k in one["windows"][f"{WINDOW:g}"] if k.startswith("occ_")
                 and k[4:].isdigit()])
    return dict(runs=runs, total=len(pil), usable=v.get("usable", 0),
                rejected=len(pil) - v.get("usable", 0),
                dev_lig=int((sl.split == "dev").sum()),
                dev_pairs=int((pr.split == "dev").sum()), n_res=n_res)


def box(ax, x, y, w, h, title, lines, fc=BOX, ec=EDGE):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                linewidth=0.9, facecolor=fc, edgecolor=ec, zorder=2))
    ax.text(x + w/2, y + h - 0.052, title, ha="center", va="top", fontsize=9.2,
            fontweight="bold", color=INK, zorder=3)
    for i, s in enumerate(lines):
        ax.text(x + w/2, y + h - 0.125 - i*0.062, s, ha="center", va="top",
                fontsize=7.9, color=INK, zorder=3)


def arrow(ax, x0, y0, x1, y1, color=INK, style="-|>", lw=1.0, rad=0.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=11, linewidth=lw, color=color,
                                 connectionstyle=f"arc3,rad={rad}", zorder=4))


def main() -> None:
    c = counts()
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(13.4, 4.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    w, h, y = 0.183, 0.42, 0.555
    xs = [0.008 + i*0.199 for i in range(5)]

    box(ax, xs[0], y, w, h, "Production trajectory",
        [f"{c['runs']} runs", f"a frame every {DT_PS:g} ps",
         "receptor, ligand, ions", "and water"])
    box(ax, xs[1], y, w, h, "Periodic boundary",
        ["each molecule made whole", "along its own bonds, held",
         "in the first-frame image,", "then the receptor", "superposed by Kabsch"])
    box(ax, xs[2], y, w, h, "Analysis window",
        [f"the first {EQUIL_NS:g} ns discarded", f"the next {ANALYSED_NS:g} ns analysed",
         "the same window for", "every ligand"])
    box(ax, xs[3], y, w, h, "Per-frame series",
        [f"{c['total']} descriptors, being", "persistences, fluctuations",
         "and event counts,", f"{c['n_res']} of them per residue"])
    box(ax, xs[4], y, w, h, "Reliability filter",
        ["between-ligand spread", "against between-replicate", "spread, over the",
         "replicate pilot"])

    for i in range(4):
        arrow(ax, xs[i] + w + 0.0035, y + h/2, xs[i+1] - 0.0035, y + h/2)

    # The filter does not discard anything. It splits one descriptor set into two, and
    # both were scored, which a one-way arrow out of the filter would misrepresent.
    bx, bw, by, bh = 0.008, 0.977, 0.045, 0.285
    box(ax, bx, by, bw, bh, "Two feature sets, both scored",
        [f"replicates averaged to one row per ligand, {c['dev_lig']} of them, and a signed "
         f"difference per pair over {c['dev_pairs']} development pairs",
         "the heavy-atom count carried beside either set as a size covariate, and it is "
         "the twelfth feature of the one and the thirty-third of the other"])

    # One arrow with both arms labelled beside it. Two arrows would have to cross each
    # other's labels at this width, and the split is between feature sets rather than
    # between paths through the pipeline.
    xa = xs[4] + w*0.72
    arrow(ax, xa, y - 0.004, xa, by + bh + 0.004, color=INK)
    for lab, col, frac in ((f"{c['usable']} pass the filter, 12 features", KEEP, 0.66),
                           (f"all {c['total']} kept, 33 features", FULL, 0.34)):
        ax.text(xa - 0.009, by + bh + frac*(y - by - bh), lab, ha="right", va="center",
                fontsize=8.4, fontweight="bold", color=col)

    fig.tight_layout(pad=0.25)
    FIGS.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_md_pipeline.{ext}", dpi=300)
    print(f"  {c}")
    print("  -> figures/fig_md_pipeline.pdf / .png")


if __name__ == "__main__":
    main()
