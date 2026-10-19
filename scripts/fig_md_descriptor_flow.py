#!/usr/bin/env python
"""Figure — how the rung-4 descriptors are obtained from a trajectory.

A corrected redraw of the supplied flowchart. The shape and the colour scheme are kept;
what changed is the content, which disagreed with the pipeline in four places.

  class sizes     19 persistence, 8 fluctuation, 3 event counts, and 2 that sit between
                  classes, rather than an even 10 / 12 / 10. Fourteen of the 32 are
                  per-residue contact occupancies, which is why persistence dominates,
                  and that imbalance is the finding rather than an untidiness to smooth.
  examples        hydrogen-bond persistence, pi-pi stacking, B-factor-like metrics,
                  hydrogen-bond events and contact events are not computed in this work.
  preprocessing   the periodic-boundary treatment and the Kabsch superposition. There is
                  no handling of missing data.
  missing stages  the analysis window, the panel of fourteen residues and the reliability
                  filter were absent, and each decides what the descriptors mean.

Drawn at the printed size and written as a vector PDF, so no raster resolution applies.
Counts are read from the files the pipeline wrote.

Usage:
    python scripts/fig_md_descriptor_flow.py
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

W_IN, H_IN = 6.30, 7.30          # 16 cm wide, printed 1:1
PT = 1/(H_IN*72)

BLUE = ("#c9dcf2", "#41618f")
PURPLE = ("#dccdef", "#6f55949")[0], "#6f5594"
GREEN = ("#d8e8c2", "#5f8a37")
TEAL = ("#c7e6e8", "#3f8d94")
ORANGE = ("#f7d6b4", "#c87a35")
ENDG = ("#cde7b8", "#5f8a37")
INK = "#15202b"
BODY = "#2f3b47"

# The three classes, with the descriptors each holds. Named here so the figure and the
# counts beside them cannot disagree.
PERSIST = ["occ_cas", "occ_pas", "occ_bridge", "bridge_longest_ns",
           "frac_lig_buried_mean"] + [f"occ_{r}" for r in
           (72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447)]
FLUCT = ["rmsd_mean", "rmsd_max", "depth_sd", "rmsf_lig", "rmsf_gorge",
         "gate_chi1_sd", "lig_dihedral_sd", "ca_rmsd_mean"]
EVENT = ["bridge_runs_per_ns", "gate_flips_per_ns", "escaped"]
BETWEEN = ["depth_mean", "n_bridge_water"]


def counts() -> dict:
    runs = len(glob.glob(str(RES/"16_rung4_per_run"/"*.json")))
    pil = pd.read_csv(RES/"18_replicate_pilot.csv")
    pil = pil[pil.window.astype(float) == WINDOW]
    scored = set(pil.descriptor)
    mine = set(PERSIST) | set(FLUCT) | set(EVENT) | set(BETWEEN)
    if mine != scored:
        raise SystemExit(f"classes do not cover the scored set\n"
                         f"  missing: {sorted(scored - mine)}\n"
                         f"  extra:   {sorted(mine - scored)}")
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pr = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    return dict(runs=runs, total=len(scored),
                usable=int((pil.verdict == "usable").sum()),
                dev_lig=int((sl.split == "dev").sum()),
                dev_pairs=int((pr.split == "dev").sum()))


def band(ax, y, h, fill, edge, title, sub=None, x=0.10, w=0.80, fs=9.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
                                linewidth=1.0, facecolor=fill, edgecolor=edge, zorder=2))
    yt = y + h/2 + (5.5*PT if sub else 0)
    ax.text(x + w/2, yt, title, ha="center", va="center", fontsize=fs,
            fontweight="bold", color=INK, zorder=3)
    if sub:
        ax.text(x + w/2, y + h/2 - 6.5*PT, sub, ha="center", va="center",
                fontsize=fs - 1.3, color=BODY, zorder=3)


def pill(ax, y, h, fill, edge, text):
    ax.add_patch(FancyBboxPatch((0.37, y), 0.26, h,
                                boxstyle="round,pad=0,rounding_size=0.014",
                                linewidth=1.0, facecolor=fill, edgecolor=edge, zorder=2))
    ax.text(0.50, y + h/2, text, ha="center", va="center", fontsize=9.5,
            fontweight="bold", color=INK, zorder=3)


def column(ax, x, w, y, h, fill, edge, head, n, gist, items):
    """Header, one line of gist, then the examples, as in the supplied chart."""
    light = fill + "80"
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
                                linewidth=1.0, facecolor="white", edgecolor=edge, zorder=2))
    hh = 26*PT
    ax.add_patch(FancyBboxPatch((x, y + h - hh), w, hh,
                                boxstyle="round,pad=0,rounding_size=0.012",
                                linewidth=0, facecolor=fill, zorder=3))
    ax.text(x + w/2, y + h - 10*PT, head, ha="center", va="center", fontsize=8.4,
            fontweight="bold", color=INK, zorder=4)
    ax.text(x + w/2, y + h - 20*PT, f"({n} descriptors)", ha="center", va="center",
            fontsize=7.6, color=INK, zorder=4)
    ax.add_patch(FancyBboxPatch((x, y + h - hh - 20*PT), w, 20*PT,
                                boxstyle="round,pad=0,rounding_size=0.002",
                                linewidth=0, facecolor=light, zorder=3))
    ax.text(x + w/2, y + h - hh - 10*PT, gist, ha="center", va="center", fontsize=7.2,
            color=BODY, zorder=4, linespacing=1.35)
    ax.text(x + 0.016, y + h - hh - 30*PT, "Examples", ha="left", va="top",
            fontsize=7.4, fontweight="bold", color=INK, zorder=4)
    for i, s in enumerate(items):
        ax.text(x + 0.020, y + h - hh - 42*PT - i*10.2*PT, "•", ha="left", va="top",
                fontsize=7.0, color=edge, zorder=4)
        ax.text(x + 0.036, y + h - hh - 42*PT - i*10.2*PT, s, ha="left", va="top",
                fontsize=7.0, color=BODY, zorder=4)


def arrow(ax, x0, y0, x1, y1, lw=1.1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9,
                                 linewidth=lw, color="#3b4654", zorder=5))


def elbow(ax, x0, y0, x1, y1, ymid):
    ax.plot([x0, x0], [y0, ymid], color="#3b4654", lw=1.1, zorder=5)
    ax.plot([x0, x1], [ymid, ymid], color="#3b4654", lw=1.1, zorder=5)
    arrow(ax, x1, ymid, x1, y1)


def main() -> None:
    c = counts()
    plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"]})
    fig, ax = plt.subplots(figsize=(W_IN, H_IN))
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    pill(ax, 0.950, 0.030, *BLUE, "Start")
    band(ax, 0.876, 0.052, *BLUE, "MD trajectory",
         f"{c['runs']} runs, a frame every 10 ps")
    band(ax, 0.784, 0.070, *PURPLE, "Trajectory preprocessing",
         "molecules made whole along their bonds, kept from jumping between\n"
         "frames, held in the first-frame image, receptor superposed by Kabsch")
    band(ax, 0.710, 0.052, *PURPLE, "Analysis window",
         "the first 5 ns discarded, the next 20 ns analysed")

    for y0, y1 in ((0.950, 0.930), (0.876, 0.856), (0.784, 0.764)):
        arrow(ax, 0.50, y0, 0.50, y1)

    cw, cy, ch = 0.300, 0.392, 0.258
    cxs = [0.016, 0.350, 0.684]
    cols = [
        (GREEN, "Persistence descriptors", len(PERSIST),
         "how much of the window a contact\nor a position is held",
         ["occupancy at Trp86 and Trp286", "bridging occupancy",
          "longest bridging episode", "buried fraction of the ligand",
          "14 per-residue occupancies"]),
        (TEAL, "Fluctuation descriptors", len(FLUCT),
         "how widely something moves\nwithin the window",
         ["ligand RMSD, mean and maximum", "ligand and gorge RMSF",
          "Tyr337 $\\chi_1$ spread", "ligand dihedral spread",
          "depth spread along the axis", "receptor C$\\alpha$ RMSD"]),
        (ORANGE, "Event counts", len(EVENT),
         "how often something happens\nover the window",
         ["bridging episodes per ns", "Tyr337 rotamer flips per ns",
          "escape from the gorge"]),
    ]
    for x, (col, head, n, gist, items) in zip(cxs, cols):
        column(ax, x, cw, cy, ch, col[0], col[1], head, n, gist, items)
        elbow(ax, 0.50, 0.710, x + cw/2, cy + ch, 0.672)
        arrow(ax, x + cw/2, cy, x + cw/2, 0.300)

    band(ax, 0.244, 0.056, *PURPLE, "Combine results",
         f"{c['total']} descriptors in total, with {' and '.join(BETWEEN)}\n"
         "lying between the classes")
    arrow(ax, 0.50, 0.244, 0.50, 0.222)
    band(ax, 0.154, 0.068, *ORANGE, "Reliability filter",
         "the between-ligand spread at least twice the between-replicate\n"
         f"spread, which {c['usable']} of the {c['total']} meet")
    arrow(ax, 0.50, 0.154, 0.50, 0.132)
    band(ax, 0.062, 0.070, *GREEN, "Two feature sets, both scored",
         f"12 features from the {c['usable']} that pass and 33 from all "
         f"{c['total']}, each with\nthe heavy-atom count, over {c['dev_lig']} ligands "
         f"and {c['dev_pairs']} development pairs")
    arrow(ax, 0.50, 0.062, 0.50, 0.042)
    pill(ax, 0.008, 0.030, *ENDG, "End")

    FIGS.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_md_descriptor_flow.{ext}", dpi=400)
    print(f"  {c}")
    print(f"  classes: {len(PERSIST)} / {len(FLUCT)} / {len(EVENT)} "
          f"+ {len(BETWEEN)} between = {len(PERSIST)+len(FLUCT)+len(EVENT)+len(BETWEEN)}")
    print("  -> figures/fig_md_descriptor_flow.pdf / .png")


if __name__ == "__main__":
    main()
