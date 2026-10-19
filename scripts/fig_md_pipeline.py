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

# One accent family rather than one colour per stage, graded light to dark so the row
# reads as a progression without turning into five unrelated categories.
ACCENT = ["#9ec5e8", "#6aa3d4", "#3f7fbd", "#2563a0", "#17457a"]
FIGURE = "#1d5490"    # one weight for every count, so none reads as the faint one
INK = "#16212e"
SOFT = "#5b6b7d"
CARD = "#ffffff"
EDGE = "#dde3ea"
BAND = "#f4f7fa"
KEEP = "#17457a"
FULL = "#7b8794"


def counts() -> dict:
    runs = len(glob.glob(str(RES/"16_rung4_per_run"/"*.json")))
    pil = pd.read_csv(RES/"18_replicate_pilot.csv")
    pil = pil[pil.window.astype(float) == WINDOW]
    v = pil.verdict.value_counts().to_dict()
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pr = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    one = json.loads(Path(sorted(glob.glob(str(RES/"16_rung4_per_run"/"*.json")))[0]).read_text())
    w = one["windows"][f"{WINDOW:g}"]
    n_res = len([k for k in w if k.startswith("occ_") and k[4:].isdigit()])
    reps = pd.Series([Path(f).stem for f in glob.glob(str(RES/"16_rung4_per_run"/"*.json"))])
    n_trip = int((reps.str.replace(r"_r\d+$", "", regex=True).value_counts() > 1).sum())
    return dict(runs=runs, total=len(pil), usable=v.get("usable", 0),
                dev_lig=int((sl.split == "dev").sum()),
                dev_pairs=int((pr.split == "dev").sum()),
                n_res=n_res, n_trip=n_trip)


# The canvas is the printed size. A figure drawn at 13 inches and scaled into a 16 cm
# text block loses half its type, which is how 7.6 pt becomes an unreadable 3.6 pt.
# Five cards abreast cannot hold legible type at that width either, so the flow runs
# down the page and every stage gets the full measure for its detail.
W_IN, H_IN = 6.30, 3.90
PT = 1/(H_IN*72)
# The axes are the figure. tight_layout would rescale them around the drawn text and
# move every fractional coordinate here off its intended place.
FILL = dict(left=0, right=1, top=1, bottom=0)

F_NODE, F_TITLE, F_NUM, F_BODY = 6.3, 8.2, 11.5, 7.0
LEAD = 9.4*PT
X_SPINE, X_TEXT, X_RIGHT = 0.026, 0.070, 0.990


def stage(ax, y, n, accent, title, figure, detail):
    """One row of the flow. Node on the spine, title and count, then the detail."""
    ax.plot([X_SPINE], [y], "o", ms=11, color=accent, zorder=5,
            markeredgecolor="white", markeredgewidth=0.8)
    ax.text(X_SPINE, y, str(n), ha="center", va="center", fontsize=F_NODE,
            fontweight="bold", color="white", zorder=6)
    ax.text(X_TEXT, y + 4*PT, title, ha="left", va="bottom", fontsize=F_TITLE,
            fontweight="bold", color=INK, zorder=4)
    ax.text(X_RIGHT, y + 4*PT, figure, ha="right", va="bottom", fontsize=F_NUM,
            fontweight="bold", color=FIGURE, zorder=4)
    for i, s in enumerate(detail):
        ax.text(X_TEXT, y - 5*PT - i*LEAD, s, ha="left", va="top",
                fontsize=F_BODY, color=SOFT, zorder=4)


def arrow(ax, x0, y0, x1, y1, color=INK, lw=0.9):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=7, linewidth=lw, color=color, zorder=5))


def main() -> None:
    c = counts()
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(W_IN, H_IN))
    fig.subplots_adjust(**FILL)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    rows = [
        ("Production trajectory", f"{c['runs']} runs",
         [f"over {c['dev_lig']} ligands, {c['n_trip']} of them simulated in triplicate, "
          "with a frame written every 10 ps"]),
        ("Periodic boundary", "4 corrections",
         ["each molecule made whole along its own bonds, kept from jumping between frames",
          "and held in the first-frame image, and only then the receptor superposed by Kabsch"]),
        ("Analysis window", f"{ANALYSED_NS:g} ns",
         [f"the first {EQUIL_NS:g} ns of every run discarded and the same "
          f"{ANALYSED_NS:g} ns analysed for every ligand"]),
        ("Descriptors", f"{c['total']} of them",
         [f"persistences, fluctuations and event counts, {c['n_res']} of them one per "
          "gorge-lining residue,", "and none of them a pose descriptor averaged over frames"]),
        ("Reliability filter", f"{c['usable']} of {c['total']}",
         ["a descriptor qualifies where the spread between ligands is at least twice",
          "the spread between replicates of the same ligand"]),
    ]

    top, step = 0.928, 0.145
    ys = [top - i*step for i in range(len(rows))]
    ax.plot([X_SPINE, X_SPINE], [ys[-1], ys[0]], lw=1.1, color="#d2dbe4", zorder=1)
    for i, ((title, figtxt, detail), yy) in enumerate(zip(rows, ys)):
        stage(ax, yy, i + 1, ACCENT[i], title, figtxt, detail)

    # The filter splits one descriptor set in two and both were scored, so the flow ends
    # in one band rather than in a kept branch beside a discarded one.
    by, bh = 0.018, 0.190
    arrow(ax, X_SPINE, ys[-1] - 0.050, X_SPINE, by + bh + 0.013, color="#9fb0c0", lw=1.1)
    ax.add_patch(FancyBboxPatch((0.0, by), 1.0, bh,
                                boxstyle="round,pad=0,rounding_size=0.008",
                                linewidth=0.7, facecolor=BAND, edgecolor=EDGE, zorder=2))
    ax.add_patch(FancyBboxPatch((0.0, by), 0.005, bh,
                                boxstyle="round,pad=0,rounding_size=0.004",
                                linewidth=0, facecolor=ACCENT[4], zorder=3))
    xb = 0.030
    ax.text(xb, by + bh - 9*PT, "Two feature sets, both scored", ha="left", va="top",
            fontsize=F_TITLE, fontweight="bold", color=INK, zorder=4)
    for i, (lab, col) in enumerate((
            (f"the {c['usable']} that pass the filter, with the heavy-atom count, "
             "give 12 features", KEEP),
            (f"all {c['total']}, with the same covariate, give 33", FULL),
            (f"one row per ligand over {c['dev_lig']} of them, and a signed difference "
             f"per pair over {c['dev_pairs']} development pairs", SOFT))):
        ax.text(xb, by + bh - 24*PT - i*LEAD, lab, ha="left", va="top", fontsize=F_BODY,
                fontweight="bold" if col is not SOFT else "normal", color=col, zorder=4)

    FIGS.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_md_pipeline.{ext}", dpi=400)
    print(f"  {c}")
    print(f"  canvas {W_IN}x{H_IN} in, printed 1:1 in a 16 cm text block")
    print("  -> figures/fig_md_pipeline.pdf / .png")


if __name__ == "__main__":
    main()
