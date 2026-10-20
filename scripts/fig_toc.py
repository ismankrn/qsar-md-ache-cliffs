#!/usr/bin/env python
"""Table-of-contents graphic: cost against what it buys.

ACS asks for a print-ready entry of 3.25 by 1.75 inches with sans-serif text, so this
is drawn at that physical size rather than scaled down from a figure panel. Text is
sized to stay legible at it.

The message is the paper's in one picture. Producing a descriptor set costs anywhere
from six milliseconds to eighteen GPU hours per ligand, seven orders of magnitude, and
every rung lands inside the band this design can resolve around the fingerprint
baseline.

Costs and scores are both measured, read from the same files the manuscript tables are
built from. Rung 3 is placed at the rung-2 cost because its input is the rung-2 pose
and the energy evaluation on top of it is not separately timed.

Usage:  python scripts/fig_toc.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/"scripts"))
import table_rungs as tr

RES, FIGS = ROOT/"results", ROOT/"figures"
RESOLUTION = 0.3          # the interval width this design can claim, from stage 17


def main() -> None:
    b = pd.read_csv(RES/"22_benchmark_full.csv")
    best = b[b.framing == "B"].groupby("rung").value.max()
    c = tr.costs()
    dock = c["dock_s"]
    pts = [
        ("2D",      c["fp_ms"]/1000.0, best["Rung 1 2D"],       True),
        ("pose",    dock,              best["Rung 2 pose"],     False),
        ("MM/GBSA", dock,              best["Rung 3 MM/GBSA"],  False),
        ("ensemble", c["ens_s"],       best["Rung 2 ensemble"], False),
        ("MD",      c["md_h"]*3600.0,  best["Rung 4 MD"],       False),
    ]

    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"]})
    fig, ax = plt.subplots(figsize=(3.25, 1.75))

    base = best["Rung 1 2D"]
    ax.axhspan(base - RESOLUTION, base + RESOLUTION, color="#d9d9d9", zorder=0)
    ax.axhline(base, color="#4d4d4d", lw=0.8, ls="--", zorder=1)

    for lab, x, y, is_base in pts:
        ax.scatter([x], [y], s=46 if is_base else 34,
                   c="#b2182b" if is_base else "#2166ac",
                   marker="o" if is_base else "^",
                   edgecolors="white", linewidths=0.6, zorder=3)
        dy = 0.085 if lab not in ("MM/GBSA",) else -0.135
        ax.annotate(lab, (x, y), textcoords="offset points",
                    xytext=(0, 9 if dy > 0 else -13), ha="center",
                    fontsize=6.4, color="#1a1a1a", zorder=4)

    ax.set_xscale("log")
    ax.set_xlim(2e-3, 4e5)
    ax.set_ylim(base - RESOLUTION - 0.06, base + RESOLUTION + 0.10)
    ax.set_xticks([1e-2, 1e0, 1e2, 1e4])
    ax.set_xticklabels(["10 ms", "1 s", "2 min", "3 h"], fontsize=6.4)
    ax.set_yticks([0.0, 0.2, 0.4])
    ax.set_yticklabels(["0.0", "0.2", "0.4"], fontsize=6.4)
    ax.set_xlabel("Cost per ligand", fontsize=7.2, labelpad=1.5)
    ax.set_ylabel("$R^2$ on $\\Delta$pChEMBL", fontsize=7.2, labelpad=1.5)
    ax.tick_params(length=2, pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)

    ax.text(0.985, 0.975,
            "seven orders of magnitude of cost,\nnothing the design can resolve",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=6.8, color="#1a1a1a", linespacing=1.25)
    ax.text(2.6e-3, base - RESOLUTION + 0.022,
            r"band $=\pm0.3$, the resolvable range",
            fontsize=5.8, color="#595959", va="bottom")

    fig.tight_layout(pad=0.25)
    FIGS.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_toc.{ext}", dpi=600)
    w, h = fig.get_size_inches()
    print(f"  -> figures/fig_toc.pdf / .png  ({w:.2f} x {h:.2f} in)")
    for lab, x, y, _ in pts:
        print(f"     {lab:9s} {x:10.3f} s   R2 {y:+.3f}")


if __name__ == "__main__":
    main()
