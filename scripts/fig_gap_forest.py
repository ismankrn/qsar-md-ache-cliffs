#!/usr/bin/env python
"""Figure — the gap between rung 4 and rung 1, with its interval, on the three framings.

Section 3.4 states three numbers and three intervals, and the claim they carry is about
one thing only, namely whether an interval clears zero. A forest plot answers that at a
glance where a sentence makes the reader hold three pairs of bounds in their head.

This plots a DIFFERENCE rather than a level, which is the quantity Section 2.7 says the
bootstrap supports. The per-rung levels carry the known inflation from resampling inside
frozen folds; the gap is affected far less, because the inflation applies to both sides.

Values come from the gap files written by rung4_model.py, so the figure cannot drift
from the numbers the subsection quotes.

Usage:
    python scripts/fig_gap_forest.py
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
WINDOW = "20"

W_IN, H_IN = 4.72, 2.45      # 12 cm wide, printed 1:1 at 0.75\textwidth
INK, SOFT, EDGE = "#15202b", "#5b6b7d", "#c8d2dc"
RESOLVED, OPEN = "#b2182b", "#41618f"


def gaps() -> pd.DataFrame:
    """One row per framing, from the files rung4_model.py wrote."""
    ac = pd.read_csv(RES/f"19_rung4_gap_ci_all_w{WINDOW}.csv").set_index("framing")
    b = pd.read_csv(RES/f"19_rung4_gap_ci_w{WINDOW}.csv").set_index("quantity")
    rows = []
    for fr, metric in (("A", "$R^2$"), ("B", "$R^2$"), ("C", "ROC AUC")):
        if fr == "B":
            r = b.loc["gap"]
            lo, med, hi, n = r["lo2.5"], r["median"], r["hi97.5"], int(r.n_pairs)
            unit = "pairs"
        else:
            r = ac.loc[fr]
            lo, med, hi, n = r.gap_lo, r.gap_med, r.gap_hi, int(r.n_units)
            unit = "ligands" if fr == "A" else "pairs"
        rows.append(dict(framing=fr, metric=metric, lo=lo, med=med, hi=hi,
                         n=n, unit=unit, resolved=(lo > 0) or (hi < 0)))
    return pd.DataFrame(rows)


def main() -> None:
    g = gaps()
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(W_IN, H_IN))
    fig.subplots_adjust(left=0.295, right=0.985, top=0.97, bottom=0.21)

    ys = range(len(g))[::-1]
    ax.axvline(0.0, color=SOFT, lw=1.0, ls="--", zorder=1)
    for y, (_, r) in zip(ys, g.iterrows()):
        col = RESOLVED if r.resolved else OPEN
        ax.plot([r.lo, r.hi], [y, y], color=col, lw=2.0, solid_capstyle="round", zorder=3)
        for x in (r.lo, r.hi):
            ax.plot([x, x], [y-0.10, y+0.10], color=col, lw=1.4, zorder=3)
        ax.plot([r.med], [y], "o", ms=6.5, color=col, markeredgecolor="white",
                markeredgewidth=0.9, zorder=4)
        # above the interval rather than beside it, which ran the longest label off the axis
        ax.text(r.med, y + 0.16, f"{r.med:+.3f}   [{r.lo:+.3f}, {r.hi:+.3f}]",
                va="bottom", ha="center", fontsize=6.8, color=col)

    ax.set_yticks(list(ys))
    ax.set_yticklabels([f"framing {r.framing}\n{r.metric}, {r.n} {r.unit}"
                        for _, r in g.iterrows()], fontsize=7.4)
    ax.set_ylim(-0.55, len(g) - 0.25)
    ax.set_xlim(-0.70, 0.40)
    ax.set_xlabel("Rung 4 minus rung 1, with the 95% bootstrap interval", fontsize=7.8)
    ax.tick_params(axis="x", labelsize=7.2)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color(EDGE)
    ax.text(-0.68, -0.48, "rung 4 worse", fontsize=6.8, color=SOFT, ha="left")
    ax.text(0.38, -0.48, "rung 4 better", fontsize=6.8, color=SOFT, ha="right")

    FIGS.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_gap_forest.{ext}", dpi=400)
    print(g.to_string(index=False))
    print("\n  -> figures/fig_gap_forest.pdf / .png")


if __name__ == "__main__":
    main()
