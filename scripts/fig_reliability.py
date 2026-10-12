#!/usr/bin/env python
"""Figure 4 — which MD-derived descriptors measure the ligand and which measure the seed.

A descriptor computed from one trajectory is only a property of the ligand if a second
trajectory of the same ligand, differing only in its initial velocities, reproduces it. The
replicate pilot measures that directly: sd_rep is the spread across replicates of one ligand,
sd_lig the spread across ligands, and their ratio says which dominates. The threshold of 2 was
fixed before the pilot ran.

The figure is organised by the three classes the descriptor set was designed around --
persistence, fluctuation, event count -- because that is where the result lies. Persistence
survives almost intact; fluctuation and event counts almost entirely do not. At 20 ns, 8 of
the 11 survivors are contact occupancies. That is actionable for anyone building descriptors
from trajectories of this length, and it is not specific to this target.

Reliability on the right axis is sd_lig^2 / (sd_lig^2 + 2 sd_rep^2), the fraction of a pair
difference's variance that is ligand rather than seed; ratio 2 corresponds to 0.67.

Usage:
    python scripts/fig_reliability.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
RATIO_USABLE = 2.0          # must match replicate_pilot.py
PRIMARY_WINDOW = 20

# The three design classes, assigned from how each quantity is computed rather than from what
# it is named: a fraction of frames is persistence, a spread or fluctuation amplitude is
# fluctuation, a count or a rate or a one-off flag is an event.
FAMILY = {
    "persistence": lambda n: n.startswith("occ_") or n in {
        "frac_lig_buried_mean", "n_bridge_water", "depth_mean"},
    "fluctuation": lambda n: n.endswith("_sd") or n.startswith("rmsf_") or n in {
        "rmsd_mean", "rmsd_max"},
    "event": lambda n: n.endswith("_per_ns") or n in {"escaped", "bridge_longest_ns"},
}
COLOUR = {"persistence": "#1b7837", "fluctuation": "#762a83", "event": "#d6604d",
          "other": "0.55"}


def family_of(name: str) -> str:
    for fam, test in FAMILY.items():
        if test(name):
            return fam
    return "other"


def main() -> None:
    d = pd.read_csv(RES/"18_replicate_pilot.csv")
    d["family"] = d.descriptor.map(family_of)
    d["reliability"] = d.sd_lig**2 / (d.sd_lig**2 + 2*d.sd_rep**2)

    other = sorted(d[d.family == "other"].descriptor.unique())
    if other:
        print(f"  unclassified (shown grey): {', '.join(other)}")

    prim = d[d.window == PRIMARY_WINDOW].set_index("descriptor")
    order = (prim.sort_values(["family", "ratio"], ascending=[True, True]).index.tolist())
    # descriptors absent at the primary window still belong on the axis
    order += [x for x in d.descriptor.unique() if x not in order]
    ypos = {name: i for i, name in enumerate(order)}

    FIGS.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 8.0),
                             gridspec_kw={"width_ratios": [2.4, 1.0]})

    # ---- (a) ratio per descriptor, one marker per window ----
    ax = axes[0]
    windows = sorted(d.window.unique())
    sizes = np.linspace(16, 46, len(windows))
    for w, s in zip(windows, sizes):
        sub = d[d.window == w]
        ax.scatter([max(r, 1e-3) for r in sub.ratio],
                   [ypos[n] for n in sub.descriptor],
                   s=s, facecolors="none",
                   edgecolors=[COLOUR[f] for f in sub.family], linewidths=0.9,
                   label=f"{w} ns")
    ax.axvline(RATIO_USABLE, color="0.3", lw=1.1, ls="--")
    ax.text(RATIO_USABLE, -0.8, f" threshold = {RATIO_USABLE:g}", fontsize=7.5,
            color="0.3", va="bottom", ha="left")
    # Descriptors whose ratio is undefined -- constant across every ligand, or non-zero in
    # too few -- plot nothing at all and would read as an omission. Mark them explicitly.
    undef = (d[~np.isfinite(d.ratio)].groupby("descriptor").verdict
             .agg(lambda s: s.mode().iat[0]))
    for name, verdict in undef.items():
        if name in ypos:
            ax.text(0.0016, ypos[name], f"{verdict} - ratio undefined", fontsize=6.5,
                    color="0.45", va="center", ha="left", style="italic")
    ax.set_xscale("log")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    for lab in ax.get_yticklabels():
        lab.set_color(COLOUR[family_of(lab.get_text())])
    ax.set_ylim(-1, len(order))
    ax.set_xlabel(r"$\mathrm{sd}_{\mathrm{ligand}} / \mathrm{sd}_{\mathrm{replicate}}$"
                  "   (log scale)", fontsize=9)
    ax.set_title("(a) Ligand signal against seed noise, by analysis window", fontsize=9.5)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="0.92", lw=0.6)
    ax.set_axisbelow(True)
    handles = [Line2D([], [], marker="o", ls="", mfc="none", mec="0.35",
                      ms=np.sqrt(s), label=f"{w} ns") for w, s in zip(windows, sizes)]
    handles += [Line2D([], [], marker="s", ls="", color=COLOUR[f], label=f)
                for f in ("persistence", "fluctuation", "event")]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="lower right", ncol=2)

    # ---- (b) verdict composition per window ----
    ax = axes[1]
    piv = (d.pivot_table(index="window", columns="verdict", values="descriptor",
                         aggfunc="count").fillna(0))
    cols = [c for c in ("usable", "seed-dominated", "sparse", "constant") if c in piv.columns]
    bottom = np.zeros(len(piv))
    shades = {"usable": "#1b7837", "seed-dominated": "#d6604d",
              "sparse": "0.70", "constant": "0.85"}
    for c in cols:
        ax.bar(range(len(piv)), piv[c].to_numpy(), bottom=bottom, width=0.62,
               color=shades[c], label=c)
        bottom += piv[c].to_numpy()
    ax.set_xticks(range(len(piv)))
    ax.set_xticklabels([f"{int(w)}" for w in piv.index], fontsize=8)
    ax.set_xlabel("Analysis window (ns)", fontsize=9)
    ax.set_ylabel("Descriptors", fontsize=9)
    ax.set_title("(b) How many survive, by window", fontsize=9.5)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    n_usable = int((d.verdict == "usable").sum())
    fig.suptitle("Of three descriptor design classes, only persistence survives a "
                 "replicate test", fontsize=11, y=0.985)
    fig.text(0.5, 0.008,
             f"{d.descriptor.nunique()} descriptors x {len(windows)} windows = {len(d)} "
             f"combinations; {n_usable} usable, {int((d.verdict=='seed-dominated').sum())} "
             f"seed-dominated. Threshold fixed before the pilot ran. Ratio {RATIO_USABLE:g} "
             f"corresponds to a reliability of "
             f"{RATIO_USABLE**2/(RATIO_USABLE**2+2):.2f}.",
             ha="center", fontsize=7.2, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig4_reliability.{ext}", dpi=300)
    print(f"  -> figures/fig4_reliability.pdf / .png")

    # per-family survival, the sentence the figure exists to support
    sv = (d.assign(ok=d.verdict.eq("usable"))
            .groupby("family").ok.agg(["sum", "count"]))
    print("\n  survival by design class (all windows):")
    for fam, r in sv.iterrows():
        print(f"    {fam:12s} {int(r['sum']):3d} / {int(r['count']):3d}")
    p20 = d[d.window == PRIMARY_WINDOW]
    occ = p20[(p20.verdict == "usable") & p20.descriptor.str.startswith("occ_")]
    print(f"\n  at {PRIMARY_WINDOW} ns: {int((p20.verdict=='usable').sum())} usable, "
          f"{len(occ)} of them contact occupancies")
    d.to_csv(RES/"18b_reliability_by_family.csv", index=False)
    print(f"  -> results/18b_reliability_by_family.csv")


if __name__ == "__main__":
    main()
