#!/usr/bin/env python
"""Figure — how the dataset is assembled, from ChEMBL records to 183 pairs.

Section 2.3 states the curation in four paragraphs and a reader has to hold nine
numbers in their head to follow it. This draws the same chain once, with the
exclusions beside the step that applies them, so the final set can be traced back to
the records it came from without re-reading the text.

Every number is recomputed here rather than copied from the prose. The record-level
chain is replayed from the cached ChEMBL download, the feasibility filter is re-run
over the whole pool, and the pair and ligand counts are read from the stage outputs.
A figure that restates hand-typed numbers would drift from the pipeline the first time
a threshold changed, and this one fails loudly instead, because `build_ligand_set`
and `chembl_screen` are imported rather than paraphrased.

The two tracks are drawn separately on purpose. Curation acts on compounds and the
cliff criterion acts on pairs, and the one structural decision in the study that is
easy to get wrong -- filtering the pool for simulation feasibility before pairing and
not after -- is only visible when the two tracks are side by side, since the arrow
into the pair track leaves the curated-compound box while the feasibility box feeds
the surviving-pairs box further down.

Usage:
    python scripts/fig_dataset_flow.py
"""
from __future__ import annotations
import re, sys, warnings
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from chembl_screen import SIM_CUTOFF, CLIFF_DELTA, SOFT_DELTA, REPLICATE_SD
from build_ligand_set import md_ok, CONTROL_MAX_DELTA, MW_MIN, MW_MAX, MAX_HEAVY, MAX_ROTB

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS, CACHE = ROOT/"results", ROOT/"figures", ROOT/"data"/"chembl_cache"
TARGET = "CHEMBL220"

BLUE, RED, GREEN, PURPLE, GREY = "#2166ac", "#b2182b", "#1b7837", "#762a83", "0.45"


def counts() -> dict:
    """Replay the curation chain and re-run the feasibility filter."""
    n = {}
    d = pd.read_csv(CACHE/f"{TARGET}.csv.gz", low_memory=False)
    n["records"] = len(d)
    d = d[(d.standard_relation == "=") & (d.standard_units == "nM")
          & d.pchembl_value.notna() & d.data_validity_comment.isna()]
    n["valid"] = len(d)
    d = d[d.potential_duplicate.fillna(0).astype(int) == 0]
    n["not_duplicate"] = len(d)
    best = d.standard_type.value_counts().idxmax()
    n["assay_type"] = best
    d = d[d.standard_type == best]
    n["single_assay"] = len(d)

    pool = pd.read_csv(RES/f"01_compounds_{TARGET}.csv")
    n["compounds"] = len(pool)
    n["p_lo"], n["p_hi"] = float(pool.pchembl.min()), float(pool.pchembl.max())
    # the aggregation step is the only one whose input count is not on disk; the
    # discarded structures are, so the pre-filter structure count follows from it
    feas = pd.read_csv(RES/"01_feasibility.csv").set_index("target_chembl_id")
    n["similar"] = int(feas.loc[TARGET, "cliff_pairs_soft"])
    n["constitutional"] = int(feas.loc[TARGET, "cliffs_constitutional"])
    n["stereochemical"] = int(feas.loc[TARGET, "cliffs_stereochemical"])
    n["strict"] = len(pd.read_csv(RES/f"01_cliffs_{TARGET}.csv"))

    chk = pool.smiles_std.map(md_ok)
    ok = [c[0] for c in chk]
    n["feasible"] = int(sum(ok))
    n["infeasible"] = len(pool) - n["feasible"]
    why = pd.Series([re.sub(r"[0-9]+", "N", c[1]) for c, k in zip(chk, ok) if not k])
    n["why"] = why.value_counts().to_dict()

    pairs = pd.read_csv(RES/f"02_pairs_{TARGET}.csv")
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
    n["cliff_pairs"] = int((pairs.group == "cliff").sum())
    n["control_pairs"] = int((pairs.group == "control").sum())
    n["cliff_ligands"] = int((lig.group == "cliff").sum())
    n["control_ligands"] = int((lig.group == "control").sum())
    n["pairs"] = len(pairs)
    n["ligands"] = len(lig)
    return n


def box(ax, x, y, w, h, text, color, *, lw=1.3, fill="white", fs=7.6, weight=None):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.008",
                                linewidth=lw, edgecolor=color, facecolor=fill, zorder=2))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color="0.12",
            zorder=3, linespacing=1.45, fontweight=weight)


def arrow(ax, a, b, color, *, lw=1.2, ls="-", conn="arc3,rad=0"):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=11,
                                 connectionstyle=conn, linewidth=lw, linestyle=ls,
                                 color=color, shrinkA=0, shrinkB=0, zorder=1))


def main() -> None:
    n = counts()
    pd.DataFrame([{k: v for k, v in n.items() if k != "why"}]).T.rename(
        columns={0: "value"}).to_csv(RES/"30_dataset_flow.csv")
    print("  -> results/30_dataset_flow.csv")
    for k, v in n.items():
        print(f"    {k:16s} {v}")

    FIGS.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(11.8, 6.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    # Exclusions sit outside both tracks so the channel between the tracks stays
    # clear for the two cross-track arrows, which are the part of the layout that
    # carries an argument rather than a count.
    EX, EW = 0.068, 0.122          # exclusions, left of the compound track
    LX, LW = 0.335, 0.262          # compound track
    RX, RW = 0.668, 0.275          # pair track
    PX, PW = 0.902, 0.128          # exclusions, right of the pair track

    ax.text(LX, 0.985, "Compound curation", ha="center", fontsize=9.5,
            fontweight="bold", color=BLUE)
    ax.text(RX, 0.985, "Pair assembly", ha="center", fontsize=9.5,
            fontweight="bold", color=RED)

    left = [
        (0.905, 0.100, f"{n['records']:,} activity records\nhuman AChE ({TARGET})\n"
                       f"binding assays, nM, pChEMBL present", None),
        (0.725, 0.085, f"{n['single_assay']:,} {n['assay_type']} records\n"
                       f"one assay type only",
         f"{n['records']-n['not_duplicate']:,} records flagged as\npotential duplicates\n\n"
         f"{n['not_duplicate']-n['single_assay']:,} records of another\nassay type"),
        (0.545, 0.105, f"{n['compounds']:,} curated compounds\n"
                       f"median over replicates\n"
                       f"pChEMBL {n['p_lo']:.2f} to {n['p_hi']:.2f}",
         f"structures whose replicates\ndisagree by more than\n"
         f"{REPLICATE_SD:g} log unit"),
        (0.355, 0.070, f"{n['feasible']:,} simulation-feasible compounds", 
         f"{n['infeasible']} excluded on elements,\nmass {MW_MIN} to {MW_MAX} Da,\n"
         f"$\\leq${MAX_HEAVY} heavy atoms,\n$\\leq${MAX_ROTB} rotatable bonds,\n"
         f"3D embedding"),
    ]
    for i, (y, h, txt, exc) in enumerate(left):
        box(ax, LX, y, LW, h, txt, BLUE)
        if i:
            arrow(ax, (LX, left[i-1][0] - left[i-1][1]/2), (LX, y + h/2), BLUE)
        if exc:
            ym = (left[i-1][0] + y)/2
            box(ax, EX, ym, EW, 0.105, exc, GREY, lw=0.9, fill="0.965", fs=6.5)
            arrow(ax, (LX - LW/2, ym), (EX + EW/2, ym), GREY, lw=0.9, ls=(0, (3, 2)))

    right = [
        (0.905, 0.105, RED, f"{n['similar']} near-identical pairs\n"
                            f"Tanimoto $\\geq$ {SIM_CUTOFF:.2f} and\n"
                            f"$|\\Delta$pChEMBL$| \\geq$ {SOFT_DELTA:.1f}", None),
        (0.725, 0.105, RED, f"{n['strict']} activity cliffs\n"
                            f"$|\\Delta$pChEMBL$| \\geq$ {CLIFF_DELTA:.1f}, "
                            f"a hundredfold\ndifference in potency",
         f"{n['stereochemical']} stereochemical pairs,\nwhich no stereo-blind\n"
         f"fingerprint can resolve"),
        (0.575, 0.055, RED, f"{n['constitutional']} constitutional cliff pairs", None),
        (0.420, 0.080, RED, f"{n['cliff_pairs']} cliff pairs\n"
                            f"{n['cliff_ligands']} compounds, both members feasible", None),
        (0.230, 0.115, GREEN, f"{n['control_pairs']} matched control pairs\n"
                              f"Tanimoto $\\geq$ {SIM_CUTOFF:.2f} and "
                              f"$|\\Delta$pChEMBL$| \\leq$ {CONTROL_MAX_DELTA:.2f}\n"
                              f"{n['control_ligands']} compounds, matched on similarity\n"
                              f"and on mean potency", None),
    ]
    for i, item in enumerate(right):
        y, h, col, txt = item[:4]
        exc = item[4]
        box(ax, RX, y, RW, h, txt, col)
        if i:
            py, ph = right[i-1][0], right[i-1][1]
            # the controls are matched to the cliff pairs rather than derived from
            # them, so that one link is dashed and labelled; a solid arrow would
            # read as 97 cliff pairs turning into 86 control pairs
            if col == GREEN:
                arrow(ax, (RX, py - ph/2), (RX, y + h/2), GREEN, ls=(0, (4, 2)))
                ax.text(RX + 0.012, (py - ph/2 + y + h/2)/2,
                        "  matched to the\n  cliff pairs", fontsize=6.4, color=GREEN,
                        ha="left", va="center", linespacing=1.35)
            else:
                arrow(ax, (RX, py - ph/2), (RX, y + h/2), RED)
        if exc:
            ym = (right[i-1][0] + y)/2
            box(ax, PX, ym, PW, 0.085, exc, GREY, lw=0.9, fill="0.965", fs=6.5)
            arrow(ax, (RX + RW/2, ym), (PX - PW/2, ym), GREY, lw=0.9, ls=(0, (3, 2)))

    # The two cross-track arrows say where each track draws from. Pairs are
    # enumerated over the curated pool, while only the feasibility-filtered pool
    # decides which pairs survive, and that order is the point.
    arrow(ax, (LX + LW/2, 0.545), (RX - RW/2, 0.905), BLUE, lw=1.1,
          conn="arc3,rad=-0.10")
    arrow(ax, (LX + LW/2, 0.355), (RX - RW/2, 0.420), BLUE, lw=1.1)

    box(ax, 0.5, 0.058, 0.66, 0.062,
        f"Final set      {n['ligands']} ligands,      {n['pairs']} pairs      "
        f"({n['cliff_pairs']} cliff, {n['control_pairs']} control)",
        PURPLE, lw=1.6, fill="#f4eef6", fs=8.6, weight="bold")
    arrow(ax, (RX, 0.230 - 0.115/2), (RX, 0.058 + 0.031), PURPLE)
    arrow(ax, (RX - RW/2, 0.400), (0.515, 0.058 + 0.031), PURPLE, lw=1.1,
          conn="angle,angleA=180,angleB=90,rad=6")

    fig.tight_layout(pad=0.3)
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_dataset_flow.{ext}", dpi=300)
    print("  -> figures/fig_dataset_flow.pdf / .png")


if __name__ == "__main__":
    main()
