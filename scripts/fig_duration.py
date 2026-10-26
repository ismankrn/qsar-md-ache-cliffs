#!/usr/bin/env python
"""Figure 5 — how much trajectory rung 4 actually needs.

The predictable objection to an MD-derived descriptor set is cost, and its honest form is a
real question: a descriptor that needs 50 ns per compound cannot be applied to a screening
library. That is answerable with a curve rather than an argument, because every run's JSON
carries the descriptors recomputed on the first 5, 10, 20 and 30 ns of analysed trajectory.

One asymmetry decides the layout. All 229 runs reach 20 ns of analysed trajectory, but only
the 55 runs from the original 50 ns design reach 30 ns, and those cover just 45 dev ligands
and 11 dev pairs. Eleven pairs cannot support framings B or C, so the 30 ns point appears for
framing A only, on its own 45-ligand subset, drawn as a separate series rather than as a
continuation of the full-data line. Joining them would compare a 217-ligand estimate with a
45-ligand one and read as a trend; the subset is also the earliest slice rather than a random
draw, so the two series are not interchangeable even at the windows they share.

Window lengths are *analysed* trajectory. The simulation each implies is EQUIL_NS (5 ns)
longer, which is the number to quote as a cost, and the upper axis shows it.

Usage:
    python scripts/fig_duration.py
"""
from __future__ import annotations
import sys, textwrap, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, roc_auc_score

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rung4_model import (REG, CLF, cv, pair_rows, sym_pairs, load_md, usable_descriptors,
                         SEED, TARGET)

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
WINDOWS = ["5", "10", "20", "30"]
EQUIL_NS = 5.0


def md_matrix(window: str, ids: list[str], cols: list[str] | None = None):
    """(matrix, ligands, columns) for the ligands in `ids` that reach this window.

    `cols=None` applies the replicate filter computed *at that window*, which is what the
    scoring pipeline does. Passing an explicit list instead holds the descriptor set fixed
    across windows, which is the only way to ask whether trajectory length matters on its own.
    """
    md = load_md(window)
    allc = [c for c in md.columns if c != "n_replicates" and md[c].dtype != object]
    keep = [c for c in (cols if cols is not None else (usable_descriptors(window) or []))
            if c in allc] or allc
    if "n_heavy" in allc and "n_heavy" not in keep:
        keep = keep + ["n_heavy"]
    have = [i for i in ids if i in md.index]
    return md.loc[have, keep].to_numpy(float), have, keep


def best(X, y, f, fams, metric, clf=False):
    out = None
    for mk in fams.values():
        try:
            v = metric(y, cv(X, y, f, mk, clf=clf))
        except Exception:
            continue
        out = v if out is None else max(out, v)
    return out


def main() -> None:
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv").merge(
        sl[["molecule_chembl_id", "split", "cv_fold"]], on="molecule_chembl_id")
    dev = lig[lig.split == "dev"].reset_index(drop=True)
    dev_ids = list(dev.molecule_chembl_id)
    dpair_all = pairs[pairs.split == "dev"].reset_index(drop=True)

    # The 45 ligands that reach 30 ns, scored at every window as their own series. Without
    # this the 30 ns point hangs alone and reads as a drop, when it is simply a different and
    # much smaller sample. A subset compared only with itself is the only honest way to ask
    # whether the curve is still rising past 20 ns.
    sub_ids = sorted(set(load_md("30").index) & set(dev_ids))

    # The descriptor set the replicate filter selects changes with the window -- 13, 13 and
    # 11 descriptors at 5, 10 and 20 ns, which with the heavy-atom covariate make 14, 14
    # and 12 features, and not the same ones. Scoring each window with its own set confounds two
    # things: how much trajectory the descriptors saw, and which descriptors were chosen. The
    # first draft of this figure reported a clean monotone rise and was titled accordingly;
    # holding the set fixed at the 20 ns selection flattens framing A entirely and makes
    # framing B peak at 10 ns. The rise was the filter, not the trajectory.
    FIXED = list(usable_descriptors("20") or [])

    rows = []
    for w in WINDOWS:
        for series, cols in (("adaptive", None), ("fixed", FIXED)):
            X, ids, keep = md_matrix(w, dev_ids, cols)
            sub = dev[dev.molecule_chembl_id.isin(ids)].set_index(
                "molecule_chembl_id").loc[ids]
            y, fold = sub.pchembl.to_numpy(), sub.cv_fold.to_numpy()
            dp = dpair_all[dpair_all.ligand_a.isin(ids) & dpair_all.ligand_b.isin(ids)
                           ].reset_index(drop=True)
            rows.append(dict(window_ns=float(w), framing="A", metric="r2",
                             value=best(X, y, fold, REG, r2_score), n_ligands=len(ids),
                             n_pairs=len(dp), n_features=len(keep), series=series))
            # Eleven pairs is not a sample; scoring B and C on it would produce numbers that
            # look like measurements and are not.
            if len(dp) >= 30:
                Xp, yp, fp, _ = pair_rows(X, y, None, ids, dp)
                rows.append(dict(window_ns=float(w), framing="B", metric="r2",
                                 value=best(Xp, yp, fp, REG, r2_score), n_ligands=len(ids),
                                 n_pairs=len(dp), n_features=len(keep), series=series))
                Xc, yc, fc = sym_pairs(X, ids, dp)
                rows.append(dict(window_ns=float(w), framing="C", metric="auc_roc",
                                 value=best(Xc, yc, fc, CLF, roc_auc_score, clf=True),
                                 n_ligands=len(ids), n_pairs=len(dp),
                                 n_features=len(keep), series=series))
            if series == "fixed":
                print(f"  window {w:>2s} ns: {len(ids):3d} ligands, {len(dp):3d} pairs, "
                      f"{len(keep):2d} fixed features"
                      + ("" if len(dp) >= 30 else "  -> framing A only"))
        # the 45 ligands that reach 30 ns, at every window, compared only with themselves
        Xs, sids, _sk = md_matrix(w, sub_ids, FIXED)
        ss = dev[dev.molecule_chembl_id.isin(sids)].set_index("molecule_chembl_id").loc[sids]
        rows.append(dict(window_ns=float(w), framing="A", metric="r2",
                         value=best(Xs, ss.pchembl.to_numpy(), ss.cv_fold.to_numpy(),
                                    REG, r2_score),
                         n_ligands=len(sids), n_pairs=0, n_features=len(_sk),
                         series="subset"))

    out = pd.DataFrame(rows)
    out.to_csv(RES/"16d_duration_curve.csv", index=False)
    print(f"\n  -> results/16d_duration_curve.csv")

    FIGS.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.2))
    titles = {"A": "Framing A — per-ligand potency",
              "B": "Framing B — signed pair difference",
              "C": "Framing C — cliff classification"}
    ylabs = {"A": "$R^2$", "B": "$R^2$", "C": "ROC AUC"}

    STYLE = {"adaptive": ("o-", "#2166ac", "per-window filter"),
             "fixed":    ("^-", "#1b7837", "filter fixed at 20 ns"),
             "subset":   ("s--", "#d6604d", None)}
    for ax, fr in zip(axes, "ABC"):
        s = out[out.framing == fr].sort_values("window_ns")
        full_n = int(s[s.series == "adaptive"].n_ligands.max())
        for series in ("adaptive", "fixed", "subset"):
            d = s[s.series == series]
            if series != "subset":
                # At 30 ns these series silently fall back to the 45 ligands that reach it,
                # so plotting them there would put a 45-ligand estimate on a 217-ligand line
                # and invite exactly the misreading the subset series exists to prevent.
                d = d[d.n_ligands == full_n]
            if not len(d):
                continue
            fmt, col, lab = STYLE[series]
            if series == "subset":
                lab = f"{int(d.n_ligands.iloc[0])}-ligand subset, fixed filter"
            ax.plot(d.window_ns, d.value, fmt, color=col, ms=5, lw=1.3, label=lab)
        ax.legend(fontsize=7.4, frameon=False, loc="best")
        ax.set_xticks([float(w) for w in WINDOWS])
        ax.set_xlabel("Analysed trajectory (ns)", fontsize=8.5)
        ax.set_ylabel(ylabs[fr], fontsize=8.5)
        # Above the secondary axis, not on it: at y=1.02 the letter sat on the leftmost
        # "Simulated per ligand" tick and panels (b) and (c) read as "(1b)" and "(1c)".
        ax.text(-0.17, 1.17, f"({'abc'['ABC'.index(fr)]})", transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="bottom", ha="left")
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        if fr == "C":
            ax.axhline(0.5, color="0.75", lw=0.8, ls=":")
        else:
            ax.axhline(0.0, color="0.75", lw=0.8)
        top = ax.secondary_xaxis("top")
        top.set_xticks([float(w) for w in WINDOWS])
        top.set_xticklabels([f"{float(w)+EQUIL_NS:.0f}" for w in WINDOWS], fontsize=7.4)
        top.set_xlabel("Simulated per ligand (ns)", fontsize=7.5, labelpad=2)
        top.tick_params(length=2)

    cap = ("Descriptors recomputed on the first 5/10/20/30 ns of analysed trajectory; frozen "
           "5-fold CV, best over algorithm families. The per-window filter selects 14/14/12 "
           "descriptors, so it confounds trajectory length with descriptor choice; the fixed "
           "series holds the 20 ns selection at every window and answers length alone. Only "
           "45 dev ligands (11 pairs) reach 30 ns, so the full-sample series stop at 20 ns.")
    fig.tight_layout(rect=(0, 0.11, 1, 0.91))
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig5_duration.{ext}", dpi=300)
    print(f"  -> figures/fig5_duration.pdf / .png")


if __name__ == "__main__":
    main()
