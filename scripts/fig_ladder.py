#!/usr/bin/env python
"""Figure 2 — the descriptor ladder across all three framings.

One point per rung per framing: the best score over that rung's descriptor arms and over the
six (five, for classification) algorithm families, on the frozen folds. Bootstrap intervals
are computed the same way for every rung, which is the point of the figure -- a ladder drawn
from point estimates alone invites "is that difference real?" and cannot answer it.

Three things this is careful about.

Every rung is scored on the same 217 dev ligands and the same 145 dev pairs. That is not an
assumption: coverage was checked and all five rungs are complete over the dev split, so no
rung is flattered by an easier subset. Stage 19b found the opposite situation while rung 4
was still running, and the correction mattered.

Rung 1 carries both of the feature sets stage 05 benchmarked, and its score is the best over
both. Which one wins depends on the framing -- binary ECFP4 takes framing B at +0.200 where
count ECFP4 + RDKit2D manages -0.007, and the order reverses on framing C. Scoring one set
only would hand rung 1 a figure well below what it achieves and inflate every later rung by
the difference.

Rung 2 appears twice, as the single rigid 4EY7 receptor and as the median over eight receptor
conformations, because stage 15 showed the ensemble median reaches +0.274 on framing B -- the
highest structure-based figure in the study, and not significantly above the rigid receptor.
Collapsing them to one point would hide the best result the structural ladder produces.

Usage:
    python scripts/fig_ladder.py                 # layout check, no intervals
    python scripts/fig_ladder.py --n-boot 300    # the figure for the manuscript
"""
from __future__ import annotations
import argparse, glob, json, os, re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, roc_auc_score

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rung4_model import (REG, CLF, cv, pair_rows, sym_pairs, boot_gap, load_md,
                         usable_descriptors, SEED, TARGET, META)

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"

# label -> (file, columns to drop); None means "built in code"
RUNGS = [
    ("Rung 1\n2D",        "r1"),
    ("Rung 2\npose",      "r2"),
    ("Rung 2\nensemble",  "r2e"),
    ("Rung 3\nMM/GBSA",   "r3"),
    ("Rung 4\nMD",        "r4"),
]
DROP = {"molecule_chembl_id", "n_seeds_used", "seed_used"}

# Panel identity only. The descriptive text lives in the caption, and draw() must not
# depend on anything main() computes, so that --from-csv can redraw without rescoring.
PANELS = [("A", "$R^2$"), ("B", "$R^2$"), ("C", "ROC AUC")]


def numeric_matrix(df: pd.DataFrame, ids: list[str]) -> np.ndarray:
    """Rows in `ids` order, numeric columns only, one row per ligand."""
    df = df.drop_duplicates("molecule_chembl_id").set_index("molecule_chembl_id")
    cols = [c for c in df.columns if c not in DROP and df[c].dtype != object]
    return df.loc[ids, cols].to_numpy(float)


def rung2_matrix(ids: list[str]) -> np.ndarray:
    """Rung 2 exactly as rung2_model.py builds it: descriptors merged with the bridging
    features, minus the per-descriptor spread columns and the seed count."""
    d2 = pd.read_csv(RES/"08_rung2_descriptors.csv")
    br = pd.read_csv(RES/"08_rung2_bridging.csv")
    desc = d2.merge(br, on="molecule_chembl_id", how="left")
    desc = desc.drop(columns=[c for c in desc.columns
                              if c.startswith("spread_") or c == "n_seeds_used"],
                     errors="ignore")
    return numeric_matrix(desc, ids)


def rung3_matrix(ids: list[str]) -> np.ndarray:
    """Rung 3 exactly as rung3_analyse.py builds it: the energy and per-residue terms only.
    n_interface_atoms and n_lig_atoms are deliberately absent -- they are size proxies."""
    d = pd.read_csv(RES/"09_rung3_mmgbsa.csv").drop_duplicates("molecule_chembl_id")
    feat = [c for c in d.columns if c.startswith(("dE_", "dG_", "res_", "dSASA"))]
    return d.set_index("molecule_chembl_id").loc[ids, feat].to_numpy(float)


def build_arms(ids: list[str], lig: pd.DataFrame, window: str):
    """{rung key: [(label, per-ligand matrix), ...]}  — all aligned to `ids`."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    from rdkit.ML.Descriptors import MoleculeDescriptors
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    names = [d[0] for d in Descriptors._descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
    smi = dict(zip(lig.molecule_chembl_id, lig.smiles_std))
    mols = [Chem.MolFromSmiles(smi[c]) for c in ids]
    cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], float)
    rdk = np.nan_to_num(np.array([calc.CalcDescriptors(m) for m in mols], float),
                        nan=0.0, posinf=0.0, neginf=0.0)
    arms = {
        "r1": [("binary ECFP4", np.array([gen.GetFingerprintAsNumPy(m) for m in mols], float)),
               ("count ECFP4 + RDKit2D", np.hstack([cnt, rdk]))],
        # Each rung uses the feature definition its own stage published, not "every numeric
        # column", so the figure agrees with the text. Getting this wrong the first time moved
        # framing A for rung 2 from +0.046 to +0.096 and framing C for rung 3 to 0.800, the
        # latter because the raw MM/GBSA table carries n_lig_atoms and n_interface_atoms --
        # size proxies that stage 10 already showed separate cliffs from controls for the
        # wrong reason, and which rung3_analyse.py therefore excludes.
        "r2": [("4EY7 + bridging", rung2_matrix(ids))],
        "r2e": [("ensemble median",
                 numeric_matrix(pd.read_csv(RES/"15_rung2_ensemble_median.csv"), ids))],
        "r3": [("MM/GBSA", rung3_matrix(ids))],
    }
    md = load_md(window)
    keep = usable_descriptors(window)
    allc = [c for c in md.columns if c not in ("n_replicates",) and md[c].dtype != object]
    keep = [c for c in (keep or []) if c in allc]
    if "n_heavy" in allc and "n_heavy" not in keep:
        keep = keep + ["n_heavy"]
    arms["r4"] = [(f"replicate-filtered ({len(keep)})", md.loc[ids, keep].to_numpy(float)),
                  (f"all descriptors ({len(allc)})", md.loc[ids, allc].to_numpy(float))]
    return arms


def best_expanded(arms, fams, metric, clf=False):
    """Best (value, arm label, algorithm, n_features) over pre-expanded arms."""
    best = None
    for label, (X, yy, ff) in arms:
        for aname, mk in fams.items():
            try:
                v = metric(yy, cv(X, yy, ff, mk, clf=clf))
            except Exception:
                continue
            if best is None or v > best[0]:
                best = (v, label, aname, X.shape[1])
    return best


def draw(out, a) -> None:
        # ---------------- figure ----------------
        FIGS.mkdir(exist_ok=True)
        fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.3))
        xs = np.arange(len(RUNGS))
        for ax, (fr, ylab) in zip(axes, PANELS):
            s = out[out.framing == fr].set_index("key").reindex([k for _, k in RUNGS])
            v = s["value"].to_numpy(float)
            ax.axhline(0.0, color="0.75", lw=0.8, zorder=1)
            # rung 1 is the thing every other rung has to beat, so draw it as a line across
            ax.axhline(v[0], color="#b2182b", lw=1.0, ls="--", zorder=1,
                       label="Rung 1 (2D) level")
            # Draw the bars whenever the data carries them. Keying this off the --n-boot
            # flag meant --from-csv silently produced a figure with no intervals, which
            # looked finished and was not.
            lo = s["ci_lo"].to_numpy(float)
            hi = s["ci_hi"].to_numpy(float)
            if np.isfinite(lo).all() and np.isfinite(hi).all():
                err = np.vstack([v - lo, hi - v])
                ax.errorbar(xs, v, yerr=err, fmt="none", ecolor="0.35", lw=1.2, capsize=3,
                            zorder=2)
            ax.plot(xs, v, "o-", color="#2166ac", ms=6, lw=1.2, zorder=3)
            ax.set_xticks(xs)
            ax.set_xticklabels([lab for lab, _ in RUNGS], fontsize=7.5)
            ax.text(0.0, 1.02, f"({'abc'['ABC'.index(fr)]})", transform=ax.transAxes,
                    fontsize=11, fontweight="bold", va="bottom", ha="left")
            ax.set_ylabel(ylab, fontsize=9)
            ax.tick_params(labelsize=8)
            ax.spines[["top", "right"]].set_visible(False)
            if fr == "C":
                ax.axhline(0.5, color="0.75", lw=0.8, ls=":", zorder=1)
                ax.set_ylim(0.45, 0.92)
                ax.text(0.02, 0.5, "chance", transform=ax.get_yaxis_transform(),
                        va="bottom", fontsize=7.4, color="0.5")
        axes[0].legend(fontsize=7.5, frameon=False, loc="lower left")
        # The first draft said "does not improve any endpoint". That is false as written: on
        # framing B the ensemble rung 2 reaches +0.274 against rung 1's +0.200. Stage 15 showed
        # the difference is not resolvable, which is the claim the data supports.
        fig.tight_layout(rect=(0, 0.03, 1, 0.94))
        for ext in ("pdf", "png"):
            fig.savefig(FIGS/f"fig2_ladder.{ext}", dpi=300)
        print(f"  -> figures/fig2_ladder.pdf / .png")




def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", default="20", help="rung-4 analysis window in ns")
    ap.add_argument("--n-boot", type=int, default=0,
                    help="bootstrap replicates per rung; 0 draws the figure without intervals")
    ap.add_argument("--from-csv", action="store_true",
                    help="redraw from results/20_ladder.csv without recomputing. The "
                         "300-replicate bootstrap takes eleven hours, so a change to the "
                         "figure's appearance must not require repeating it.")
    a = ap.parse_args()

    if a.from_csv:
        out = pd.read_csv(RES/"20_ladder.csv")
        draw(out, a)
        return

    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv").merge(
        sl[["molecule_chembl_id", "split", "cv_fold"]], on="molecule_chembl_id")
    dev = lig[lig.split == "dev"].reset_index(drop=True)
    ids = list(dev.molecule_chembl_id)
    y = dev.pchembl.to_numpy()
    fold = dev.cv_fold.to_numpy()
    dpair = pairs[pairs.split == "dev"].reset_index(drop=True)
    dpair = dpair[dpair.ligand_a.isin(ids) & dpair.ligand_b.isin(ids)].reset_index(drop=True)

    print("=" * 74)
    print(f"LADDER — {len(ids)} dev ligands, {len(dpair)} dev pairs, rung-4 window {a.window} ns")
    print("=" * 74)
    arms = build_arms(ids, lig, a.window)

    # framing -> (make_xy, families, metric, clf, units)
    def xy_A(M):
        return M, y, fold

    def xy_B(M):
        Xp, yp, fp, _ = pair_rows(M, y, None, ids, dpair)
        return Xp, yp, fp

    def xy_C(M):
        return sym_pairs(M, ids, dpair)

    nB = len(xy_B(arms["r1"][0][1])[1]) // 2
    nC = len(xy_C(arms["r1"][0][1])[1])
    FRAMINGS = [
        ("A", "Framing A\nper-ligand potency", "$R^2$", xy_A, REG, r2_score, False,
         [np.array([i]) for i in range(len(ids))]),
        ("B", "Framing B\nsigned pair difference", "$R^2$", xy_B, REG, r2_score, False,
         [np.array([2*i, 2*i+1]) for i in range(nB)]),
        ("C", "Framing C\ncliff classification", "ROC AUC", xy_C, CLF, roc_auc_score, True,
         [np.array([i]) for i in range(nC)]),
    ]

    rows = []
    for fr, _title, _ylab, make_xy, fams, metric, clf, units in FRAMINGS:
        print(f"\n  --- framing {fr} ---")
        # Expand every arm once. Framings B and C walk the pair table to build their rows, and
        # doing that inside the bootstrap loop would repeat 145 iterrows calls for every arm of
        # every replicate -- the resampling only ever needs to index rows of a fixed matrix.
        exp = {key: [(lab, make_xy(M)) for lab, M in arms[key]] for _, key in RUNGS}
        for label, key in RUNGS:
            b = best_expanded(exp[key], fams, metric, clf=clf)
            if b is None:
                continue
            val, arm, alg, nfeat = b
            lo = hi = np.nan
            if a.n_boot:
                mk = lambda XYF: (lambda r, XYF=XYF: (XYF[0][r], XYF[1][r], XYF[2][r]))
                q, _, _ = boot_gap(units, [(mk(XYF), fams) for _, XYF in exp[key]],
                                   None, metric, a.n_boot, SEED, clf=clf)
                lo, hi = q[0], q[2]
            rows.append(dict(framing=fr, rung=label.replace("\n", " "), key=key, value=val,
                             arm=arm, algorithm=alg, n_features=nfeat, ci_lo=lo, ci_hi=hi,
                             n_ligands=len(ids), n_pairs=len(dpair)))
            ci = f"  [{lo:+.3f}, {hi:+.3f}]" if a.n_boot else ""
            print(f"    {label.replace(chr(10),' '):28s} {val:+.3f}  ({alg}, {arm}){ci}")

    out = pd.DataFrame(rows)
    out.to_csv(RES/"20_ladder.csv", index=False)
    print(f"\n  -> results/20_ladder.csv")

    draw(out, a)


if __name__ == "__main__":
    main()
