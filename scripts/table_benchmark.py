#!/usr/bin/env python
"""Table 2 — every rung, every framing, every algorithm.

Figure 2 shows one number per rung: the best over its descriptor arms and over the algorithm
families. That is the right summary and it hides something a reader is entitled to see, which
is how much of each rung's score is the descriptors and how much is the choice of model. The
full grid is 5 rungs x up to 2 arms x 3 framings x 6 (or 5) algorithms, which belongs in the
supplement; the main text gets the best per rung per framing with the winning algorithm named.

Descriptor arms are imported from fig_ladder.py rather than rebuilt, so the two cannot drift.
Each rung uses the feature definition its own stage published -- rung 2 with its bridging
features merged in, rung 3 restricted to the energy and per-residue terms so the size proxies
in the raw MM/GBSA table stay out, and rung 1 carrying both benchmarked feature sets because
the winner depends on the framing.

Outputs:
    results/22_benchmark_full.csv   every cell
    results/22_benchmark_main.tex   best per rung per framing, for the main text
    results/22_benchmark_si.tex     the full grid, for the supplement

Usage:
    python scripts/table_benchmark.py
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import (r2_score, mean_absolute_error, roc_auc_score,
                             average_precision_score, matthews_corrcoef,
                             balanced_accuracy_score)
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rung4_model import REG, CLF, cv, pair_rows, sym_pairs, SEED, TARGET
from fig_ladder import build_arms, RUNGS

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
LABEL = {k: lab.replace("\n", " ") for lab, k in RUNGS}


def esc(s: str) -> str:
    return str(s).replace("_", r"\_").replace("&", r"\&")


def main() -> None:
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv").merge(
        sl[["molecule_chembl_id", "split", "cv_fold"]], on="molecule_chembl_id")
    dev = lig[lig.split == "dev"].reset_index(drop=True)
    ids = list(dev.molecule_chembl_id)
    y, fold = dev.pchembl.to_numpy(), dev.cv_fold.to_numpy()
    dpair = pairs[(pairs.split == "dev") & pairs.ligand_a.isin(ids)
                  & pairs.ligand_b.isin(ids)].reset_index(drop=True)
    print(f"  {len(ids)} dev ligands, {len(dpair)} dev pairs")

    arms = build_arms(ids, lig, "20")
    rows = []
    for lab, key in RUNGS:
        for arm_name, M in arms[key]:
            # framing A
            for alg, mk in REG.items():
                p = cv(M, y, fold, mk)
                rows.append(dict(rung=LABEL[key], arm=arm_name, framing="A", algorithm=alg,
                                 metric="r2", value=r2_score(y, p),
                                 rho=spearmanr(y, p).statistic,
                                 mae=mean_absolute_error(y, p), n_features=M.shape[1]))
            # framing B
            Xp, yp, fp, _ = pair_rows(M, y, None, ids, dpair)
            for alg, mk in REG.items():
                p = cv(Xp, yp, fp, mk)
                rows.append(dict(rung=LABEL[key], arm=arm_name, framing="B", algorithm=alg,
                                 metric="r2", value=r2_score(yp, p),
                                 rho=spearmanr(yp, p).statistic,
                                 mae=mean_absolute_error(yp, p), n_features=Xp.shape[1]))
            # framing C
            Xc, yc, fc = sym_pairs(M, ids, dpair)
            for alg, mk in CLF.items():
                pr = cv(Xc, yc, fc, mk, clf=True)
                hard = (pr >= 0.5).astype(int)
                rows.append(dict(rung=LABEL[key], arm=arm_name, framing="C", algorithm=alg,
                                 metric="auc_roc", value=roc_auc_score(yc, pr),
                                 auc_pr=average_precision_score(yc, pr),
                                 mcc=matthews_corrcoef(yc, hard),
                                 bal_acc=balanced_accuracy_score(yc, hard),
                                 n_features=Xc.shape[1]))
        print(f"    {LABEL[key]:22s} done")

    out = pd.DataFrame(rows)
    out.to_csv(RES/"22_benchmark_full.csv", index=False)
    print(f"\n  -> results/22_benchmark_full.csv  ({len(out)} cells)")

    # ---- main text: best per rung per framing ----
    best = (out.sort_values("value", ascending=False)
               .groupby(["rung", "framing"], as_index=False).first())
    order = [LABEL[k] for _, k in RUNGS]
    tex = [r"\begin{tabular}{lccc}", r"\toprule",
           r"& \framingA{} & \framingB{} & \framingC{} \\",
           r"Rung & $R^2$ & $R^2$ & ROC AUC \\", r"\midrule"]
    for r in order:
        cells = []
        for fr in "ABC":
            s = best[(best.rung == r) & (best.framing == fr)]
            if s.empty:
                cells.append("--")
                continue
            v, alg = float(s.value.iloc[0]), s.algorithm.iloc[0]
            num = f"{v:+.3f}" if fr != "C" else f"{v:.3f}"
            # Braces around \scriptsize, not a bare \\ -- the first version emitted
            # "\\scriptsize", which LaTeX reads as a line break followed by the word.
            cells.append(num + r" {\scriptsize (" + esc(alg) + ")}")
        tex.append(f"{esc(r)} & " + " & ".join(cells) + r" \\")
    tex += [r"\bottomrule", r"\end{tabular}"]
    (RES/"22_benchmark_main.tex").write_text("\n".join(tex) + "\n")
    print(f"  -> results/22_benchmark_main.tex")

    # ---- supplement: two grids, because the algorithm families differ ----
    # Framings A and B are regressions and C is a classification, so they do not share a
    # model suite. Forcing them into one grid printed "--" for Ridge and PLS on framing C and
    # implied those models had been tried and failed, when they were never applicable.
    def grid(path, framings, fams, headers, fmts):
        lines = [r"\begin{tabular}{lll" + "r"*len(framings) + "}", r"\toprule",
                 r"Rung & Descriptor arm & Algorithm & " + " & ".join(headers) + r" \\",
                 r"\midrule"]
        for r in order:
            sub = out[out.rung == r]
            for arm in sub.arm.unique():
                s = sub[sub.arm == arm]
                for i_a, alg in enumerate([a for a in fams if a in set(s.algorithm)]):
                    lead = (f"{esc(r)} & {esc(arm)}" if i_a == 0 else " & ")
                    cells = []
                    for fr, fmt in zip(framings, fmts):
                        q = s[(s.framing == fr) & (s.algorithm == alg)]
                        cells.append(fmt.format(float(q.value.iloc[0])) if len(q) else "--")
                    lines.append(f"{lead} & {esc(alg)} & " + " & ".join(cells) + r" \\")
                lines.append(r"\addlinespace")
        lines += [r"\bottomrule", r"\end{tabular}"]
        (RES/path).write_text("\n".join(lines) + "\n")
        print(f"  -> results/{path}")

    grid("22_benchmark_si_reg.tex", ["A", "B"], list(REG),
         [r"\framingA{} $R^2$", r"\framingB{} $R^2$"], ["{:+.3f}", "{:+.3f}"])
    grid("22_benchmark_si_clf.tex", ["C"], list(CLF),
         [r"\framingC{} ROC AUC"], ["{:.3f}"])

    print("\n  best per rung per framing:")
    piv = best.pivot(index="rung", columns="framing", values="value").reindex(order)
    print(piv.round(3).to_string())


if __name__ == "__main__":
    main()
