#!/usr/bin/env python
"""Table 3 — how much a partial subset inflates an expensive descriptor set's apparent gain.

An expensive descriptor set arrives incrementally. The tempting comparison, while it is still
arriving, is the subset you have against the baseline's published figure over the full set.
This study ran that comparison for real and it was wrong by more than the effect it reported:
at 59 of 145 pairs rung 4 read +0.478 against rung 1's published +0.200, a gain of +0.278 and
a headline. Recomputing rung 1 on exactly those pairs removed all of it.

The table reconstructs that, and separates two explanations for it.

Completion order is real, not simulated. Each run's JSON carries the mtime it was written on
the machine that produced it, preserved through rsync -a, so the prefixes below are the
ligands this project actually held at each point in the campaign. The alternative -- ordering
by size and calling it a proxy -- would assume the mechanism the table is meant to test.

Random subsets of the same size are the control. If the inflation were only "small samples
are noisy", random subsets would show it too. They do not, so what matters is *which* ligands
arrive first, not how many.

What makes them different is not one thing, and the table records all three rather than
asserting the convenient one. The first 111 ligands to finish carry a median of 31 heavy atoms
against 35 for the rest (Mann-Whitney p = 0.0005), they hold every multi-replicate run in the
study, and they average 28.8 ns of analysed trajectory against 20.0 -- because the campaign
began on a 50 ns design. Smaller systems run faster, replicate averaging lowers descriptor
noise, and longer runs were the early protocol; all three travel together and n = 217 cannot
separate them. The claim this table supports is that the early subset was systematically
different and that comparing it against a published full-set baseline manufactured a gain,
not that ligand size alone caused it.

Outputs results/21_subset_inflation.csv and a LaTeX fragment for the manuscript.

Usage:
    python scripts/table_subset_inflation.py --n-random 40
"""
from __future__ import annotations
import argparse, glob, os, re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import r2_score

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rung4_model import REG, cv, pair_rows, load_md, usable_descriptors, SEED, TARGET

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
PUBLISHED_RUNG1_B = 0.200      # 05_rung1_benchmark.csv, framing B, binary ECFP4 + GB, 145 pairs
TARGET_PAIRS = [60, 85, 120, 145]


def completion_order(dev_ids: set[str]) -> list[str]:
    """Dev ligands in the order their first run finished, by JSON mtime."""
    seen: dict[str, float] = {}
    for f in glob.glob(str(RES/"16_rung4_per_run"/"*.json")):
        lig = re.sub(r"_r\d+$", "", os.path.basename(f)[:-5])
        if lig not in dev_ids:
            continue
        t = os.path.getmtime(f)
        if lig not in seen or t < seen[lig]:
            seen[lig] = t
    return [k for k, _ in sorted(seen.items(), key=lambda kv: kv[1])]


def best_r2(X, y, f) -> float:
    vals = []
    for mk in REG.values():
        try:
            vals.append(r2_score(y, cv(X, y, f, mk)))
        except Exception:
            pass
    return max(vals) if vals else float("nan")


def score_subset(ids, dev, dpair, md_keep_cols, md, X1s, dev_index):
    """(rung4 R2, matched rung1 R2, n_pairs) on framing B for this ligand subset."""
    dp = dpair[dpair.ligand_a.isin(ids) & dpair.ligand_b.isin(ids)].reset_index(drop=True)
    if len(dp) < 20:
        return np.nan, np.nan, len(dp)
    sub = dev.set_index("molecule_chembl_id").loc[ids]
    y, _f = sub.pchembl.to_numpy(), sub.cv_fold.to_numpy()
    ids = list(ids)
    X4 = md.loc[ids, md_keep_cols].to_numpy(float)
    Xp, yp, fp, _ = pair_rows(X4, y, None, ids, dp)
    r4 = best_r2(Xp, yp, fp)
    # rung 1 takes the best over both published feature sets, as stage 05 reported. Scoring
    # one set only understates it by up to 0.21 on this framing and would inflate rung 4 again
    # through a different door.
    rowsel = [dev_index[i] for i in ids]
    r1 = max(best_r2(*pair_rows(X[rowsel], y, None, ids, dp)[:3]) for X in X1s)
    return r4, r1, len(dp)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-random", type=int, default=40,
                    help="random subsets drawn per size for the control")
    ap.add_argument("--window", default="20")
    a = ap.parse_args()

    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv").merge(
        sl[["molecule_chembl_id", "split", "cv_fold"]], on="molecule_chembl_id")
    dev = lig[lig.split == "dev"].reset_index(drop=True)
    all_ids = list(dev.molecule_chembl_id)
    dpair = pairs[(pairs.split == "dev") & pairs.ligand_a.isin(all_ids)
                  & pairs.ligand_b.isin(all_ids)].reset_index(drop=True)

    md = load_md(a.window)
    keep = [c for c in (usable_descriptors(a.window) or []) if c in md.columns]
    if "n_heavy" in md.columns and "n_heavy" not in keep:
        keep = keep + ["n_heavy"]

    # 2D arms, built once over all dev ligands and indexed per subset
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    from rdkit.ML.Descriptors import MoleculeDescriptors
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    names = [d[0] for d in Descriptors._descList]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
    smi = dict(zip(lig.molecule_chembl_id, lig.smiles_std))
    mols = [Chem.MolFromSmiles(smi[c]) for c in all_ids]
    cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], float)
    rdk = np.nan_to_num(np.array([calc.CalcDescriptors(m) for m in mols], float),
                        nan=0.0, posinf=0.0, neginf=0.0)
    dev_index = {c: i for i, c in enumerate(all_ids)}
    X1s = [np.array([gen.GetFingerprintAsNumPy(m) for m in mols], float), np.hstack([cnt, rdk])]

    heavy = dict(zip(md.index, md["n_heavy"])) if "n_heavy" in md.columns else {}
    # replicate count and native trajectory length per ligand, recorded because they travel
    # with arrival order and would otherwise be invisible co-explanations
    nrep, nsa_acc = {}, {}
    for f in glob.glob(str(RES/"16_rung4_per_run"/"*.json")):
        import json as _j
        d_ = _j.loads(Path(f).read_text())
        l_ = re.sub(r"_r\d+$", "", d_["run"])
        nrep[l_] = nrep.get(l_, 0) + 1
        nsa_acc.setdefault(l_, []).append(float(d_.get("ns_analysed", np.nan)))
    nsmean = {k: float(np.nanmean(v)) for k, v in nsa_acc.items()}
    order = completion_order(set(all_ids))
    print(f"  completion order recovered for {len(order)} of {len(all_ids)} dev ligands")

    rng = np.random.default_rng(SEED)
    rows = []
    for tgt in TARGET_PAIRS:
        # smallest prefix of the real completion order reaching this many pairs
        n_lig = next((k for k in range(5, len(order) + 1)
                      if len(dpair[dpair.ligand_a.isin(order[:k])
                                   & dpair.ligand_b.isin(order[:k])]) >= tgt), len(order))
        ids = order[:n_lig]
        r4, r1, npair = score_subset(ids, dev, dpair, keep, md, X1s, dev_index)
        rep = [nrep.get(i, 1) for i in ids]
        nsa = [nsmean.get(i, np.nan) for i in ids]
        cov = [heavy[i] for i in ids if i in heavy]
        unc = [heavy[i] for i in all_ids if i not in set(ids) and i in heavy]
        p = stats.mannwhitneyu(cov, unc).pvalue if unc else np.nan
        rows.append(dict(kind="arrival order", target_pairs=tgt, n_ligands=n_lig,
                         n_pairs=npair, rung4=r4, rung1_matched=r1,
                         gain_naive=r4 - PUBLISHED_RUNG1_B, gain_matched=r4 - r1,
                         median_heavy_covered=float(np.median(cov)) if cov else np.nan,
                         median_heavy_rest=float(np.median(unc)) if unc else np.nan,
                         size_bias_p=p, mean_replicates=float(np.mean(rep)),
                         mean_ns_analysed=float(np.nanmean(nsa)), n_draws=1))
        print(f"    arrival  {npair:3d} pairs / {n_lig:3d} ligands: rung4 {r4:+.3f}, "
              f"matched rung1 {r1:+.3f}, naive gain {r4-PUBLISHED_RUNG1_B:+.3f}, "
              f"matched gain {r4-r1:+.3f}")

        if n_lig < len(all_ids):
            g4, g1, gn, gm, ps = [], [], [], [], []
            for _ in range(a.n_random):
                pick = list(rng.choice(all_ids, n_lig, replace=False))
                b4, b1, bn = score_subset(pick, dev, dpair, keep, md, X1s, dev_index)
                if not np.isfinite(b4) or not np.isfinite(b1):
                    continue
                g4.append(b4); g1.append(b1)
                gn.append(b4 - PUBLISHED_RUNG1_B); gm.append(b4 - b1)
                cv_ = [heavy[i] for i in pick if i in heavy]
                un_ = [heavy[i] for i in all_ids if i not in set(pick) and i in heavy]
                ps.append(stats.mannwhitneyu(cv_, un_).pvalue)
            rows.append(dict(kind="random", target_pairs=tgt, n_ligands=n_lig,
                             n_pairs=np.nan, rung4=np.mean(g4), rung1_matched=np.mean(g1),
                             gain_naive=np.mean(gn), gain_matched=np.mean(gm),
                             median_heavy_covered=np.nan, median_heavy_rest=np.nan,
                             size_bias_p=float(np.median(ps)), mean_replicates=np.nan,
                             mean_ns_analysed=np.nan, n_draws=len(g4)))
            print(f"    random   mean over {len(g4)} draws: rung4 {np.mean(g4):+.3f}, "
                  f"matched rung1 {np.mean(g1):+.3f}, naive gain {np.mean(gn):+.3f}, "
                  f"matched gain {np.mean(gm):+.3f}")

    out = pd.DataFrame(rows)
    out.to_csv(RES/"21_subset_inflation.csv", index=False)
    print(f"\n  -> results/21_subset_inflation.csv")

    tex = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
           r"Subset & Ligands & Pairs & Rung 4 & Rung 1 & Gain vs & Gain vs \\",
           r" & & & $R^2$ & matched $R^2$ & published & matched \\", r"\midrule"]
    for _, r in out.iterrows():
        pr = "--" if not np.isfinite(r.n_pairs) else f"{int(r.n_pairs)}"
        tex.append(f"{r['kind']} & {int(r.n_ligands)} & {pr} & {r.rung4:+.3f} & "
                   f"{r.rung1_matched:+.3f} & {r.gain_naive:+.3f} & {r.gain_matched:+.3f} \\\\")
    tex += [r"\bottomrule", r"\end{tabular}"]
    (RES/"21_subset_inflation.tex").write_text("\n".join(tex) + "\n")
    print(f"  -> results/21_subset_inflation.tex")


if __name__ == "__main__":
    main()
