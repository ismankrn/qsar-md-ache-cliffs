"""Stage 28 — the sealed test split, opened once.

The development split chose everything: which descriptor arm each rung uses, which algorithm
family, and the replicate filter that decides which dynamics descriptors are admissible. The
test split has been sealed since the beginning and has never been scored. This opens it.

Until 2026-10-25 it could not be opened honestly. Four of the five rungs covered all 271
ligands, but rung 4 had descriptors for **0 of the 54 test ligands** — the molecular dynamics
had only ever been run on the development set. Scoring the ladder on the test split would
have meant scoring four rungs and omitting the one the study is about. The 54 test ligands
were simulated for 25 ns each to close exactly that gap.

WHAT THIS IS NOT ALLOWED TO DO
------------------------------
Two rules, both fixed before the test ligands were simulated, and both of which make this a
weaker-looking but honest number rather than a stronger-looking dishonest one:

  * **The replicate filter is not recomputed.** It was fitted on the development replicates
    before any of this existed. Refitting it now on data that includes the test split would
    let the held-out ligands help choose the descriptors they are then used to judge. The
    descriptors that passed at the 20 ns window stay exactly as they were.

  * **No model is retrained or reselected.** Each rung keeps the descriptor arm and the
    algorithm family that won on development, as recorded in `results/20_ladder.csv`, refits
    that one combination on all 217 development ligands, and is scored once. There is no
    search here, no best-of, and no second look.

Both rules are enforced by construction: the arms and families are read from the ladder, and
nothing in this file chooses between them.

WHAT TO EXPECT
--------------
Intervals about twice as wide as the development ones. Framing B rests on 38 test pairs
against 145 development pairs, and an interval widens roughly as the square root of that
ratio. A test split this size cannot resolve a small difference; it can only confirm or
contradict a large one, which is the same limitation the development result already has and
the reason the detection limit was measured in advance.

Usage:  python scripts/test_split_eval.py [--window 20] [--n-boot 300]
"""
from __future__ import annotations
import argparse, glob, json, re, sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.cross_decomposition import PLSRegression
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.svm import SVR, SVC
from sklearn.ensemble import (RandomForestRegressor, RandomForestClassifier,
                              HistGradientBoostingRegressor, HistGradientBoostingClassifier)
from sklearn.metrics import r2_score, roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET = "CHEMBL220"
SEED = 20260902            # the project's frozen-split seed, reused for the bootstrap

# --------------------------------------------------------------------------- estimators
# Copied verbatim from rung4_model.py, which itself copied them from rung2_model.py. The
# point of the ladder is that the rungs differ in their descriptors and in nothing else, so
# re-choosing hyperparameters here would quietly break the comparison this script exists to
# make.
def pipe(est, scale=True):
    steps = [("imp", SimpleImputer(strategy="median")), ("vt", VarianceThreshold())]
    if scale:
        steps.append(("sc", StandardScaler()))
    steps.append(("est", est))
    return Pipeline(steps)

REG = {
 "Ridge (linear)":    lambda: pipe(Ridge(alpha=10.0)),
 "PLS":               lambda: pipe(PLSRegression(n_components=5)),
 "k-NN (k=5)":        lambda: pipe(KNeighborsRegressor(5)),
 "SVM (RBF)":         lambda: pipe(SVR(C=10.0, gamma="scale")),
 "Random Forest":     lambda: pipe(RandomForestRegressor(n_estimators=500,
                                    random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting": lambda: pipe(HistGradientBoostingRegressor(random_state=SEED), False),
}
CLF = {
 "Logistic regression": lambda: pipe(LogisticRegression(max_iter=3000)),
 "k-NN (k=5)":          lambda: pipe(KNeighborsClassifier(n_neighbors=5)),
 "SVM (RBF)":           lambda: pipe(SVC(C=10.0, probability=True, random_state=SEED)),
 "Random Forest":       lambda: pipe(RandomForestClassifier(n_estimators=500,
                            random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting":   lambda: pipe(HistGradientBoostingClassifier(random_state=SEED), False),
}


# --------------------------------------------------------------------------- descriptor arms
def arm_rung1(which: str, ids: list[str]) -> pd.DataFrame:
    """2D fingerprints, recomputed from SMILES exactly as stage 05 does."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdFingerprintGenerator
    from rdkit.ML.Descriptors import MoleculeDescriptors
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv").set_index("molecule_chembl_id")
    smiles = lig.loc[ids, "smiles_std"].tolist()
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    if which == "binary ECFP4":
        X = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], float)
    else:
        names = [d[0] for d in Descriptors._descList]
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(names)
        cnt = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], float)
        rdk = np.nan_to_num(np.array([calc.CalcDescriptors(m) for m in mols], float),
                            nan=0.0, posinf=0.0, neginf=0.0)
        X = np.hstack([cnt, rdk])
    return pd.DataFrame(X, index=ids)


def arm_rung2() -> pd.DataFrame:
    d2 = pd.read_csv(RES/"08_rung2_descriptors.csv")
    br = pd.read_csv(RES/"08_rung2_bridging.csv")
    d = d2.merge(br, on="molecule_chembl_id", how="left").set_index("molecule_chembl_id")
    return d.drop(columns=[c for c in d.columns
                           if c.startswith("spread_") or c == "n_seeds_used"], errors="ignore")


def arm_rung2e() -> pd.DataFrame:
    d = pd.read_csv(RES/"15_rung2_ensemble_median.csv").set_index("molecule_chembl_id")
    return d.drop(columns=[c for c in d.columns if c.startswith("spread_")], errors="ignore")


def arm_rung3() -> pd.DataFrame:
    d = (pd.read_csv(RES/"09_rung3_mmgbsa.csv")
           .drop_duplicates("molecule_chembl_id").set_index("molecule_chembl_id"))
    feat = [c for c in d.columns if c.startswith(("dE_", "dG_", "res_", "dSASA"))]
    return d[feat]


def _rung4_module():
    """rung4_model.py itself, imported rather than reimplemented.

    An earlier version of this file rebuilt the rung-4 table from the JSONs with its own
    list of metadata keys. The list was a guess, it missed `dt_ps`, and that key was then
    carried into the model as if it were a dynamic descriptor — 34 features where the
    development result had 33. The n_features check against the ladder caught it, but only
    because the ladder happened to record the number.

    Importing the real loader removes the whole class of error: the development result and
    the test result are now produced by the same code, so they cannot drift apart silently.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("rung4_model", ROOT/"scripts"/"rung4_model.py")
    mod = importlib.util.module_from_spec(spec)
    argv, sys.argv = sys.argv, ["rung4_model.py"]
    try:
        spec.loader.exec_module(mod)
    except SystemExit:
        pass                       # its __main__ guard, not an error
    finally:
        sys.argv = argv
    return mod


def arm_rung4(window: str, filtered: bool) -> pd.DataFrame:
    """One row per ligand, replicates averaged, built by rung4_model's own loader.

    `filtered` reproduces the development selection exactly: the descriptors the replicate
    pilot cleared at this window, plus the heavy-atom count as a size covariate. That
    covariate is the difference between "11 descriptors" and the 12 features the ladder
    records, and leaving it out would quietly score a different model than the one selected.
    """
    m = _rung4_module()
    md = m.load_md(window)
    allcols = [c for c in md.columns if c not in ("n_replicates",) and md[c].dtype != object]
    if not filtered:
        return md[allcols]
    keep = m.usable_descriptors(window)
    if not keep:
        sys.exit("uji replika tidak tersedia: filter replikat tidak boleh dihitung ulang di sini")
    keep = [c for c in keep if c in allcols]
    if m.SIZE in allcols and m.SIZE not in keep:
        keep = keep + [m.SIZE]
    return md[keep]


def build_arm(arm: str, ids: list[str], window: str) -> pd.DataFrame:
    if arm in ("binary ECFP4", "count ECFP4 + RDKit2D"):
        return arm_rung1(arm, ids)
    if arm == "4EY7 + bridging":
        return arm_rung2()
    if arm == "ensemble median":
        return arm_rung2e()
    if arm == "MM/GBSA":
        return arm_rung3()
    if arm.startswith("all descriptors"):
        return arm_rung4(window, filtered=False)
    if arm.startswith("replicate-filtered"):
        return arm_rung4(window, filtered=True)
    sys.exit(f"arm tidak dikenal: {arm}")


# --------------------------------------------------------------------------- pair matrices
def anti_pairs(X: np.ndarray, ids: list[str], dpair: pd.DataFrame):
    """Framing B: signed difference, written both ways round.

    Both orderings are kept so the model cannot learn a sign convention from the order the
    pairs happen to be written in. Identical to stage 15 and stage 19.
    """
    idx = {c: i for i, c in enumerate(ids)}
    xs, ys, pid = [], [], []
    for _, r in dpair.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx:
            continue
        a, b = idx[r.ligand_a], idx[r.ligand_b]
        xs += [X[a]-X[b], X[b]-X[a]]
        ys += [r.delta, -r.delta]
        pid += [r.pair_id, r.pair_id]
    return np.array(xs), np.array(ys), np.array(pid)


def sym_pairs(X: np.ndarray, ids: list[str], dpair: pd.DataFrame):
    """Framing C: order-invariant features and the cliff label."""
    idx = {c: i for i, c in enumerate(ids)}
    xs, ys, pid = [], [], []
    for _, r in dpair.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx:
            continue
        da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        xs.append(np.concatenate([da + db, np.abs(da - db)]))
        ys.append(1 if r.group == "cliff" else 0)
        pid.append(r.pair_id)
    return np.array(xs), np.array(ys), np.array(pid)


# --------------------------------------------------------------------------- bootstrap
def boot_ci(units, per_unit, metric, n_boot, rng):
    """Percentile interval over units — ligands for framing A, pairs for B and C.

    Resampling units rather than rows is what keeps the two antisymmetric rows of a pair
    together. Splitting them would treat one pair as two independent observations and
    narrow every interval by about sqrt(2) for nothing.
    """
    u = np.array(sorted(set(units)))
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(u, size=len(u), replace=True)
        rows = np.concatenate([per_unit[p] for p in pick])
        v = metric(rows)
        if v is not None and np.isfinite(v):
            vals.append(v)
    if len(vals) < 20:
        return np.nan, np.nan
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", default="20", help="jendela ablasi rung 4, ns")
    ap.add_argument("--n-boot", type=int, default=300)
    a = ap.parse_args()

    ladder = pd.read_csv(RES/"20_ladder.csv")
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id", "pchembl"]]
    lig = lig.merge(sl[["molecule_chembl_id", "split"]], on="molecule_chembl_id")

    dev_ids = lig[lig.split == "dev"].molecule_chembl_id.tolist()
    test_ids = lig[lig.split == "test"].molecule_chembl_id.tolist()
    y_all = dict(zip(lig.molecule_chembl_id, lig.pchembl))
    dpair = pairs[pairs.split == "dev"].reset_index(drop=True)
    tpair = pairs[pairs.split == "test"].reset_index(drop=True)

    print(f"pengembangan : {len(dev_ids)} ligan, {len(dpair)} pasangan")
    print(f"test         : {len(test_ids)} ligan, {len(tpair)} pasangan "
          f"({(tpair.group=='cliff').sum()} cliff / {(tpair.group=='control').sum()} kontrol)")
    print(f"bootstrap    : {a.n_boot} ulangan atas unit test\n")

    rng = np.random.default_rng(SEED)
    out = []
    preds = {}          # (framing, key) -> predictions on the test units, for the gaps below

    for framing in ("A", "B", "C"):
        print("=" * 78)
        print({"A": "FRAMING A — regresi per ligan, y = pChEMBL",
               "B": "FRAMING B — regresi per pasangan, y = ΔpChEMBL bertanda",
               "C": "FRAMING C — klasifikasi cliff / kontrol, AUC"}[framing])
        sub = ladder[ladder.framing == framing]
        for _, row in sub.iterrows():
            arm, algo = row.arm, row.algorithm
            D = build_arm(arm, dev_ids + test_ids, a.window)

            ids = [c for c in dev_ids + test_ids if c in D.index]
            dev_here = [c for c in dev_ids if c in D.index]
            test_here = [c for c in test_ids if c in D.index]
            X = D.loc[ids].to_numpy(float)
            yv = np.array([y_all[c] for c in ids])

            if framing == "A":
                tr = np.array([c in set(dev_here) for c in ids])
                Xtr, ytr = X[tr], yv[tr]
                Xte, yte = X[~tr], yv[~tr]
                units = np.array(ids)[~tr]
                make = REG[algo]
            else:
                Xtr, ytr, _ = (anti_pairs if framing == "B" else sym_pairs)(X, ids, dpair)
                Xte, yte, units = (anti_pairs if framing == "B" else sym_pairs)(X, ids, tpair)
                make = (REG if framing == "B" else CLF)[algo]

            n_feat_expected = int(row.n_features)
            if Xtr.shape[1] != n_feat_expected:
                print(f"  ! {row.rung}: {Xtr.shape[1]} fitur, ladder mencatat "
                      f"{n_feat_expected} -- konstruksi arm tidak cocok")

            m = make()
            m.fit(Xtr, ytr)
            pred = (m.predict_proba(Xte)[:, 1] if framing == "C"
                    else np.asarray(m.predict(Xte)).ravel())

            if framing == "C":
                value = roc_auc_score(yte, pred)
                metric = lambda rows: (roc_auc_score(yte[rows], pred[rows])
                                       if len(set(yte[rows])) > 1 else None)
            else:
                value = r2_score(yte, pred)
                metric = lambda rows: r2_score(yte[rows], pred[rows])

            per_unit = {u: np.where(units == u)[0] for u in set(units)}
            lo, hi = boot_ci(units, per_unit, metric, a.n_boot, rng)

            label = "AUC" if framing == "C" else "R²"
            print(f"  {row.rung:<17} {arm:<24} {algo:<19} "
                  f"{label} {value:6.3f}  [{lo:6.3f}, {hi:6.3f}]  "
                  f"n={len(set(units))}")
            preds[(framing, row.key)] = {"pred": pred, "y": yte, "units": units}
            out.append({"framing": framing, "rung": row.rung, "key": row.key,
                        "arm": arm, "algorithm": algo,
                        "n_features": Xtr.shape[1],
                        "dev_value": row.value, "test_value": value,
                        "test_ci_lo": lo, "test_ci_hi": hi,
                        "n_test_units": len(set(units))})
        print()

    df = pd.DataFrame(out)
    df.to_csv(RES/"28_test_split.csv", index=False)

    # ---- the gap against rung 1, which is the only quantity comparable across splits ----
    # Levels are not. Every rung scores higher on the test pairs than on the development
    # folds, rung 3 moving from -0.06 to +0.51, which says the test sample is easier rather
    # than that the ladder suddenly works. The difference between two rungs measured on the
    # SAME 38 pairs is insensitive to that shift, and it is what the development result
    # reports too, so it is what can be set beside it.
    #
    # Paired by unit: each bootstrap replicate resamples pairs once and scores both rungs on
    # that same resample. Resampling them independently would add variance that is not there
    # and widen every interval for nothing.
    print("=" * 78)
    print("SELISIH TERHADAP RUNG 1, pada unit test yang sama")
    gaps = []
    for framing in ("A", "B", "C"):
        base = preds[(framing, "r1")]
        for key in ("r2", "r2e", "r3", "r4"):
            if (framing, key) not in preds:
                continue
            cur = preds[(framing, key)]
            yte, units = cur["y"], cur["units"]
            if framing == "C":
                mt = lambda rows, p: (roc_auc_score(yte[rows], p[rows])
                                      if len(set(yte[rows])) > 1 else None)
            else:
                mt = lambda rows, p: r2_score(yte[rows], p[rows])
            per_unit = {u: np.where(units == u)[0] for u in set(units)}
            obs = mt(np.arange(len(yte)), cur["pred"]) - mt(np.arange(len(yte)), base["pred"])
            u = np.array(sorted(per_unit))
            vals = []
            for _ in range(a.n_boot):
                pick = rng.choice(u, size=len(u), replace=True)
                rows = np.concatenate([per_unit[p] for p in pick])
                va, vb = mt(rows, cur["pred"]), mt(rows, base["pred"])
                if va is not None and vb is not None and np.isfinite(va) and np.isfinite(vb):
                    vals.append(va - vb)
            lo, hi = ((np.percentile(vals, 2.5), np.percentile(vals, 97.5))
                      if len(vals) >= 20 else (np.nan, np.nan))
            dev_gap = (ladder[(ladder.framing == framing) & (ladder.key == key)].value.iloc[0]
                       - ladder[(ladder.framing == framing) & (ladder.key == "r1")].value.iloc[0])
            # Two-sided bootstrap p, floored at 1/B: a p of 0 would only mean the bootstrap
            # was not run long enough to see the other tail, and reporting it as 0 would turn
            # a resolution limit into a claim.
            v = np.asarray(vals)
            p = 2.0 * min((v <= 0).mean(), (v >= 0).mean())
            p = float(min(1.0, max(p, 1.0/len(v))))
            mark = "  interval tidak memuat nol" if (lo > 0 or hi < 0) else ""
            print(f"  {framing}  {key:<4} selisih test {obs:+6.3f}  [{lo:+6.3f}, {hi:+6.3f}]"
                  f"   selisih dev {dev_gap:+6.3f}   p={p:.4f}{mark}")
            gaps.append({"framing": framing, "key": key, "test_gap": obs,
                         "ci_lo": lo, "ci_hi": hi, "dev_gap": dev_gap, "p_boot": p})

    # ---- multiplicity ---------------------------------------------------------------
    # Twelve gaps are tested, so about 0.6 of them would clear a 5% threshold by chance
    # alone. Two corrections rather than one, because they answer different questions and
    # disagreeing is informative: Holm controls the chance of ANY false positive among the
    # twelve, BH controls the expected share of false positives among those declared.
    #
    # Both are applied to tests that are NOT independent -- all twelve share the same 38
    # pairs, and the rungs share pose information with each other. Holm is valid under any
    # dependence and is therefore conservative here; BH is valid under the positive
    # dependence these tests plausibly have. Neither is exact, and that is stated rather
    # than hidden by picking whichever one gives the nicer answer.
    g = pd.DataFrame(gaps).sort_values("p_boot").reset_index(drop=True)
    m = len(g)
    holm, run_max = [], 0.0
    for i, p in enumerate(g.p_boot):
        run_max = max(run_max, min(1.0, (m - i) * p))
        holm.append(run_max)
    g["p_holm"] = holm
    bh, run_min = [0.0]*m, 1.0
    for i in range(m - 1, -1, -1):
        run_min = min(run_min, min(1.0, m * g.p_boot[i] / (i + 1)))
        bh[i] = run_min
    g["p_bh"] = bh
    g["lolos_holm"] = g.p_holm < 0.05
    g["lolos_bh"] = g.p_bh < 0.05

    print()
    print("=" * 78)
    print(f"KOREKSI MULTIPLISITAS atas {m} selisih ({a.n_boot} ulangan bootstrap)")
    print(f"  {'':3}{'':5}{'selisih':>9}{'p':>9}{'p Holm':>9}{'p BH':>8}  keputusan")
    for _, r in g.iterrows():
        verdict = ("bertahan (Holm & BH)" if r.lolos_holm else
                   "bertahan hanya BH" if r.lolos_bh else "tidak bertahan")
        print(f"  {r.framing:<3}{r.key:<5}{r.test_gap:>+9.3f}{r.p_boot:>9.4f}"
              f"{r.p_holm:>9.4f}{r.p_bh:>8.4f}  {verdict}")
    print(f"\n  lolos Holm: {int(g.lolos_holm.sum())} dari {m}"
          f"   |   lolos BH: {int(g.lolos_bh.sum())} dari {m}")
    g.to_csv(RES/"28_test_split_gaps.csv", index=False)

    # ---- jackknife: is a surviving gap the work of two or three pairs? ----------------
    # Surviving a multiplicity correction says the gap is unlikely to be noise. It says
    # nothing about whether it rests on the whole test set or on a handful of pairs, and
    # with 38 pairs a couple of extreme ones can carry a difference this size on their own.
    # Each pair is dropped in turn and the gap recomputed. The models are not refitted: the
    # question here is how much of the measured difference depends on single observations,
    # not how the fit would have changed, and refitting would confound the two.
    jk = []
    print()
    print("=" * 78)
    print("JACKKNIFE, satu pasangan dibuang bergiliran")
    for _, r in g.iterrows():
        framing, key = r.framing, r.key
        cur, base = preds[(framing, key)], preds[(framing, "r1")]
        yte, units = cur["y"], cur["units"]
        if framing == "C":
            mt = lambda rows, p: (roc_auc_score(yte[rows], p[rows])
                                  if len(set(yte[rows])) > 1 else None)
        else:
            mt = lambda rows, p: r2_score(yte[rows], p[rows])
        per_unit = {u: np.where(units == u)[0] for u in set(units)}
        vals = {}
        for u in per_unit:
            rows = np.concatenate([ix for k2, ix in per_unit.items() if k2 != u])
            va, vb = mt(rows, cur["pred"]), mt(rows, base["pred"])
            if va is not None and vb is not None:
                vals[u] = va - vb
        if not vals:
            continue
        arr = np.array(list(vals.values()))
        worst = min(vals, key=vals.get)
        rec = {"framing": framing, "key": key, "gap": r.test_gap,
               "jk_min": arr.min(), "jk_max": arr.max(), "jk_median": float(np.median(arr)),
               "n_below_zero": int((arr <= 0).sum()), "n_pairs": len(arr),
               "most_influential": worst, "gap_without_it": vals[worst]}
        jk.append(rec)
        if r.lolos_holm or r.lolos_bh:
            print(f"  {framing} {key}  selisih {r.test_gap:+.3f}")
            print(f"     rentang jackknife : [{arr.min():+.3f}, {arr.max():+.3f}]  "
                  f"median {np.median(arr):+.3f}")
            print(f"     jatuh ke <= 0     : {int((arr<=0).sum())} dari {len(arr)} pasangan")
            print(f"     paling berpengaruh: {worst} -> selisih jadi {vals[worst]:+.3f}")
    jkd = pd.DataFrame(jk)
    jkd.to_csv(RES/"28_test_split_jackknife.csv", index=False)

    # ---- the manuscript table -------------------------------------------------------
    # Written from the same objects that were just printed, so the table cannot drift from
    # the numbers in the text the way a hand-copied one does.
    jkmap = {(r.framing, r.key): r for r in jkd.itertuples()} if len(jkd) else {}
    RUNG = {"r2": "Rung 2 pose", "r2e": "Rung 2 ensemble",
            "r3": "Rung 3 MM/GBSA", "r4": "Rung 4 MD"}
    order = [(f, k) for f in ("A", "B", "C") for k in ("r2", "r2e", "r3", "r4")]
    gi = g.set_index(["framing", "key"])
    lines = [r"\begin{tabular}{llr@{\,}lrrr}", r"\hline",
             r"Framing & Rung & \multicolumn{2}{c}{Test gap (95\% CI)} "
             r"& $p_{\text{Holm}}$ & $p_{\text{BH}}$ & Dev.\ gap \\", r"\hline"]
    for f, k in order:
        if (f, k) not in gi.index:
            continue
        r = gi.loc[(f, k)]
        star = r"^{*}" if r.p_holm < 0.05 else ""
        lines.append(
            f"{f} & {RUNG[k]} & ${r.test_gap:+.3f}{star}$ & "
            f"$[{r.ci_lo:+.3f}, {r.ci_hi:+.3f}]$ & "
            f"{r.p_holm:.3f} & {r.p_bh:.3f} & ${r.dev_gap:+.3f}$ \\\\")
    lines += [r"\hline", r"\end{tabular}"]
    (RES/"28_test_split_gaps.tex").write_text("\n".join(lines) + "\n")

    print(f"\n-> {RES.name}/28_test_split_jackknife.csv")
    print(f"-> {RES.name}/28_test_split_gaps.tex")
    print(f"\n-> {RES.name}/28_test_split.csv")
    print(f"-> {RES.name}/28_test_split_gaps.csv")
    return df


if __name__ == "__main__":
    main()
