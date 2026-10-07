"""Stage 19 — does rung 4 (MD dynamics) beat rung 1 (2D) on the pair endpoint?

Scored to be comparable by construction with scripts/rung2_ensemble_model.py: the same
frozen folds (seed 20260902), the same six algorithms, the same antisymmetric pair
construction for framing B, the same preprocessing pipeline.

Three things this must get right, and each has already bitten the project once:

**One window for every run.** Production was 50 ns until 2026-10-02 and 25 ns after, so
`ns_analysed` is 45 for some runs and 20 for others. Scoring the full trajectory of each
would mean comparing ligands on unequal evidence. Every JSON carries the ablation windows,
so the table is built on the common 20 ns window for all of them. Runs that have only
`windows["20"]` and runs that also have 30 and 45 contribute the same quantity.

**The replicate filter.** The pilot (results/18_replicate_pilot.md) found 15 of ~30
descriptors measure the random seed more than the ligand. Feeding all thirty to a model
trained on ~100 pairs is an invitation to fit noise. The filtered set is the primary
result; the unfiltered one is reported beside it so the filter's effect is visible rather
than assumed.

**A matched baseline.** Rung 1 scored +0.200 over 145 dev pairs. MD does not yet cover all
217 ligands, so rung 4 is scored on whatever subset is complete. Comparing that number to
0.200 would be comparing two different problems. Rung 1 is therefore recomputed here on
exactly the pairs rung 4 can see, and that is the number to compare against.

Usage:  python scripts/rung4_model.py [--window 20] [--min-pairs 30]
"""
from __future__ import annotations
import argparse, glob, json, os, re, sys, warnings

# The 300-replicate bootstrap on 2026-10-08 wrote 2.37 million lines, 355 MB, of which 2.34
# million were one repeated sklearn UserWarning about sklearn.utils.parallel.delayed, and 49
# were the output. GitHub refused the push.
#
# I first added a message-specific filter here and said in the commit that it fixed this. It
# does not. Running 3 replicates with and without it gives 57 lines either way, because the
# blanket filterwarnings("ignore") a few lines below was already in force during the run that
# produced the 355 MB -- so the filter was never what was missing, and the warnings escaped it
# by a route I have not identified. Most likely they come from joblib worker processes, which do
# not inherit the parent's filters, but PYTHONWARNINGS=ignore changed nothing either, so that
# is a guess and not a finding.
#
# The protection is therefore not a filter: results/19_rung4_w20_*.log is gitignored, and the
# reading itself is recorded in 19b_rung4_85pairs.md. Check a log's size before adding it.
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.cross_decomposition import PLSRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.ensemble import (RandomForestRegressor, HistGradientBoostingRegressor,
                              RandomForestClassifier, HistGradientBoostingClassifier)
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.metrics import (r2_score, mean_absolute_error, roc_auc_score,
                             average_precision_score, matthews_corrcoef,
                             balanced_accuracy_score)
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
TARGET, SEED, NF = "CHEMBL220", 20260902, 5
# Carried for the size adjustment, never as a dynamic descriptor in its own right.
SIZE = "n_heavy"
META = {"run", "n_frames", "ns_analysed", "dt_ps", "windows"}


def pipe(est, scale=True):
    steps = [("imp", SimpleImputer(strategy="median")), ("vt", VarianceThreshold())]
    if scale: steps.append(("sc", StandardScaler()))
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

# Copied from rung2_model.py rather than re-chosen, so framing C at rung 4 is comparable with
# rung 2 by construction and not only by intention. Five families, not six: there is no
# classification counterpart to PLS in the rung-2 table either.
CLF = {
 "Logistic regression": lambda: pipe(LogisticRegression(max_iter=3000)),
 "k-NN (k=5)":          lambda: pipe(KNeighborsClassifier(n_neighbors=5)),
 "SVM (RBF)":           lambda: pipe(SVC(C=10.0, probability=True, random_state=SEED)),
 "Random Forest":       lambda: pipe(RandomForestClassifier(n_estimators=500,
                            random_state=SEED, n_jobs=-1), False),
 "Gradient Boosting":   lambda: pipe(HistGradientBoostingClassifier(random_state=SEED), False),
}


def load_md(window: str) -> pd.DataFrame:
    """One row per ligand, descriptors taken from the given ablation window.

    Replicates are averaged. That is the right estimator for a ligand property and it is
    also what the replicate pilot measured the noise of: averaging n replicates divides the
    seed variance by n, so ligands with three runs carry less noise than ligands with one.
    """
    rows = []
    for f in sorted(glob.glob(str(RES/"16_rung4_per_run"/"*.json"))):
        d = json.loads(Path(f).read_text())
        w = (d.get("windows") or {}).get(window)
        if w is None:
            continue                       # run too short to reach this window
        rec = {k: v for k, v in w.items() if k not in META}
        rec["molecule_chembl_id"] = re.sub(r"_r\d+$", "", d["run"])
        rec["replicate"] = int(re.search(r"_r(\d+)$", d["run"]).group(1))
        rows.append(rec)
    if not rows:
        sys.exit(f"no runs reach the {window} ns window")
    df = pd.DataFrame(rows)
    n_rep = df.groupby("molecule_chembl_id").size()
    agg = df.drop(columns=["replicate"]).groupby("molecule_chembl_id").mean(numeric_only=True)
    agg["n_replicates"] = n_rep
    return agg


def usable_descriptors(window: str) -> list[str] | None:
    """Descriptors the replicate pilot cleared at this window, or None if it has not run."""
    p = RES/"18_replicate_pilot.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    d = d[(d.window.astype(float) == float(window)) & (d.verdict == "usable")]
    return sorted(set(d.descriptor)) or None


def pair_rows(X, y, _unused, ids, dpair):
    """Antisymmetric pair construction, identical to stage 15."""
    idx = {c: i for i, c in enumerate(ids)}
    xs, ys, fs, used = [], [], [], []
    for _, r in dpair.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx:
            continue
        a, b = idx[r.ligand_a], idx[r.ligand_b]
        xs += [X[a]-X[b], X[b]-X[a]]
        ys += [y[a]-y[b], y[b]-y[a]]
        fs += [r.cv_fold]*2
        used.append(r.pair_id if "pair_id" in r else f"{r.ligand_a}|{r.ligand_b}")
    return np.array(xs), np.array(ys), np.array(fs), used


def cv(X, y, f, make, clf=False):
    out = np.zeros(len(y))
    for k in range(NF):
        tr, te = f != k, f == k
        if te.sum() == 0 or tr.sum() == 0:
            continue
        if clf and len(np.unique(y[tr])) < 2:
            continue                      # a fold with one class cannot be fitted
        m = make(); m.fit(X[tr], y[tr])
        out[te] = (m.predict_proba(X[te])[:, 1] if clf
                   else np.asarray(m.predict(X[te])).ravel())
    return out


def sym_pairs(X, ids, dpair):
    """Framing C rows: order-invariant pair features and the cliff label.

    [da+db, |da-db|] is the rung-2 construction, kept verbatim. It has to be symmetric --
    "is this pair a cliff" does not depend on which ligand is written first, and the
    antisymmetric rows used for framing B would make the label ambiguous.
    """
    idx = {c: i for i, c in enumerate(ids)}
    xs, ys, fs = [], [], []
    for _, r in dpair.iterrows():
        if r.ligand_a not in idx or r.ligand_b not in idx:
            continue
        da, db = X[idx[r.ligand_a]], X[idx[r.ligand_b]]
        xs.append(np.concatenate([da + db, np.abs(da - db)]))
        ys.append(1 if r.group == "cliff" else 0)
        fs.append(r.cv_fold)
    return np.array(xs), np.array(ys), np.array(fs)


def score_clf(tag, X, y, f, rows_out, framing="C"):
    best = None
    print(f"\n  {tag}   ({X.shape[0]} baris, {X.shape[1]} fitur, "
          f"{int(y.sum())} cliff / {len(y)-int(y.sum())} kontrol)")
    print(f"  {'algoritma':22s}{'AUC':>9s}{'AP':>9s}{'MCC':>9s}{'bal.acc':>9s}")
    for a, mk in CLF.items():
        pr = cv(X, y, f, mk, clf=True)
        hard = (pr >= 0.5).astype(int)
        auc = roc_auc_score(y, pr); ap = average_precision_score(y, pr)
        mcc = matthews_corrcoef(y, hard); ba = balanced_accuracy_score(y, hard)
        print(f"  {a:22s}{auc:9.3f}{ap:9.3f}{mcc:+9.3f}{ba:9.3f}")
        rows_out.append(dict(framing=framing, descriptors=tag, algorithm=a, auc_roc=auc,
                             auc_pr=ap, mcc=mcc, bal_acc=ba,
                             n_rows=X.shape[0], n_features=X.shape[1]))
        if best is None or auc > best[1]:
            best = (a, auc)
    return best


def score_reg_A(tag, X, y, f, rows_out):
    best = None
    print(f"\n  {tag}   ({X.shape[0]} ligan, {X.shape[1]} fitur)")
    print(f"  {'algoritma':22s}{'R2':>9s}{'rho':>9s}{'MAE':>9s}")
    for a, mk in REG.items():
        p = cv(X, y, f, mk)
        r2 = r2_score(y, p); rho = spearmanr(y, p).statistic
        mae = mean_absolute_error(y, p)
        print(f"  {a:22s}{r2:9.3f}{rho:+9.3f}{mae:9.3f}")
        rows_out.append(dict(framing="A", descriptors=tag, algorithm=a, r2=r2, rho=rho,
                             mae=mae, n_rows=X.shape[0], n_features=X.shape[1]))
        if best is None or r2 > best[1]:
            best = (a, r2)
    return best


def boot_gap(unit_rows, build, metric, n_units, n_boot, seed, clf=False):
    """CI on a best-of-family gap, resampling whole units with replacement.

    `unit_rows` maps a unit index to the row indices it owns, so a pair contributes both of
    its antisymmetric rows together and a ligand contributes its single row. `build` returns
    (X, y, fold) for a list of row indices, for each of the two rungs in turn.

    One limitation, stated rather than hidden: resampling observations and then running the
    frozen 5-fold CV over the resample puts duplicated units on both sides of the train/test
    boundary, which inflates the absolute scores. It is applied identically to both rungs, so
    the gap -- the quantity actually read -- is far less affected than the levels. The levels
    reported from a bootstrap median should not be quoted as the study's performance; the
    point estimates from the unresampled fit are those.
    """
    rng = np.random.default_rng(seed)
    out = {"a": [], "b": [], "gap": []}
    for _ in range(n_boot):
        take = rng.integers(0, n_units, n_units)
        rows = np.concatenate([unit_rows[u] for u in take])
        best = []
        for which in (0, 1):
            X, y, f = build(which, rows)
            vals = []
            for mk in (CLF if clf else REG).values():
                try:
                    p = cv(X, y, f, mk, clf=clf)
                    vals.append(metric(y, p))
                except Exception:
                    vals.append(-np.inf)
            best.append(max(vals))
        out["a"].append(best[0]); out["b"].append(best[1])
        out["gap"].append(best[0] - best[1])
    q = lambda v: np.percentile(v, [2.5, 50, 97.5])
    return q(out["a"]), q(out["b"]), q(out["gap"])


def score(tag, X, y, f, rows_out):
    best = None
    print(f"\n  {tag}   ({X.shape[0]} baris, {X.shape[1]} fitur)")
    print(f"  {'algoritma':22s}{'R2':>9s}{'rho':>9s}{'MAE':>9s}")
    for a, mk in REG.items():
        p = cv(X, y, f, mk)
        r2 = r2_score(y, p); rho = spearmanr(y, p).statistic
        mae = mean_absolute_error(y, p)
        print(f"  {a:22s}{r2:9.3f}{rho:+9.3f}{mae:9.3f}")
        rows_out.append(dict(descriptors=tag, algorithm=a, r2_B=r2, rho_B=rho, mae_B=mae,
                             n_rows=X.shape[0], n_features=X.shape[1]))
        if best is None or r2 > best[1]:
            best = (a, r2)
    return best


# This lived below the __main__ block until 2026-10-08, which meant it was never defined
# while main ran and never called from anywhere -- the one function that answers the study's
# actual question, "is the gap's lower bound above zero", was dead code through the whole
# preliminary reading. A function placed after the block that uses it is not a late
# definition, it is an absent one.

def bootstrap_gap(md, keep, ids, y, dpair, X1, n_boot=300, seed=SEED):
    """CI on the rung4 - rung1 gap, resampling PAIRS with replacement.

    Two things this does deliberately. It resamples pairs, not rows, because the two
    antisymmetric rows of a pair are one observation. And it re-selects the best algorithm
    inside every replicate for both rungs, because picking the best of six is part of the
    procedure being measured -- fixing the winner chosen on the full sample would report a
    CI for a procedure nobody ran.
    """
    rng = np.random.default_rng(seed)
    Xm = md.loc[ids, keep].to_numpy(float)
    idx = {c: i for i, c in enumerate(ids)}
    have = set(ids)
    cov = dpair[dpair.ligand_a.isin(have) & dpair.ligand_b.isin(have)].reset_index(drop=True)
    gaps, r4s, r1s = [], [], []
    for _ in range(n_boot):
        take = cov.iloc[rng.integers(0, len(cov), len(cov))]
        best = []
        for X in (Xm, X1):
            xs, ys, fs = [], [], []
            for _, r in take.iterrows():
                a, b = idx[r.ligand_a], idx[r.ligand_b]
                xs += [X[a]-X[b], X[b]-X[a]]
                ys += [y[a]-y[b], y[b]-y[a]]
                fs += [r.cv_fold]*2
            Xp, yp, fp = np.array(xs), np.array(ys), np.array(fs)
            r2s = []
            for mk in REG.values():
                try:
                    r2s.append(r2_score(yp, cv(Xp, yp, fp, mk)))
                except Exception:
                    r2s.append(-np.inf)
            best.append(max(r2s))
        r4s.append(best[0]); r1s.append(best[1]); gaps.append(best[0]-best[1])
    q = lambda v: np.percentile(v, [2.5, 50, 97.5])
    return q(r4s), q(r1s), q(gaps)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", default="20", help="ablation window in ns (default 20)")
    ap.add_argument("--min-pairs", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=300,
                    help="bootstrap replicates for the gap CI; 0 skips")
    ap.add_argument("--framings", default="ABC",
                    help="which framings to score, e.g. B or AC (default ABC)")
    a = ap.parse_args()

    CI_ROWS = []
    sl = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    lig = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
    lig = lig.merge(sl[["molecule_chembl_id", "split", "cv_fold"]], on="molecule_chembl_id")
    dpair = pairs[pairs.split == "dev"].reset_index(drop=True)

    md = load_md(a.window)
    dev = lig[(lig.split == "dev") & (lig.molecule_chembl_id.isin(md.index))].reset_index(drop=True)
    ids = list(dev.molecule_chembl_id)
    y = dev.pchembl.to_numpy()

    have = set(ids)
    covered = dpair[dpair.ligand_a.isin(have) & dpair.ligand_b.isin(have)]
    print("=" * 74)
    print(f"RUNG 4 — jendela {a.window} ns")
    print("=" * 74)
    print(f"  ligan dev dengan MD   : {len(dev)} / {(lig.split=='dev').sum()}")
    print(f"  replika per ligan     : {md.loc[ids,'n_replicates'].value_counts().to_dict()}")
    print(f"  pasangan dev lengkap  : {len(covered)} / {len(dpair)}")
    if len(covered) < a.min_pairs:
        sys.exit(f"\n  BERHENTI: {len(covered)} pasangan terlalu sedikit untuk angka yang berarti "
                 f"(ambang {a.min_pairs}).\n  Jalankan lagi setelah lebih banyak run selesai.")

    keep = usable_descriptors(a.window)
    allcols = [c for c in md.columns if c not in ("n_replicates",) and md[c].dtype != object]
    if keep:
        keep = [c for c in keep if c in allcols]
        if SIZE in allcols and SIZE not in keep:
            keep = keep + [SIZE]           # size is a covariate, not a dynamic descriptor
    print(f"  deskriptor tersedia   : {len(allcols)}")
    print(f"  lolos saringan replika: {len(keep) if keep else 'uji replika belum ada'}")

    rows = []
    sets = {"rung4 (tersaring)": keep, "rung4 (semua deskriptor)": allcols}
    results = {}
    for tag, cols in sets.items():
        if not cols:
            continue
        X = md.loc[ids, cols].to_numpy(float)
        Xp, yp, fp, _ = pair_rows(X, y, None, ids, dpair)
        results[tag] = score(tag, Xp, yp, fp, rows)

    # ---- matched rung-1 baseline, on exactly these pairs ----
    X1 = None
    # Comparing rung 4 on a subset against rung 1's published 0.200 over 145 pairs would be
    # comparing two different problems. The baseline is recomputed on the same pairs.
    try:
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
        X1 = np.hstack([cnt, rdk])
        Xp1, yp1, fp1, _ = pair_rows(X1, y, None, ids, dpair)
        results["rung1 2D (pasangan yang sama)"] = score(
            "rung1 2D (pasangan yang sama)", Xp1, yp1, fp1, rows)
    except Exception as e:
        print(f"\n  baseline rung 1 dilewati: {type(e).__name__}: {e}")

    out = pd.DataFrame(rows)
    out.to_csv(RES/f"19_rung4_benchmark_w{a.window}.csv", index=False)
    print("\n" + "=" * 74)
    print(f"RINGKASAN — framing B, {len(covered)} pasangan, jendela {a.window} ns")
    print("=" * 74)
    for tag, (alg, r2) in results.items():
        print(f"  {tag:30s} R2 = {r2:+.3f}   ({alg})")
    if "rung4 (tersaring)" in results and "rung1 2D (pasangan yang sama)" in results:
        d = results["rung4 (tersaring)"][1] - results["rung1 2D (pasangan yang sama)"][1]
        print(f"\n  selisih rung4 - rung1 : {d:+.3f}")
        print(f"  (analisis daya stage 17: studi ini memisahkan R2 yang berjarak ~0,30;")
        print(f"   selisih di bawah itu tidak bisa dibedakan dari nol pada jumlah pasangan penuh,")
        print(f"   apalagi pada {len(covered)})")
        if a.n_boot and keep and X1 is not None:
            print(f"\n  bootstrap atas pasangan, {a.n_boot} replika, algoritma dipilih ulang")
            print("  di dalam setiap replika:")
            q4, q1, qg = bootstrap_gap(md, keep, ids, y, dpair, X1,
                                       n_boot=a.n_boot, seed=SEED)
            for nm, q in (("rung4", q4), ("rung1 2D", q1), ("selisih", qg)):
                print(f"    {nm:9s} median {q[1]:+.3f}   95% CI [{q[0]:+.3f}, {q[2]:+.3f}]")
            lo, hi = qg[0], qg[2]
            if lo > 0:
                verdict = "batas bawah DI ATAS nol -- rung 4 menambah sesuatu"
            elif hi < 0:
                verdict = "batas atas DI BAWAH nol -- rung 4 lebih buruk"
            else:
                verdict = ("selang masih memuat nol -- pada jumlah pasangan ini studi tidak "
                           "bisa memisahkan keduanya")
            print(f"\n  -> {verdict}")
            pd.DataFrame({"quantity": ["rung4", "rung1_2d", "gap"],
                          "lo2.5": [q4[0], q1[0], qg[0]],
                          "median": [q4[1], q1[1], qg[1]],
                          "hi97.5": [q4[2], q1[2], qg[2]],
                          "n_pairs": len(covered), "n_boot": a.n_boot,
                          "window_ns": a.window}).to_csv(
                RES/f"19_rung4_gap_ci_w{a.window}.csv", index=False)
    print(f"\n  -> {RES.name}/19_rung4_benchmark_w{a.window}.csv")
    # ================================================================================
    # Framings A and C. Until 2026-10-08 rung 4 had only ever been scored on framing B,
    # the signed per-pair difference -- which is both the primary endpoint and the hardest
    # of the three. A null on B is not a null on A or C, and reporting it as one would be
    # a straightforwardly false claim about what was tested.
    # ================================================================================
    fold = dev.cv_fold.to_numpy()
    md_keep = md.loc[ids, keep].to_numpy(float) if keep else None
    md_all = md.loc[ids, allcols].to_numpy(float)

    if "A" in a.framings.upper():
        print("\n" + "=" * 74)
        print(f"FRAMING A — regresi per ligan, {len(ids)} ligan, jendela {a.window} ns")
        print("=" * 74)
        rowsA, resA = [], {}
        if md_keep is not None:
            resA["rung4 (tersaring)"] = score_reg_A("rung4 (tersaring)", md_keep, y, fold, rowsA)
        resA["rung4 (semua deskriptor)"] = score_reg_A("rung4 (semua deskriptor)",
                                                      md_all, y, fold, rowsA)
        if X1 is not None:
            resA["rung1 2D (ligan yang sama)"] = score_reg_A("rung1 2D (ligan yang sama)",
                                                            X1, y, fold, rowsA)
        pd.DataFrame(rowsA).to_csv(RES/f"19_rung4_framingA_w{a.window}.csv", index=False)
        print("\n  ringkasan framing A")
        for tag, (alg, v) in resA.items():
            print(f"    {tag:30s} R2 = {v:+.3f}   ({alg})")
        if md_keep is not None and X1 is not None:
            gap = resA["rung4 (tersaring)"][1] - resA["rung1 2D (ligan yang sama)"][1]
            print(f"    selisih rung4 - rung1 : {gap:+.3f}")
            if a.n_boot:
                units = [np.array([i]) for i in range(len(ids))]
                mats = (md_keep, X1)
                qa, qb, qg = boot_gap(units,
                                      lambda w, r: (mats[w][r], y[r], fold[r]),
                                      r2_score, len(units), a.n_boot, SEED)
                for nm, q in (("rung4", qa), ("rung1 2D", qb), ("selisih", qg)):
                    print(f"    {nm:9s} median {q[1]:+.3f}   95% CI [{q[0]:+.3f}, {q[2]:+.3f}]")
                CI_ROWS.append(dict(framing="A", metric="r2", n_units=len(units),
                                    gap_lo=qg[0], gap_med=qg[1], gap_hi=qg[2],
                                    rung4_med=qa[1], rung1_med=qb[1]))

    if "C" in a.framings.upper() and X1 is not None:
        print("\n" + "=" * 74)
        print(f"FRAMING C — klasifikasi cliff per pasangan, jendela {a.window} ns")
        print("=" * 74)
        rowsC, resC = [], {}
        if md_keep is not None:
            Xc, yc, fc = sym_pairs(md_keep, ids, dpair)
            resC["rung4 (tersaring)"] = score_clf("rung4 (tersaring)", Xc, yc, fc, rowsC)
        Xca, yca, fca = sym_pairs(md_all, ids, dpair)
        resC["rung4 (semua deskriptor)"] = score_clf("rung4 (semua deskriptor)",
                                                    Xca, yca, fca, rowsC)
        Xc1, yc1, fc1 = sym_pairs(X1, ids, dpair)
        resC["rung1 2D (pasangan yang sama)"] = score_clf("rung1 2D (pasangan yang sama)",
                                                         Xc1, yc1, fc1, rowsC)
        pd.DataFrame(rowsC).to_csv(RES/f"19_rung4_framingC_w{a.window}.csv", index=False)
        print("\n  ringkasan framing C")
        for tag, (alg, v) in resC.items():
            print(f"    {tag:30s} AUC = {v:.3f}   ({alg})")
        if md_keep is not None:
            gap = resC["rung4 (tersaring)"][1] - resC["rung1 2D (pasangan yang sama)"][1]
            print(f"    selisih AUC rung4 - rung1 : {gap:+.3f}")
            if a.n_boot:
                units = [np.array([i]) for i in range(len(yc))]
                mats = (Xc, Xc1)
                qa, qb, qg = boot_gap(units,
                                      lambda w, r: (mats[w][r], yc[r], fc[r]),
                                      roc_auc_score, len(units), a.n_boot, SEED, clf=True)
                for nm, q in (("rung4", qa), ("rung1 2D", qb), ("selisih", qg)):
                    print(f"    {nm:9s} median {q[1]:+.3f}   95% CI [{q[0]:+.3f}, {q[2]:+.3f}]")
                CI_ROWS.append(dict(framing="C", metric="auc_roc", n_units=len(units),
                                    gap_lo=qg[0], gap_med=qg[1], gap_hi=qg[2],
                                    rung4_med=qa[1], rung1_med=qb[1]))

    if CI_ROWS:
        pd.DataFrame(CI_ROWS).to_csv(RES/f"19_rung4_gap_ci_all_w{a.window}.csv", index=False)
        print(f"\n  -> {RES.name}/19_rung4_gap_ci_all_w{a.window}.csv")
