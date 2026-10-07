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

# sklearn raises this from inside its own parallel helper, once per joblib task. A 300-replicate
# bootstrap over six algorithms emitted it 100,000+ times: 2.37 million lines, a 355 MB log that
# GitHub refused, and 49 lines of actual output buried in it. Silencing it loses nothing -- it
# reports an internal call convention, not anything about the fit.
warnings.filterwarnings("ignore", message=r".*sklearn\.utils\.parallel\.delayed.*")
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
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.metrics import r2_score, mean_absolute_error
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


def cv(X, y, f, make):
    out = np.zeros(len(y))
    for k in range(NF):
        tr, te = f != k, f == k
        if te.sum() == 0 or tr.sum() == 0:
            continue
        m = make(); m.fit(X[tr], y[tr])
        out[te] = np.asarray(m.predict(X[te])).ravel()
    return out


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
    a = ap.parse_args()

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