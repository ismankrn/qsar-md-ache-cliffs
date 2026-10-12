"""Stage 17b — power analysis through the ACTUAL pipeline, not around it.

The first attempt generated predictions directly and compared them to a fixed 0.200. Both
shortcuts flatter the result: real out-of-fold predictions carry the variance of fitting a
model on 4/5 of the data, and rung 1's R² is itself an estimate, so comparing against it as
though it were a constant overstates how sharply the two can be separated.

This version generates synthetic DESCRIPTOR MATRICES, pushes them through the same SVM and
the same frozen 5-fold CV that every rung was scored with, and compares each against a
synthetic baseline tuned to rung 1's 0.200 using a paired bootstrap over the 145 pairs.
Whatever power it reports is what this study actually has.
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
SEED, NF, NSIM, NBOOT, NFEAT = 20260902, 5, 150, 800, 4
rng = np.random.default_rng(SEED)

sl = pd.read_csv(RES/"03_splits_ligands_CHEMBL220.csv")
pairs = pd.read_csv(RES/"03_splits_pairs_CHEMBL220.csv")
dev_pairs = pairs[pairs.split == "dev"].reset_index(drop=True)
n_pairs = len(dev_pairs)

y1 = dev_pairs.delta.to_numpy()
y = np.concatenate([y1, -y1])                       # framing B, antisymmetric
fold = np.concatenate([dev_pairs.cv_fold.to_numpy()]*2)
pair_id = np.concatenate([np.arange(n_pairs)]*2)
rows_of = {p: np.where(pair_id == p)[0] for p in range(n_pairs)}

def mk(): return Pipeline([("sc", StandardScaler()), ("est", SVR(C=10.0, gamma="scale"))])

def cv_predict(X):
    out = np.zeros(len(y))
    for k in range(NF):
        tr, te = fold != k, fold == k
        if te.sum() == 0: continue
        m = mk(); m.fit(X[tr], y[tr]); out[te] = np.asarray(m.predict(X[te])).ravel()
    return out

def r2(yt, yp):
    return 1 - ((yt-yp)**2).sum()/((yt-yt.mean())**2).sum()

def synth_X(rho):
    """Antisymmetric descriptor block whose signal strength is set by rho."""
    z = (y - y.mean())/y.std()
    sig = rho*z[:, None] + np.sqrt(max(1e-9, 1-rho**2))*rng.standard_normal((len(y), NFEAT))
    return sig

# calibrate rho -> achieved out-of-fold R2, so results are reported against what the
# pipeline actually scores rather than a nominal target
print("kalibrasi rho -> R2 out-of-fold (rata-rata 25 ulangan)")
grid = {}
# Fewer features and a lower grid: with 12 correlated features the SVM averages them
# and the weakest setting already scored 0.309, overshooting the 0.200 baseline and
# leaving the decision-relevant 0.25-0.40 band untested.
for rho in (0.20, 0.26, 0.32, 0.38, 0.44, 0.52):
    vals = [r2(y, cv_predict(synth_X(rho))) for _ in range(25)]
    grid[rho] = float(np.mean(vals))
    print(f"  rho {rho:.2f} -> R2 {grid[rho]:+.3f}")

BASE_RHO = min(grid, key=lambda r: abs(grid[r]-0.200))
print(f"\nbaseline rung 1 disimulasikan pada rho {BASE_RHO:.2f} (R2 ~ {grid[BASE_RHO]:+.3f})")
print()
rows = []   # collected so fig_power.py reads the same numbers this prints
print(f"{'R2 kandidat':>13s}{'selisih rata2':>15s}{'daya (CI>0)':>14s}{'CI lebar':>11s}")
for rho, target in sorted(grid.items()):
    if rho == BASE_RHO: continue
    hits, diffs, widths = 0, [], []
    for _ in range(NSIM):
        pa = cv_predict(synth_X(rho))          # kandidat, mis. rung 4
        pb = cv_predict(synth_X(BASE_RHO))     # baseline, rung 1
        d = r2(y, pa) - r2(y, pb); diffs.append(d)
        b = np.empty(NBOOT)
        for k in range(NBOOT):
            samp = rng.integers(0, n_pairs, n_pairs)
            ii = np.concatenate([rows_of[p] for p in samp])
            b[k] = r2(y[ii], pa[ii]) - r2(y[ii], pb[ii])
        lo, hi = np.percentile(b, [2.5, 97.5])
        widths.append(hi-lo)
        if lo > 0: hits += 1
    rows.append(dict(rho=rho, candidate_r2=target, baseline_r2=grid[BASE_RHO],
                     mean_gain=float(np.mean(diffs)), power=hits / NSIM,
                     ci_width=float(np.mean(widths)), n_pairs=n_pairs,
                     n_sim=NSIM, n_boot=NBOOT))
    print(f"{target:>13.3f}{np.mean(diffs):>15.3f}{100*hits/NSIM:>13.0f}%{np.mean(widths):>11.3f}")
print(f"\n  NSIM={NSIM}, NBOOT={NBOOT}; baseline ikut disimulasikan ulang tiap kali,")
print(f"  jadi ketidakpastian rung 1 ikut diperhitungkan.")

pd.DataFrame(rows).to_csv(RES/"17_power_curve.csv", index=False)
print(f"\n  -> results/17_power_curve.csv")
