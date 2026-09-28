"""Stage 17 — can this study detect a rung-4 improvement at all, and what is the ceiling?

Two questions that decide whether the full MD design is worth buying, both answerable from
data already in hand:

1. NOISE CEILING. pChEMBL aggregates assays from different labs. If a compound's potency is
   measured with standard deviation sigma, then Δp carries noise from both members of the
   pair, and no descriptor can explain that part. The ceiling is the R² a perfect predictor
   of the TRUE Δp would score against the OBSERVED Δp.

2. POWER. The bootstrap in stage 15 gave a 95% CI about ±0.12 wide on R² over 145 pairs.
   This asks it directly: for a descriptor whose true R² is X, how often would this study
   call it better than the rung-1 baseline of 0.200?

A null result means "no effect large enough for 145 pairs to see", not "no effect". Which of
those two sentences the paper gets to write is decided here, before any money is spent.
"""
from pathlib import Path
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
SEED = 20260902
rng = np.random.default_rng(SEED)

cmp_ = pd.read_csv(RES/"01_compounds_CHEMBL220.csv")
pairs = pd.read_csv(RES/"03_splits_pairs_CHEMBL220.csv")
dev = pairs[pairs.split == "dev"].reset_index(drop=True)

print("=" * 72)
print("1. LANGIT-LANGIT DERAU PENGUKURAN")
print("=" * 72)
rep = cmp_[(cmp_.n_meas > 1) & cmp_.sd.notna()]
print(f"  senyawa dengan >1 pengukuran : {len(rep)} dari {len(cmp_)} "
      f"({100*len(rep)/len(cmp_):.1f}%)")
print(f"  sd antar-assay               : median {rep.sd.median():.3f}, "
      f"rata2 {rep.sd.mean():.3f}, p90 {rep.sd.quantile(.9):.3f} satuan log")
sigma = float(rep.sd.median())

# standard error of each compound's reported pchembl = sd / sqrt(n_meas);
# compounds measured once get the population sigma
info = cmp_.set_index("molecule_chembl_id")
def se_of(cid):
    if cid not in info.index: return sigma
    r = info.loc[cid]
    if pd.notna(r.sd) and r.n_meas and r.n_meas > 1:
        return float(r.sd) / np.sqrt(float(r.n_meas))
    return sigma

se_a = np.array([se_of(c) for c in dev.ligand_a])
se_b = np.array([se_of(c) for c in dev.ligand_b])
noise_var = (se_a**2 + se_b**2).mean()
obs_var = float(np.var(dev.delta, ddof=1))
ceiling = max(0.0, 1.0 - noise_var/obs_var)
print(f"\n  varians Δp teramati          : {obs_var:.3f}")
print(f"  varians derau pengukuran     : {noise_var:.3f}")
print(f"  LANGIT-LANGIT R²             : {ceiling:.3f}")
print(f"    (prediktor sempurna atas Δp SEBENARNYA hanya mencapai ini)")
print(f"\n  rung 1 (2D) R² = 0.200  ->  {100*0.200/ceiling:.0f}% dari yang mungkin")
# Only 6% of compounds carry a repeat measurement, and the spread is skewed (median 0.064,
# mean 0.184, p90 0.560), so the ceiling depends on which summary is taken. Report all three
# rather than the flattering one.
print(f"\n  sensitivitas langit-langit terhadap pilihan sigma:")
for lbl, sg in (("median", rep.sd.median()), ("rata-rata", rep.sd.mean()),
                ("p90 (pesimistis)", rep.sd.quantile(.9))):
    nv = 2*float(sg)**2
    print(f"    sigma={float(sg):.3f} ({lbl:16s}) -> langit-langit R² {max(0,1-nv/obs_var):.3f}")

print()
print("=" * 72)
print("2. DAYA: bisakah 145 pasangan membedakan rung 4 dari rung 1?")
print("=" * 72)
n_pairs = len(dev)
y = dev.delta.to_numpy()
y = np.concatenate([y, -y])                      # framing B antisimetris
pair_id = np.concatenate([np.arange(n_pairs)]*2)
rows_of = {p: np.where(pair_id == p)[0] for p in range(n_pairs)}
BASE, NBOOT, NSIM = 0.200, 1500, 400

def r2_of(yt, yp):
    ss = ((yt-yp)**2).sum(); st = ((yt-yt.mean())**2).sum()
    return 1 - ss/st

print(f"  {'R² sebenarnya':>14s}{'R² terukur':>13s}{'lebar CI 95%':>15s}{'daya vs 0.200':>16s}")
for true_r2 in (0.20, 0.30, 0.40, 0.50, 0.60):
    detect, widths, ests = 0, [], []
    for _ in range(NSIM):
        # The synthetic predictor must be CALIBRATED, i.e. shrunk toward the mean the way a
        # fitted model is. An uncalibrated signal with correlation sqrt(r2) scores
        # R² = 2*sqrt(r2) - 1, which is -0.106 at r2 = 0.20 -- a systematic penalty, not
        # sampling error, and it made the study look far weaker than it is.
        #   pred = r2*ys + sqrt(r2*(1-r2))*eps  ->  var = r2, cov = r2, R² = r2
        ys = (y - y.mean())/y.std()
        pred = true_r2*ys + np.sqrt(true_r2*(1-true_r2))*rng.standard_normal(len(ys))
        pred = pred*y.std() + y.mean()
        ests.append(r2_of(y, pred))
        b = np.empty(NBOOT)
        for k in range(NBOOT):
            samp = rng.integers(0, n_pairs, n_pairs)
            ii = np.concatenate([rows_of[p] for p in samp])
            b[k] = r2_of(y[ii], pred[ii])
        lo, hi = np.percentile(b, [2.5, 97.5])
        widths.append(hi-lo)
        if lo > BASE: detect += 1
    print(f"  {true_r2:>14.2f}{np.mean(ests):>13.3f}{np.mean(widths):>15.3f}"
          f"{100*detect/NSIM:>15.0f}%")

print(f"\n  (daya = berapa persen simulasi yang batas bawah CI-nya melewati 0.200)")
