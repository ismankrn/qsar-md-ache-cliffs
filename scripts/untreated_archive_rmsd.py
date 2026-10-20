#!/usr/bin/env python
"""What an RMSD computed straight off a deposited archive returns.

The archives in md/reduced are wrapped molecule-wise and the capped receptor is three
separate chains, so a reader who opens one and computes an RMSD against its first frame
without unwrapping or superposing gets a number dominated by the box, not by the
structure. The data-availability statement warns about this and needs a figure behind
it rather than an impression.

Measured exactly as the naive reader would: RMSD of the C-alpha atoms, and of the
ligand heavy atoms, against frame 0 of the same archive, with no periodic-boundary
treatment and no superposition. Compare with results/25_md_stability.csv, which
reports the same quantities after treatment.

Usage:  python scripts/untreated_archive_rmsd.py [--limit N]
"""
from __future__ import annotations
import argparse, glob, os, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda

ROOT = Path(__file__).resolve().parent.parent
RED, RES = ROOT/"md"/"reduced", ROOT/"results"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    runs = sorted(os.path.basename(x)[:-4] for x in glob.glob(str(RED/"*.xtc"))
                  if not os.path.basename(x).startswith("#"))
    if a.limit:
        runs = runs[:a.limit]
    rows = []
    for i, r in enumerate(runs, 1):
        try:
            u = mda.Universe(str(RED/f"{r}.gro"), str(RED/f"{r}.xtc"))
            ca = u.select_atoms("protein and name CA")
            lig = u.select_atoms("resname UNK and not name H*")
            u.trajectory[0]
            c0, l0 = ca.positions.copy(), lig.positions.copy()
            cv, lv = [], []
            for _ in u.trajectory:
                cv.append(np.sqrt(((ca.positions - c0)**2).sum(1).mean()))
                lv.append(np.sqrt(((lig.positions - l0)**2).sum(1).mean()))
            rows.append(dict(run=r, n_frames=len(cv),
                             ca_rmsd_mean=float(np.mean(cv)), ca_rmsd_max=float(np.max(cv)),
                             lig_rmsd_mean=float(np.mean(lv)), lig_rmsd_max=float(np.max(lv))))
        except Exception as e:
            rows.append(dict(run=r, n_frames=0, error=str(e)[:80]))
        if i % 25 == 0:
            print(f"  {i}/{len(runs)}", flush=True)
    d = pd.DataFrame(rows)
    ok = d[d.n_frames > 0]
    d.to_csv(RES/"25b_untreated_archive_rmsd.csv", index=False)
    print(f"\n  runs read {len(ok)} of {len(d)}")
    print(f"  receptor CA RMSD, median of per-run means {ok.ca_rmsd_mean.median():.1f} A, "
          f"largest single frame {ok.ca_rmsd_max.max():.1f} A")
    print(f"  ligand      RMSD, median of per-run means {ok.lig_rmsd_mean.median():.1f} A, "
          f"largest single frame {ok.lig_rmsd_max.max():.1f} A")
    print(f"  -> results/25b_untreated_archive_rmsd.csv")


if __name__ == "__main__":
    main()
