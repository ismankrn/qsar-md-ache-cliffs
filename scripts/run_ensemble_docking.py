#!/usr/bin/env python
"""Stage 11c — dock every ligand into all eight ensemble receptors.

Single docking into 4EY7 failed for ligands whose binding mode differs from
donepezil's, and the cause was measured: the Tyr337 side chain moves 2.53 A
between structures and blocks the true pose. This asks whether that is why
rungs 2 and 3 found nothing, by giving every ligand eight receptor
conformations to choose from.

One seed per receptor rather than three. Stage 07 measured the seed-to-seed
spread at 0.14 kcal/mol median, while the ensemble members differ by 1.51 A
in gorge side-chain position -- the conformational variation is far larger
than the sampling noise, so the budget is better spent on receptors.

Ligand PDBQT files prepared in stage 07 are reused. Resumable: existing
outputs are skipped and scores are flushed as they are produced.
"""
from __future__ import annotations
import argparse, csv, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from vina import Vina

ROOT = Path(__file__).resolve().parent.parent
RES, DOCK = ROOT/"results", ROOT/"docking"
POSES, ENS_POSES = DOCK/"poses", DOCK/"ens_poses"
BOX, N_KEEP = 30.0, 3

ap = argparse.ArgumentParser()
ap.add_argument("--exhaustiveness", type=int, default=32)
ap.add_argument("--limit", type=int, default=0)
a = ap.parse_args()

ENS_POSES.mkdir(parents=True, exist_ok=True)
receptors = (DOCK/"ensemble_ready.txt").read_text().split()
center = [float(c) for c in np.load(DOCK/f"4ey7_box_center.npy")]
ligands = sorted(pd.read_csv(RES/"07_docking_scores.csv").molecule_chembl_id.unique())
if a.limit:
    ligands = ligands[:a.limit]
print(f"ligands   : {len(ligands)}")
print(f"receptors : {len(receptors)}  ({' '.join(r.upper() for r in receptors)})")
print(f"runs      : {len(ligands)*len(receptors)}\n", flush=True)

out_csv = RES/"11_ensemble_docking.csv"
done = set()
if out_csv.exists():
    prev = pd.read_csv(out_csv)
    done = {(r.molecule_chembl_id, r.receptor) for r in prev.itertuples()}
    print(f"resuming: {len(done)} runs recorded\n", flush=True)

fh = open(out_csv, "a", newline="")
w = csv.writer(fh)
if not done:
    w.writerow(["molecule_chembl_id", "receptor", "score_1", "score_2",
                "score_3", "n_poses", "seconds"])
    fh.flush()

# maps are expensive to build, so the outer loop is over receptors
t0, n_run, n_skip, n_fail = time.time(), 0, 0, 0
for ri, rec in enumerate(receptors, 1):
    v = Vina(sf_name="vina", cpu=0, verbosity=0, seed=1)
    v.set_receptor(str(DOCK/f"ens_{rec}_rec.pdbqt"))
    v.compute_vina_maps(center=center, box_size=[BOX]*3)
    print(f"[{ri}/{len(receptors)}] {rec.upper()}", flush=True)
    for i, cid in enumerate(ligands, 1):
        if (cid, rec) in done:
            n_skip += 1; continue
        lig = POSES/f"{cid}_lig.pdbqt"
        if not lig.exists():
            n_fail += 1; continue
        out = ENS_POSES/f"{cid}_{rec}.pdbqt"
        t1 = time.time()
        try:
            v.set_ligand_from_file(str(lig))
            v.dock(exhaustiveness=a.exhaustiveness, n_poses=N_KEEP)
            v.write_poses(str(out), n_poses=N_KEEP, overwrite=True)
            e = v.energies(n_poses=N_KEEP)
            sc = [float(e[k][0]) if k < len(e) else None for k in range(3)]
            w.writerow([cid, rec, sc[0], sc[1], sc[2], len(e),
                        round(time.time()-t1, 1)])
            fh.flush(); n_run += 1
        except Exception as ex:
            n_fail += 1
            print(f"    FAIL {cid}: {type(ex).__name__}: {str(ex)[:60]}", flush=True)
        if i % 100 == 0:
            el = (time.time()-t0)/3600
            total = len(ligands)*len(receptors)
            frac = (n_run+n_skip)/total
            print(f"    {i}/{len(ligands)}  overall {100*frac:.1f}%  "
                  f"{el:.2f} h  ETA {el/max(frac,1e-9)*(1-frac):.1f} h", flush=True)

fh.close()
print(f"\ndone: {n_run} new, {n_skip} skipped, {n_fail} failed, "
      f"{(time.time()-t0)/3600:.2f} h")
d = pd.read_csv(out_csv)
print(f"  rows: {len(d)}   ligands: {d.molecule_chembl_id.nunique()}"
      f"   receptors: {d.receptor.nunique()}")
best = d.groupby("molecule_chembl_id").score_1.min()
single = d[d.receptor == "4ey7"].set_index("molecule_chembl_id").score_1
gain = (single - best.reindex(single.index)).dropna()
print(f"  improvement from the ensemble over 4EY7 alone: "
      f"median {gain.median():.2f}, max {gain.max():.2f} kcal/mol")
