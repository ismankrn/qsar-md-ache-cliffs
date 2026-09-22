#!/usr/bin/env python
"""Stage 07 — production docking of the full ligand set into 4EY7.

Every ligand is docked with three independent seeds, because Vina is
stochastic: during protocol validation galantamine gave 5.89 A on one seed
and 1.61 A on four others with identical settings. A single run cannot be
trusted.

The top three poses per run are kept rather than only the best, because
validation also showed the correct pose is sometimes ranked below first.
Carrying several poses into the next stage lets a short MD pick the stable one.

Resumable by design: a run whose output already exists is skipped, and
scores are appended to the CSV as they are produced, so an interrupted job
loses at most one docking.

    python run_docking.py [--seeds 3] [--exhaustiveness 32] [--limit N]
"""
from __future__ import annotations
import argparse, csv, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from vina import Vina
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES, DOCK = ROOT / "results", ROOT / "docking"
POSES = DOCK / "poses"
TARGET = "CHEMBL220"
RECEPTOR = "4ey7"
BOX = 30.0
N_KEEP = 3

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--exhaustiveness", type=int, default=32)
ap.add_argument("--limit", type=int, default=0, help="dock only the first N ligands")
a = ap.parse_args()

POSES.mkdir(parents=True, exist_ok=True)
out_csv = RES / "07_docking_scores.csv"

# ---- ligand list: study set + continuous-delta validation set ----------
study = pd.read_csv(RES / f"02_ligand_set_{TARGET}.csv")
valid = pd.read_csv(RES / f"04_validation_ligands_{TARGET}.csv")
cols = ["molecule_chembl_id", "smiles_md", "pchembl"]
lig = (pd.concat([study[cols].assign(origin="study"),
                  valid[cols].assign(origin="validation")])
         .drop_duplicates("molecule_chembl_id")
         .reset_index(drop=True))
if a.limit:
    lig = lig.head(a.limit)
print(f"ligands to dock : {len(lig)}", flush=True)
print(f"seeds per ligand: {a.seeds}   ->  {len(lig)*a.seeds} docking runs", flush=True)

# ---- receptor, prepared once and reused ------------------------------
center = np.load(DOCK / f"{RECEPTOR}_box_center.npy")
print(f"receptor {RECEPTOR}, box {BOX} A at {np.round(center,2)}, "
      f"exhaustiveness {a.exhaustiveness}\n", flush=True)

done = set()
if out_csv.exists():
    prev = pd.read_csv(out_csv)
    done = {(r.molecule_chembl_id, int(r.seed)) for r in prev.itertuples()}
    print(f"resuming: {len(done)} runs already recorded\n", flush=True)

fh = open(out_csv, "a", newline="")
w = csv.writer(fh)
if not done:
    w.writerow(["molecule_chembl_id", "origin", "pchembl", "seed",
                "score_1", "score_2", "score_3", "n_poses", "seconds"])
    fh.flush()

t_start = time.time()
n_run = n_skip = n_fail = 0
for i, r in lig.iterrows():
    # ligand PDBQT, prepared once per ligand and reused across seeds
    lig_pdbqt = POSES / f"{r.molecule_chembl_id}_lig.pdbqt"
    if not lig_pdbqt.exists():
        try:
            m = Chem.AddHs(Chem.MolFromSmiles(r.smiles_md))
            p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
            if AllChem.EmbedMolecule(m, p) != 0:
                raise RuntimeError("embedding failed")
            AllChem.MMFFOptimizeMolecule(m, maxIters=500)
            s, ok, err = PDBQTWriterLegacy.write_string(
                MoleculePreparation().prepare(m)[0])
            if not ok:
                raise RuntimeError(f"meeko: {err}")
            lig_pdbqt.write_text(s)
        except Exception as e:
            print(f"  PREP FAIL {r.molecule_chembl_id}: {type(e).__name__}: "
                  f"{str(e)[:60]}", flush=True)
            n_fail += 1
            continue

    for seed in range(1, a.seeds + 1):
        out = POSES / f"{r.molecule_chembl_id}_seed{seed}.pdbqt"
        if (r.molecule_chembl_id, seed) in done and out.exists():
            n_skip += 1
            continue
        t0 = time.time()
        try:
            v = Vina(sf_name="vina", cpu=0, verbosity=0, seed=seed)
            v.set_receptor(str(DOCK / f"{RECEPTOR}_rec.pdbqt"))
            v.compute_vina_maps(center=[float(c) for c in center],
                                box_size=[BOX] * 3)
            v.set_ligand_from_file(str(lig_pdbqt))
            v.dock(exhaustiveness=a.exhaustiveness, n_poses=N_KEEP)
            v.write_poses(str(out), n_poses=N_KEEP, overwrite=True)
            e = v.energies(n_poses=N_KEEP)
            sc = [float(e[k][0]) if k < len(e) else None for k in range(3)]
            w.writerow([r.molecule_chembl_id, r.origin, r.pchembl, seed,
                        sc[0], sc[1], sc[2], len(e), round(time.time()-t0, 1)])
            fh.flush()
            n_run += 1
        except Exception as ex:
            print(f"  DOCK FAIL {r.molecule_chembl_id} seed {seed}: "
                  f"{type(ex).__name__}: {str(ex)[:60]}", flush=True)
            n_fail += 1

    if (i + 1) % 10 == 0 or i == len(lig) - 1:
        el = time.time() - t_start
        rate = n_run / el * 3600 if n_run else 0
        remain = (len(lig) - i - 1) * a.seeds
        eta = remain / rate if rate else 0
        print(f"  [{i+1:>3}/{len(lig)}] {n_run} run, {n_skip} skip, "
              f"{n_fail} fail | {rate:.0f} run/jam | ETA {eta:.1f} jam",
              flush=True)

fh.close()
print(f"\nselesai: {n_run} run baru, {n_skip} dilewati, {n_fail} gagal, "
      f"{(time.time()-t_start)/3600:.2f} jam")
print(f"skor    : {out_csv}")
print(f"pose    : {POSES}")
