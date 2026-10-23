"""Stage 14 — lay out the 379 production runs and give each one a recorded velocity seed.

One directory per trajectory, not per system: 217 systems yield 379 runs because the 81
cliff-pair ligands get three replicates each. The system files are symlinked, not copied,
so the replicates provably share one structure and 217 x ~9 MB is not written four times.

Each replicate is equilibrated independently, with its own velocity seed at the NVT
velocity-generation step, rather than branching three production runs off one equilibrated
state. Equilibration is cheap (200 ps against 50 ns of production), and independent
equilibration makes the replicates genuinely independent samples -- which is the whole
point of running them.

Seeds are derived deterministically from the ligand id, the replicate index and one base
seed, then written into the manifest. A random `gen-seed = -1` would make the
reproducibility test unreproducible, which defeats it.
"""
import argparse
import hashlib, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES, MD = ROOT/"results", ROOT/"md"
SYS, RUNS, MDP = MD/"systems", MD/"runs", MD/"mdp"
BASE_SEED = 20260902          # the project's frozen-split seed, reused for traceability

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", default="13_md_ligands.csv",
                help="ligand manifest under results/")
ap.add_argument("--out", default="14_md_runs.csv",
                help="run list to write under results/")
ap.add_argument("--slices", type=int, default=0,
                help="also write md/slice_NN.txt, one per instance")
a = ap.parse_args()
manifest = pd.read_csv(RES/a.manifest)
RUNS.mkdir(parents=True, exist_ok=True)

def seed_for(cid, rep):
    h = hashlib.sha256(f"{cid}|r{rep}|{BASE_SEED}".encode()).hexdigest()
    return int(h[:8], 16) % (2**31 - 1)

rows, missing = [], []
for _, r in manifest.iterrows():
    cid, nrep = r.molecule_chembl_id, int(r.n_replicates)
    src = SYS/cid
    if not (src/"system.top").exists() or not (src/"system.gro").exists():
        missing.append(cid); continue
    for rep in range(1, nrep+1):
        name = f"{cid}_r{rep}"
        d = RUNS/name
        d.mkdir(exist_ok=True)
        for f in ("system.top", "system.gro"):
            link = d/f
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(Path("../../systems")/cid/f)
        seed = seed_for(cid, rep)
        # only nvt.mdp differs per run; the rest are shared and symlinked
        nvt = (MDP/"nvt.mdp").read_text().replace("gen-seed                = -1",
                                                 f"gen-seed                = {seed}")
        (d/"nvt.mdp").write_text(nvt)
        for f in ("em.mdp", "npt.mdp", "prod.mdp"):
            link = d/f
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(Path("../../mdp")/f)   # from md/runs/<name>/, mdp is two levels up
        rows.append({"run": name, "molecule_chembl_id": cid, "replicate": rep,
                     "n_replicates": nrep, "gen_seed": seed,
                     "in_cliff_pair": bool(r.in_cliff_pair)})

df = pd.DataFrame(rows)
df.to_csv(RES/a.out, index=False)
print(f"sistem hilang        : {len(missing)}" + (f"  {missing[:5]}" if missing else ""))
print(f"direktori run dibuat : {len(df)}")
if df.empty:
    raise SystemExit(
        "\nTidak ada sistem yang ditemukan. Jalankan build_md_systems.py lebih dulu,\n"
        "dari lingkungan qsar-md dan bukan dari base:\n"
        f"    conda run -n qsar-md python scripts/build_md_systems.py --manifest {a.manifest}")
print(f"  replika 1 saja     : {int((df.n_replicates==1).sum())}")
print(f"  bagian dari 3      : {int((df.n_replicates==3).sum())}")
print(f"ligan unik           : {df.molecule_chembl_id.nunique()}")

# 2. splitting the list into per-instance slices used to be a shell step, with
#    `split -n l/N`, which only GNU split understands; macOS rejects it outright.
if a.slices:
    import math
    per = math.ceil(len(df) / a.slices)
    for i in range(a.slices):
        chunk = df.run.iloc[i*per:(i+1)*per]
        if chunk.empty:
            continue
        f = RUNS.parent/f"slice_{i+1:02d}.txt"
        f.write_text("\n".join(chunk) + "\n")
        print(f"  {f.relative_to(ROOT)}: {len(chunk)} run")
print(f"seed unik            : {df.gen_seed.nunique()} / {len(df)}")
ns = int((ROOT/"md"/"NSTEPS_PROD").read_text().strip()) * 0.002 / 1000 \
     if (ROOT/"md"/"NSTEPS_PROD").exists() else 50
print(f"ns total @ {ns:g} ns     : {len(df)*ns:,.0f}")
