"""Stage 27 -- check the 54 test-split descriptor JSONs before they are used for anything.

There is one way this campaign can fail silently, and it is the expensive one. vast_worker.sh
accepts an mdrun as finished when `prod.gro` or `prod.xtc` exists, which is true of a
trajectory that stopped early as well as one that ran to the end. A short trajectory extracts
perfectly happily, and the JSON it produces then blocks the real run forever, because the
worker skips any run that already has a result. The descriptors would simply be computed over
a shorter window than the 217 development runs they exist to be compared against.

The check is exact rather than approximate: a 25 ns production run with the project's output
interval gives 2001 frames and 20.01 ns analysed after the 5 ns equilibration discard, and
173 of the development runs carry precisely those two numbers. Anything else is a short run,
not a rounding difference. (The development set also holds 39 runs at 4501 frames and 17 at
~3600 from the campaign's original 50 ns design; the test split is 25 ns throughout, so only
2001 is admissible here.)

Usage: python scripts/verify_test_json.py [--dir results/16_rung4_per_run]
"""
import argparse, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
N_FRAMES, NS_ANALYSED = 2001, 20.01

ap = argparse.ArgumentParser()
ap.add_argument("--manifest", default="results/14m_md_runs_test.csv")
ap.add_argument("--dir", default="results/16_rung4_per_run",
                help="direktori JSON hasil, setelah pull_vast.sh")
a = ap.parse_args()

runs = pd.read_csv(ROOT/a.manifest).run.tolist()
d = ROOT/a.dir
missing, short, ok = [], [], []
for r in runs:
    f = d/f"{r}.json"
    if not f.exists():
        missing.append(r); continue
    j = json.loads(f.read_text())
    nf, ns = j.get("n_frames"), round(float(j.get("ns_analysed", 0)), 2)
    if nf != N_FRAMES or ns != NS_ANALYSED:
        short.append((r, nf, ns))
    else:
        ok.append(r)

print(f"manifest     : {len(runs)} run")
print(f"lengkap      : {len(ok)}  (n_frames={N_FRAMES}, ns_analysed={NS_ANALYSED})")
print(f"belum ada    : {len(missing)}" + (f"  {missing[:6]}" if missing else ""))
print(f"TIDAK SESUAI : {len(short)}")
for r, nf, ns in short:
    print(f"   {r}: n_frames={nf} ns_analysed={ns}")
if short:
    print("\nJSON ini dihitung dari trayektori yang tidak selesai dan TIDAK boleh dipakai.")
    print("Hapus JSON-nya di instance dan di results/, lalu jalankan ulang worker:")
    for r, _, _ in short:
        print(f"   rm -f {a.dir}/{r}.json   # dan /root/qsar/out/{r}.json di instance")
    raise SystemExit(1)
if missing:
    raise SystemExit(f"\n{len(missing)} run belum selesai; kampanye masih berjalan.")
print("\nSemua JSON lengkap dan sepanjang jendela yang sama dengan 217 run dev.")
