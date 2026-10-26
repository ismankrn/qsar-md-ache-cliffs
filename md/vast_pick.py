"""Rank Vast.ai offers for this project's MD workload, against criteria that were measured
rather than assumed.

Offers are rented within minutes, so a list of ids goes stale faster than it can be acted on.
This queries the public search API and re-ranks live, which is the only form the answer keeps.

WHAT THE RANKING IS BUILT ON
----------------------------
Measured on the three RTX 3090 instances run on 2026-10-23, on the real production system
(96,260 atoms, GPU-resident mdrun, NTOMP=4):

  * 273 ns/day, identical on all three, so it is the configuration and not one lucky card
  * 2.11 h per 25 ns run end to end, including equilibration and descriptor extraction
  * 8.2-8.9 GB of disk in use with everything installed, one trajectory in flight included

Only the RTX 3090 number is a measurement. Every other entry in SPEED below is an estimate
from published GROMACS scaling for a system of this size, and is printed with a ~ to say so.
If one is wrong the calendar moves; the money barely does, because cost per run turns out to
be nearly flat across GPU classes -- a cheap slow card simply spends more hours.

WHAT IS DELIBERATELY NOT A CRITERION
------------------------------------
  * dlperf: a deep-learning benchmark. It predicts nothing for GROMACS.
  * disk bandwidth: measured write demand is ~0.25 MB/s. The slowest offer on the market is
    three orders of magnitude above that. Over-weighting this once sent a recommendation
    badly wrong and had to be retracted.
  * PCIe bandwidth: on the campus GPU node an x4 card beat an x8 card by 5%, i.e. noise.

Usage:  python md/vast_pick.py [-n 6] [--max-hours 2.45] [--reliability 0.98]
"""
import argparse, json, urllib.parse, urllib.request

API = "https://console.vast.ai/api/v0/bundles/"
H_3090 = 2.11          # measured hours per 25 ns run, end to end
DISK_GB = 32           # measured peak use is 8.2-8.9 GB; see module docstring

# Relative GROMACS throughput, ~96k atoms, GPU-resident. 1.00 = RTX 3090 = measured.
SPEED = {
    "RTX 4090": 1.45, "RTX 4080S": 1.25, "RTX 4080": 1.20, "RTX 3090 Ti": 1.08,
    "RTX 3090": 1.00, "RTX A6000": 0.95, "RTX 3080 Ti": 0.95, "A40": 0.90,
    "RTX 4070S Ti": 0.90, "RTX 3080": 0.88, "RTX 4070 Ti": 0.85, "RTX 4070S": 0.80,
    "RTX 4070": 0.72, "Tesla V100S": 0.72, "Tesla V100": 0.68, "Titan RTX": 0.62,
    "RTX 2080 Ti": 0.60, "RTX 4060 Ti": 0.60, "Titan V": 0.60, "A10": 0.58,
    "Q RTX 8000": 0.58, "RTX 3070": 0.58, "Q RTX 6000": 0.55, "RTX 4060": 0.50,
    "RTX A4000": 0.50, "L4": 0.45, "RTX 3060": 0.42, "Tesla T4": 0.28,
}
# Architectures the conda-forge GROMACS 2023 CUDA build predates. A 2023 build carries no
# kernels for Blackwell (sm_120) and nothing here has tested whether it falls back cleanly,
# so these are refused rather than ranked low. Mining cards (CMP) have cut-down compute and
# erratic driver support; laptop parts throttle over runs that last hours.
AVOID = ("RTX 5090", "RTX 5080", "RTX 5070", "B200", "CMP ", "laptop")

ap = argparse.ArgumentParser()
ap.add_argument("-n", type=int, default=6, help="berapa instance yang akan disewa")
ap.add_argument("--runs", type=int, default=49, help="run yang tersisa")
ap.add_argument("--reliability", type=float, default=0.98)
ap.add_argument("--max-hours", type=float, default=2.45,
                help="batas jam/run; armada harus seragam karena slice dibagi rata")
ap.add_argument("--min-cores", type=float, default=8.0,
                help="core efektif; NTOMP=4 dan pair search jalan di CPU, 4 adalah batas mati")
a = ap.parse_args()

q = {"verified": {"eq": True}, "rentable": {"eq": True}, "num_gpus": {"eq": 1},
     "disk_space": {"gte": DISK_GB}, "cuda_max_good": {"gte": 12.4},
     "reliability2": {"gte": a.reliability}, "direct_port_count": {"gte": 8},
     "order": [["dph_total", "asc"]], "limit": 500, "type": "on-demand"}
url = API + "?q=" + urllib.parse.quote(json.dumps(q))
with urllib.request.urlopen(url, timeout=60) as r:
    offers = json.load(r)["offers"]

rows = []
for x in offers:
    g = x["gpu_name"]
    if any(b.lower() in g.lower() for b in AVOID):  continue
    if g not in SPEED:                              continue
    if x["cpu_cores_effective"] < a.min_cores:      continue
    if x["duration"] / 86400 < 7:                   continue
    x["hourly"] = x["dph_total"] + DISK_GB * x["storage_cost"] / 730
    x["h_run"] = H_3090 / SPEED[g]
    x["per_run"] = x["hourly"] * x["h_run"]
    x["measured"] = (g == "RTX 3090")
    rows.append(x)

band = sorted((x for x in rows if x["h_run"] <= a.max_hours), key=lambda x: x["hourly"])

# One offer per host: a fleet sharing a host shares its outages, and the whole point of
# spreading the manifest is that no single failure takes more than one slice with it.
fleet, seen = [], set()
for x in band:
    if x["host_id"] in seen: continue
    seen.add(x["host_id"]); fleet.append(x)
    if len(fleet) == a.n: break

print(f"{len(offers)} penawaran dari API, {len(rows)} lolos syarat keras, "
      f"{len(band)} di pita <= {a.max_hours} jam/run\n")
hdr = (f"{'gpu':<14} {'offer id':>9} {'$/jam':>6} {'jam/run':>8} {'$/run':>6} "
       f"{'rely':>6} {'core':>5} {'GHz':>4}  lokasi")
print(hdr); print("-" * len(hdr))
for x in band[:20]:
    mark = " " if x["measured"] else "~"
    print(f"{x['gpu_name']:<14} {x['id']:>9} {x['hourly']:>6.3f} {mark}{x['h_run']:>7.2f} "
          f"{x['per_run']:>6.3f} {x['reliability2']*100:>5.1f}% "
          f"{x['cpu_cores_effective']:>5.1f} {x['cpu_ghz']:>4.1f}  {x['geolocation']}")
print("\n~ = jam/run diestimasi dari penskalaan GROMACS, bukan diukur di proyek ini.")
print("  Hanya RTX 3090 yang terukur (273 ns/hari, 2.11 jam/run).")

if len(fleet) < a.n:
    print(f"\nHanya {len(fleet)} host berbeda tersedia di pita ini, diminta {a.n}.")
    print("Longgarkan --max-hours atau --reliability, lalu jalankan lagi.")
else:
    per = a.runs / a.n
    cost = sum(per * x["h_run"] * x["hourly"] for x in fleet)
    cal = max(per * x["h_run"] for x in fleet)
    print(f"\n=== pilihan {a.n} instance, host berbeda semua ===")
    for i, x in enumerate(fleet, 1):
        print(f"  {i}. {x['gpu_name']:<13} {x['id']:>9}  ${x['hourly']:.3f}/jam  "
              f"{x['h_run']:.2f} jam/run  {x['reliability2']*100:.1f}%  {x['geolocation']}")
    print(f"\n  {a.runs} run, {per:.1f} per instance")
    print(f"  total   : ${cost:.2f}")
    print(f"  kalender: {cal:.1f} jam, ditentukan instance paling lambat")
