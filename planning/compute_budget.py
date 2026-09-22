#!/usr/bin/env python3
"""Compute + cost estimator for the activity-cliff QSAR/MD study.

All throughput figures are ESTIMATES for a ~50,000-atom solvated
protein-ligand system on a single GPU, GROMACS/OpenMM in GPU-resident mode.
Re-run with your own numbers once you have a real benchmark.
"""

# ---- study design -------------------------------------------------
N_LIGANDS      = 60     # compounds in the QSAR series
N_REPLICAS     = 3      # independent velocity seeds per ligand
NS_PER_REPLICA = 30     # production nanoseconds per replica
APO_NS         = 300    # apo-protein control simulations
EQUIL_OVERHEAD = 0.15   # minimisation + NVT/NPT, as fraction of production

# ---- throughput (ns/day, one GPU) ---------------------------------
THROUGHPUT = {
    "RTX 4090, 2 fs":        200,
    "RTX 4090, 4 fs (HMR)":  350,
    "RTX 3090, 2 fs":        120,
    "A100 40GB, 2 fs":       210,
    "AWS A10G (g5), 2 fs":   100,
    "Tesla T4 (Kaggle free)": 45,
}

# ---- price per GPU-hour (USD, approx. late 2026) -------------------
PRICES = {
    "Vast.ai 4090 (interruptible)": 0.35,
    "RunPod 4090 Community":        0.34,
    "RunPod 4090 Secure":           0.74,
    "Lambda A100 80GB":             2.06,
    "AWS g5.xlarge on-demand":      1.01,
    "AWS g5.xlarge spot":           0.35,
}

prod_ns  = N_LIGANDS * N_REPLICAS * NS_PER_REPLICA + APO_NS
total_ns = prod_ns * (1 + EQUIL_OVERHEAD)

print(f"Production:        {prod_ns:,.0f} ns")
print(f"Total incl. equil: {total_ns:,.0f} ns  ({total_ns/1000:.1f} us)\n")

print(f"{'hardware / timestep':<24} {'GPU-days':>9} {'GPU-hours':>10}")
print("-" * 46)
hours = {}
for name, nsday in THROUGHPUT.items():
    days = total_ns / nsday
    hours[name] = days * 24
    print(f"{name:<24} {days:>9.1f} {days*24:>10.0f}")

print(f"\nCost at each provider, assuming RTX 4090 @ 2 fs "
      f"({hours['RTX 4090, 2 fs']:.0f} GPU-hours):\n")
h = hours["RTX 4090, 2 fs"]
for prov, rate in sorted(PRICES.items(), key=lambda kv: kv[1]):
    print(f"  {prov:<32} ${h * rate:>8,.0f}")

h_hmr = hours["RTX 4090, 4 fs (HMR)"]
best  = min(PRICES.values())
print(f"\nWith hydrogen mass repartitioning at 4 fs on the cheapest "
      f"option: ${h_hmr * best:,.0f}")
