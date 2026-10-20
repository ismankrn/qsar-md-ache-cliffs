#!/usr/bin/env python
"""The replicate test applied to ProLIF typed occupancy, for comparison with the
distance-only occ_cas and occ_pas of stage 18.

Section 3.5 contrasts the two measures on reproducibility rather than on chemistry, and
that contrast needs both sides computed the same way. This mirrors replicate_pilot.py
exactly: sd_rep is the pooled within-ligand standard deviation, the square root of the
mean within-ligand variance rather than the mean of the within-ligand standard
deviations, and sd_lig is the standard deviation of the per-ligand means taken over
every ligand, not only the replicated ones.

Typed occupancy at a residue is the maximum over interaction types, which is what
Figure 11(b) plots and what "engagement of any type" means.

Usage:  python scripts/typed_replicate_ratio.py
"""
from __future__ import annotations
import re, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
SITES = {86: "Trp86", 286: "Trp286"}


def main() -> None:
    d = pd.read_csv(RES / "26_dynamic_interactions.csv")
    tot = (d[d.resid_4ey7.isin(SITES)]
           .groupby(["run", "resid_4ey7"]).occupancy.max().unstack())
    tot["ligand"] = [re.sub(r"_r\d+$", "", r) for r in tot.index]

    rows = []
    for resid, name in SITES.items():
        va = tot[[resid, "ligand"]].dropna()
        reps = va.groupby("ligand").filter(lambda g: len(g) >= 2)
        var_within = reps.groupby("ligand")[resid].var(ddof=1).dropna()
        sd_rep = float(np.sqrt(var_within.mean()))
        sd_lig = float(va.groupby("ligand")[resid].mean().std(ddof=1))
        rows.append(dict(site=name, resid=resid, measure="ProLIF typed (any type)",
                         n_ligands=va.ligand.nunique(),
                         n_replicated=var_within.size,
                         sd_rep=sd_rep, sd_lig=sd_lig,
                         ratio=sd_lig / sd_rep if sd_rep else np.inf))

    pilot = pd.read_csv(RES / "18_replicate_pilot.csv")
    w = pilot[pilot.window == 20]
    for resid, name, col in ((86, "Trp86", "occ_cas"), (286, "Trp286", "occ_pas")):
        r = w[w.descriptor == col]
        if len(r):
            r = r.iloc[0]
            rows.append(dict(site=name, resid=resid, measure=f"distance-only {col}",
                             n_ligands=int(r.n_lig_all), n_replicated=int(r.n_lig_rep),
                             sd_rep=float(r.sd_rep), sd_lig=float(r.sd_lig),
                             ratio=float(r.ratio)))

    out = pd.DataFrame(rows)
    out.to_csv(RES / "26c_typed_replicate_ratio.csv", index=False)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\n  -> results/26c_typed_replicate_ratio.csv")


if __name__ == "__main__":
    main()
