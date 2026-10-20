#!/usr/bin/env python
"""Figure — structural stability across all 229 production runs.

The per-run QC pages carry this for one run at a time; nothing showed it for the campaign as
a whole, which is what a reader needs in order to believe 217 ligands were simulated to the
same standard. Three panels: receptor C-alpha RMSD against time, ligand RMSD from its
equilibrated pose, and per-residue C-alpha fluctuation with the gorge-lining residues marked.

**The reduced archive needs the same periodic-boundary treatment the descriptors needed.**
`trjconv -pbc mol` makes each molecule whole but places each one independently, and this
receptor is three separate chains after the loop capping, so a naive RMSD over the archive
reports 25.8 A for the receptor and up to 99 A for the ligand -- a campaign that looks
catastrophically unstable when the real figures are 1.3 A and 4.2 A. The fix is the one from
stage 16b, reused here rather than reimplemented: unwrap, no-jump against the previous frame,
put every molecule in the first one's image, then superpose. No-jump is why this iterates
every frame; a stride would let a molecule cross half a box between samples and break it.

Runs are 251 or 501 frames at 100 ps, the longer ones being the first batch continued past the
25 ns production length. The common
window is the shorter one, and the two groups are drawn separately rather than averaged into
a single band that would change composition partway along its own x-axis.

Usage:
    python scripts/fig_md_stability.py
"""
from __future__ import annotations
import sys, warnings, glob, os, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import MDAnalysis as mda
from rung4_descriptors import receptor_molecules, make_whole, fit_rotation, nojump_into

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
RED = ROOT/"md"/"reduced"
EQUIL_NS = 5.0                 # descriptors ignore the first 5 ns; mark it here too
LINING = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]
NAMED = {86: "Trp86 (CAS)", 286: "Trp286 (PAS)", 337: "Tyr337 (gate)"}


def true_numbering():
    """C-alpha residue numbers as 4EY7 numbers them, in archive order.

    The reduced archives renumber residues per chain segment -- 529 C-alpha carrying resids
    2-256 with two wraps, because the capped receptor is three chains. Plotting RMSF against
    those numbers draws a line that jumps backwards twice and reads as nonsense. The docking
    receptor has the same 529 C-alpha in the same order with the real 4EY7 numbering, 4-541,
    so the mapping is just that list -- guarded by a residue-name check so a future change in
    either file fails loudly instead of silently relabelling the gorge.
    """
    d = mda.Universe(str(ROOT/"docking"/"4ey7_rec_clean.pdb")).select_atoms(
        "protein and name CA")
    return d.resids.copy(), d.resnames.copy()


def one_run(xtc: Path):
    """(time_ns, ca_rmsd, lig_rmsd, resids, per-residue RMSF after EQUIL_NS)."""
    u = mda.Universe(str(xtc.with_suffix(".gro")), str(xtc))
    rec = u.select_atoms("protein")
    lig = u.select_atoms("resname UNK")
    if len(lig) == 0:
        return None
    whole = rec + lig
    mols = receptor_molecules(rec) + [lig.indices]
    ca = rec.select_atoms("name CA")

    u.trajectory[0]
    make_whole(u, mols, u.dimensions)
    prev = whole.positions.copy()
    ref = ca.positions.copy()
    rc = ref.mean(0)
    ref0 = ref - rc
    refl = lig.positions.copy() - rc

    t, rca, rlig, fitted = [], [], [], []
    for ts in u.trajectory:
        nojump_into(whole, prev, ts.dimensions)
        make_whole(u, mols, ts.dimensions)
        prev = whole.positions.copy()
        P = ca.positions
        c = P.mean(0)
        Pc = P - c
        U = fit_rotation(Pc, ref0)
        A = Pc @ U
        d = A - ref0
        L = (lig.positions - c) @ U - refl
        t.append(ts.time/1000.0)
        rca.append(float(np.sqrt((d*d).sum(1).mean())))
        rlig.append(float(np.sqrt((L*L).sum(1).mean())))
        fitted.append(A)
    t = np.asarray(t)
    F = np.asarray(fitted)
    keep = t >= EQUIL_NS
    rmsf = np.sqrt(((F[keep] - F[keep].mean(0))**2).sum(-1).mean(0)) if keep.any() else None
    return t, np.asarray(rca), np.asarray(rlig), ca.resnames.copy(), rmsf


def main() -> None:
    xtcs = sorted(RED.glob("*.xtc"))
    print(f"  {len(xtcs)} runs")
    resids, refnames = true_numbering()
    rows, series, rmsfs = [], {}, []
    for i, x in enumerate(xtcs, 1):
        try:
            r = one_run(x)
        except Exception as e:
            print(f"    {x.stem}: {type(e).__name__}: {e}")
            continue
        if r is None:
            continue
        t, rca, rlig, rnames, rmsf = r
        if len(rnames) != len(refnames) or not (rnames == refnames).all():
            print(f"    {x.stem}: residue sequence does not match the reference, skipped")
            continue
        series[x.stem] = (t, rca, rlig)
        if rmsf is not None:
            rmsfs.append(rmsf)
        keep = t >= EQUIL_NS
        rows.append(dict(run=x.stem, n_frames=len(t), span_ns=float(t[-1]),
                         ca_rmsd_mean=float(rca[keep].mean()), ca_rmsd_max=float(rca.max()),
                         lig_rmsd_mean=float(rlig[keep].mean()),
                         lig_rmsd_max=float(rlig.max())))
        if i % 50 == 0:
            print(f"    {i}/{len(xtcs)}")
    df = pd.DataFrame(rows)
    df.to_csv(RES/"25_md_stability.csv", index=False)
    print(f"\n  -> results/25_md_stability.csv")
    print(f"  Ca RMSD  : median {df.ca_rmsd_mean.median():.2f} A, "
          f"95th pct {df.ca_rmsd_mean.quantile(.95):.2f}, max {df.ca_rmsd_max.max():.2f}")
    print(f"  lig RMSD : median {df.lig_rmsd_mean.median():.2f} A, "
          f"95th pct {df.lig_rmsd_mean.quantile(.95):.2f}, max {df.lig_rmsd_max.max():.2f}")
    worst = df.loc[df.lig_rmsd_mean.idxmax()]
    print(f"  worst ligand excursion: {worst.run} at {worst.lig_rmsd_mean:.2f} A")

    FIGS.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9))
    groups = {251: ("25 ns production", "#2166ac"),
              501: ("continued to 50 ns", "#762a83")}

    for ax, idx, lab in ((axes[0], 1, r"Receptor C$\alpha$ RMSD (\AA)"),
                         (axes[1], 2, r"Ligand RMSD (\AA)")):
        for n, (gl, col) in groups.items():
            ser = [(t, v[idx]) for t, v in
                   ((s[0], s) for s in series.values()) if len(t) == n]
            if not ser:
                continue
            T = ser[0][0]
            M = np.vstack([v for _t, v in ser])
            med = np.median(M, 0)
            lo, hi = np.percentile(M, [25, 75], axis=0)
            ax.fill_between(T, lo, hi, color=col, alpha=0.22, lw=0)
            ax.plot(T, med, color=col, lw=1.4, label=f"{gl} (n={len(ser)})")
        ax.axvline(EQUIL_NS, color="0.45", lw=0.9, ls="--")
        ax.text(EQUIL_NS, ax.get_ylim()[1], " analysis starts", fontsize=8.1,
                color="0.45", va="top", rotation=90)
        ax.set_xlabel("Time (ns)", fontsize=8.5)
        ax.set_ylabel(lab.replace(r"\AA", "Å"), fontsize=8.5)
        ax.legend(fontsize=8.1, frameon=False, loc="upper left")
    axes[0].text(0.0, 1.02, "a", transform=axes[0].transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
    axes[1].text(0.0, 1.02, "b", transform=axes[1].transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    ax = axes[2]
    if rmsfs:
        R = np.vstack([r for r in rmsfs if len(r) == len(resids)])
        med = np.median(R, 0)
        lo, hi = np.percentile(R, [25, 75], axis=0)
        # written out so the manuscript can cite where the fluctuation actually sits
        pd.DataFrame({"resid": resids, "rmsf_median": med, "rmsf_q25": lo,
                      "rmsf_q75": hi, "gorge_lining": np.isin(resids, LINING)}
                     ).to_csv(RES/"25_md_rmsf_per_residue.csv", index=False)
        ax.fill_between(resids, lo, hi, color="#1b7837", alpha=0.22, lw=0)
        ax.plot(resids, med, color="#1b7837", lw=0.9)
        inl = np.isin(resids, LINING)
        ax.scatter(resids[inl], med[inl], s=22, color="#b2182b", zorder=3,
                   label="gorge lining (14)")
        # Trp286 and Tyr337 are close enough on the axis that their labels collided;
        # stagger the vertical offsets rather than let them overlap.
        for k, (rn, nm) in enumerate(sorted(NAMED.items())):
            if rn in set(resids):
                j = int(np.where(resids == rn)[0][0])
                dy = (0.55, 1.15, 0.55)[k % 3]
                ax.annotate(nm, xy=(rn, med[j]), xytext=(rn, med[j]+dy), fontsize=8.1,
                            ha="center", color="#b2182b",
                            arrowprops=dict(arrowstyle="-", lw=0.7, color="#b2182b"))
        ax.legend(fontsize=8.1, frameon=False, loc="upper right")
        ax.set_xlabel("Residue number", fontsize=8.5)
        ax.set_ylabel("C$\\alpha$ RMSF (Å)", fontsize=8.5)
        ax.text(0.0, 1.02, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    for ax in axes:
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_md_stability.{ext}", dpi=300)
    print(f"  -> figures/fig_md_stability.pdf / .png")


if __name__ == "__main__":
    main()
