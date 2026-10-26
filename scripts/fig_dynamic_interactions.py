#!/usr/bin/env python
"""Typed interaction occupancy over the trajectories, and what the rung-4 descriptors miss.

Figure 5 shows the reference ligand's interactions in a single pose. This asks the same
question of the dynamics, for every run, with ProLIF's interaction typing: how often is each
contact present, and of what kind. The rung-4 descriptors are distance-only occupancies --
`occ_cas`, `occ_pas`, `occ_86` and the rest count frames within 4.5 A with no distinction
between a hydrogen bond, a pi-stack and a passing van der Waals brush. This says what that
simplification keeps and what it discards.

Panel (b) is the point: ProLIF's occupancy at Trp86 and Trp286 -- that of the most
persistent interaction type there -- against the project's own occ_cas and occ_pas, one
point per run. If they agree, the distance-only descriptor is a
faithful summary and the typing adds interpretation rather than information. If they diverge,
the descriptor set has a gap worth naming.

Strictly descriptive, like Figure 5: nothing here enters any model. A typed descriptor family
introduced after the ladder was scored would be post hoc.

Three things the reduced archives force, all of which produce silently wrong answers if
skipped:

* **Periodic images.** `trjconv -pbc mol` places each molecule independently and this
  receptor is three chains, so the ligand's centre of mass wanders 4.6 to 88.9 A from the
  receptor's across a single run. ProLIF would report no interactions for most frames. The
  stage 16b `make_whole` is applied as an on-the-fly transformation.
* **Elements.** A `.gro` carries no element column, and the RDKit converter needs one.
* **Residue identity.** The archives renumber residues per chain segment, so the three chains
  carry overlapping resids; ProLIF groups by identifier and merges atoms from different
  chains into one residue, which it then reports as fragmented. Renumbering every residue
  uniquely fixes it, and a map back to 4EY7 numbering -- verified against the docking
  receptor by residue name -- keeps the labels meaningful.

Usage:
    python scripts/fig_dynamic_interactions.py --stride 10
"""
from __future__ import annotations
import argparse, re, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types
import prolif
from rung4_descriptors import receptor_molecules, make_whole

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
RED = ROOT/"md"/"reduced"
EQUIL_NS = 5.0
NAMED = {86: "Trp86 (CAS)", 286: "Trp286 (PAS)", 337: "Tyr337 (gate)"}
LINING = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]


class Whole:
    """Stateless per-frame minimum-image fix, applied as an MDAnalysis transformation."""

    def __init__(self, u, mols):
        self.u, self.mols = u, mols

    def __call__(self, ts):
        make_whole(self.u, self.mols, ts.dimensions)
        return ts


def prepare(xtc: Path, dock_ca):
    u = mda.Universe(str(xtc.with_suffix(".gro")), str(xtc))
    u.add_TopologyAttr("elements", guess_types(u.atoms.names))
    rec = u.select_atoms("protein")
    segs = receptor_molecules(rec)
    u.residues.resids = np.arange(1, len(u.residues) + 1)
    rec = u.select_atoms("protein")
    lig = u.select_atoms("resname UNK")
    if len(lig) == 0:
        return None
    ca = rec.select_atoms("name CA")
    if len(ca) != len(dock_ca) or not (ca.resnames == dock_ca.resnames).all():
        return None
    true_of = {int(a): int(b) for a, b in zip(ca.resids, dock_ca.resids)}
    u.atoms.guess_bonds()
    u.trajectory.add_transformations(Whole(u, segs + [lig.indices]))
    return u, rec, lig, true_of


def resid_of(label: str) -> int:
    m = re.match(r"[A-Z]{3}(\d+)", str(label))
    return int(m.group(1)) if m else -1


def figure(d, a) -> None:
        # ---------------- figure ----------------
        FIGS.mkdir(exist_ok=True)
        fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2),
                                 gridspec_kw={"width_ratios": [1.5, 1.0]})

        # (a) median occupancy per residue per type, gorge residues only
        g = d[d.resid_4ey7.isin(LINING)]
        # A median over one run looked exactly like a median over two hundred: the
        # Glu334 salt bridge and the Trp286 hydrogen bond were each a single
        # trajectory, printed as 87 and 52 beside cells backed by the whole campaign.
        # Cells supported by fewer than MIN_SHARE of runs are therefore blanked.
        MIN_SHARE = 0.10
        n_runs = d.run.nunique()
        piv = g.pivot_table(index="resid_4ey7", columns="kind", values="occupancy",
                            aggfunc="median").fillna(0.0)
        sup = (g.pivot_table(index="resid_4ey7", columns="kind", values="run",
                             aggfunc="nunique").reindex_like(piv).fillna(0) / n_runs)
        piv = piv.where(sup >= MIN_SHARE, 0.0)
        ax = axes[0]
        im = ax.imshow(piv.to_numpy()*100, aspect="auto", cmap="YlGnBu", vmin=0, vmax=100)
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels(piv.columns, rotation=30, ha="right", fontsize=7.8)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([f"{r}" + (f"  {NAMED[r].split()[0]}" if r in NAMED else "")
                            for r in piv.index], fontsize=7.8)
        for yi, r in enumerate(piv.index):
            for xi, c in enumerate(piv.columns):
                v = piv.iloc[yi, xi]*100
                if v >= 5:
                    ax.text(xi, yi, f"{v:.0f}", ha="center", va="center", fontsize=7.8,
                            color="white" if v > 55 else "0.25")
        plt.colorbar(im, ax=ax, label="median occupancy (%)", fraction=0.035, pad=0.02)
        ax.text(0.0, 1.02, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
        ax.set_ylabel("Gorge-lining residue (4EY7 numbering)", fontsize=8.5)

        # (b) ProLIF total occupancy vs the project's distance-only occupancy
        ax = axes[1]
        import json, glob as _g
        own = {}
        for f in _g.glob(str(RES/"16_rung4_per_run"/"*.json")):
            j = json.loads(Path(f).read_text())
            w = (j.get("windows") or {}).get("20")
            if w:
                own[j["run"]] = (w.get("occ_cas"), w.get("occ_pas"))
        # The max over interaction types, not their union: a frame may carry two types
        # at one residue, so this is a lower bound on "an interaction of some kind was
        # typed". The axis is labelled for what is plotted rather than for the union,
        # which the per-frame table would be needed to compute.
        tot = (d[d.resid_4ey7.isin([86, 286])]
               .groupby(["run", "resid_4ey7"]).occupancy.max().unstack())
        pts = []
        for run, r in tot.iterrows():
            if run in own:
                for rid, col, lab in ((86, 0, "Trp86 / occ_cas"), (286, 1, "Trp286 / occ_pas")):
                    if rid in r and np.isfinite(r[rid]) and own[run][col] is not None:
                        pts.append(dict(run=run, site=lab, prolif=float(r[rid]),
                                        own=float(own[run][col])))
        P = pd.DataFrame(pts)
        for lab, col in (("Trp86 / occ_cas", "#b2182b"), ("Trp286 / occ_pas", "#2166ac")):
            s = P[P.site == lab]
            if len(s):
                rho = s.prolif.corr(s.own, method="spearman")
                ax.scatter(s.own*100, s.prolif*100, s=14, alpha=0.6, color=col,
                           edgecolors="none", label=f"{lab}  (rho={rho:.2f}, n={len(s)})")
        ax.plot([0, 100], [0, 100], color="0.6", lw=0.9, ls="--")
        ax.set_xlabel("Distance-only occupancy, this study (%)", fontsize=8.5)
        ax.set_ylabel("ProLIF typed occupancy, most persistent type (%)", fontsize=8.5)
        ax.text(0.0, 1.02, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
        ax.legend(fontsize=7.8, frameon=False, loc="lower right")
        P.to_csv(RES/"26b_occupancy_agreement.csv", index=False)

        # The replicate ratio for the typed occupancy, computed exactly as
        # replicate_pilot.py computes it for the distance-only descriptors: a pooled
        # within-ligand SD over the six ligands that carry replicates, against the SD of
        # the per-ligand means over every ligand with a result. Section 3.7 quotes these
        # two numbers against 103 and 52, and an ad hoc calculation in a terminal is not
        # a source a reader can check, so they are written out with the figure.
        tot2 = tot.copy()
        tot2["ligand"] = [r.rsplit("_", 1)[0] for r in tot2.index]
        reps = [l for l, grp in tot2.groupby("ligand") if len(grp) >= 3]
        rr = []
        for rid, site in ((86, "Trp86"), (286, "Trp286")):
            v = pd.to_numeric(tot2[rid], errors="coerce")
            sub = tot2[tot2.ligand.isin(reps)]
            var_within = pd.to_numeric(sub[rid], errors="coerce").groupby(
                sub.ligand).var(ddof=1).dropna()
            sd_rep = float(np.sqrt(var_within.mean()))
            lig_means = v.groupby(tot2.ligand).mean().dropna()
            sd_lig = float(lig_means.std(ddof=1))
            rr.append(dict(site=site, n_ligands=int(lig_means.size),
                           n_replicate_ligands=len(reps), sd_rep=sd_rep, sd_lig=sd_lig,
                           ratio=sd_lig/sd_rep if sd_rep > 0 else float("inf")))
        R = pd.DataFrame(rr)
        R.to_csv(RES/"26c_typed_replicate_ratio.csv", index=False)
        print("\n  replicate ratio of the typed occupancy (threshold 2):")
        for r in R.itertuples():
            print(f"    {r.site:7s} sd_rep = {r.sd_rep:.4f}  sd_lig = {r.sd_lig:.4f}  "
                  f"ratio = {r.ratio:.2f}  ({r.n_ligands} ligands, "
                  f"{r.n_replicate_ligands} with replicates)")

        for ax in axes:
            ax.tick_params(labelsize=8)
        axes[1].spines[["top", "right"]].set_visible(False)

        fig.tight_layout(rect=(0, 0.04, 1, 0.94))
        for ext in ("pdf", "png"):
            fig.savefig(FIGS/f"fig_dynamic_interactions.{ext}", dpi=300)
        print(f"  -> figures/fig_dynamic_interactions.pdf / .png")
        if len(P):
            print("\n  agreement (Spearman):")
            for lab in P.site.unique():
                s = P[P.site == lab]
                print(f"    {lab:20s} rho = {s.prolif.corr(s.own, method='spearman'):+.3f}  "
                      f"median typed {s.prolif.median()*100:.0f}%  "
                      f"distance-only {s.own.median()*100:.0f}%")




def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=10)
    ap.add_argument("--limit", type=int, default=0, help="debug: only this many runs")
    ap.add_argument("--from-csv", action="store_true",
                    help="redraw from results/26_dynamic_interactions.csv without rerunning "
                         "ProLIF. The full pass over 229 trajectories takes about four "
                         "hours, so a change to the figure must not require repeating it.")
    a = ap.parse_args()

    if a.from_csv:
        d = pd.read_csv(RES/"26_dynamic_interactions.csv")
        print(f"  redrawing from {len(d)} rows, {d.run.nunique()} runs")
        figure(d, a)
        return

    dock_ca = mda.Universe(str(ROOT/"docking"/"4ey7_rec_clean.pdb")).select_atoms(
        "protein and name CA")
    xtcs = sorted(RED.glob("*.xtc"))
    if a.limit:
        xtcs = xtcs[:a.limit]
    print(f"  {len(xtcs)} runs, stride {a.stride}")

    rows = []
    for i, x in enumerate(xtcs, 1):
        try:
            prep = prepare(x, dock_ca)
            if prep is None:
                print(f"    {x.stem}: skipped (topology mismatch)")
                continue
            u, rec, lig, true_of = prep
            start = int(EQUIL_NS * 1000 / u.trajectory.dt)
            fp = prolif.Fingerprint()
            fp.run(u.trajectory[start::a.stride], lig, rec, progress=False, n_jobs=1)
            df = fp.to_dataframe()
        except Exception as e:
            print(f"    {x.stem}: {type(e).__name__}: {e}")
            continue
        if df.empty:
            continue
        occ = df.mean()
        for (_l, res, kind), v in occ.items():
            tr = true_of.get(resid_of(res), -1)
            rows.append(dict(run=x.stem, resid_4ey7=tr, resname=str(res)[:3],
                             kind=kind, occupancy=float(v), n_frames=len(df)))
        if i % 20 == 0:
            print(f"    {i}/{len(xtcs)}")

    d = pd.DataFrame(rows)
    d.to_csv(RES/"26_dynamic_interactions.csv", index=False)
    print(f"\n  -> results/26_dynamic_interactions.csv  ({len(d)} rows, "
          f"{d.run.nunique()} runs)")

    figure(d, a)


if __name__ == "__main__":
    main()
