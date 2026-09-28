"""Stage 16c — QC plots for the rung-4 production runs.

Why this exists: stage 16 shipped with an RMSD that measured the protein tumbling in the box
rather than the ligand moving in the pocket, and it reported 38 A of fluctuation for buried
sidechains for two completed runs before anyone looked. A single plotted line would have made
that impossible to miss. See 16b_rung4_pbc_fix.md.

So this is not decoration. It is the check that the extraction has to pass, and it is built to
be unable to pass it for the wrong reason: the unwrapping, the periodic-image placement and the
Ca superposition are IMPORTED from rung4_descriptors.py, not reimplemented here. A QC plot with
its own copy of that geometry could agree with a broken extraction, which is the one failure
mode a QC plot must not have.

Two outputs:

  results/16c_qc/<run>.png    per-run time series, from the reduced archive
  results/16c_qc_summary.png  fleet-wide distributions, from the per-run JSONs
  results/16c_qc_flags.csv    every run that trips a threshold, with the reason

The per-run panels come from md/reduced/ (protein + ligand, 100 ps). That is coarser than the
10 ps the descriptors use and holds no water, so the numbers recomputed here are not expected
to match the JSON to the digit -- the panel prints both side by side precisely so that a gross
disagreement is visible instead of averaged away.

Usage:
  python scripts/rung4_qc_plots.py                 # everything available
  python scripts/rung4_qc_plots.py --only RUN      # one run's panel
  python scripts/rung4_qc_plots.py --fleet-only    # skip the per-run panels
"""
from __future__ import annotations
import argparse, json, sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rung4_descriptors as r4          # the geometry is shared, deliberately

ROOT = Path(__file__).resolve().parent.parent
REDUCED = ROOT/"md"/"reduced"
JSONS = ROOT/"results"/"16_rung4_per_run"
OUTDIR = ROOT/"results"/"16c_qc"
SUMMARY = ROOT/"results"/"16c_qc_summary.png"
FLAGS = ROOT/"results"/"16c_qc_flags.csv"

# Thresholds. The first two are physical impossibilities, not preferences: AChE is a rigid
# folded protein and its buried sidechains do not fluctuate by nanometres. They exist to catch
# a regression of the 16b defect, whose signature was exactly rmsf_gorge ~ 38 A.
CA_RMSD_BROKEN = 8.0       # A -- receptor not reassembled, or the fit failed
RMSF_GORGE_BROKEN = 5.0    # A -- buried sidechains cannot do this
# These are judgements about the simulation, not about the code.
CA_RMSD_HIGH = 3.0         # A -- receptor drifting further than AChE should
LIG_RMSD_HIGH = 5.0        # A -- ligand no longer in the pose it was docked in
NS_MIN = 40.0              # ns of the expected 45 after discarding equilibration


def series_from_archive(run: str):
    """Ligand and receptor time series from one reduced archive, in the receptor's own frame.

    Mirrors analyse() step for step, using its helpers, so the two cannot drift apart.
    """
    gro, xtc = REDUCED/f"{run}.gro", REDUCED/f"{run}.xtc"
    if not (gro.exists() and xtc.exists()):
        raise FileNotFoundError(f"reduced archive missing for {run}")
    u = mda.Universe(str(gro), str(xtc))
    rec = mda.Universe(str(r4.RECEPTOR))
    nr = rec.atoms.n_atoms
    if not (u.atoms.resnames[:nr] == rec.atoms.resnames).all():
        raise RuntimeError(f"{run}: receptor residue names do not match the start of the "
                           f"archive -- atom order changed, indices are not usable")

    def by_index(sel):
        ix = rec.select_atoms(sel).indices
        if len(ix) == 0:
            raise RuntimeError(f"empty receptor selection: {sel}")
        return u.atoms[ix]

    lig = u.select_atoms("resname UNK and not name H*")
    ca = by_index("protein and name CA")
    site = {k: by_index(f"resid {v} and not backbone and not name H*")
            for k, v in r4.SITES.items()}
    # Sidechain-only, exactly as the descriptor's gorge group is defined. Gly121 has no
    # sidechain heavy atom and so contributes nothing -- selecting it on its own raises, which
    # is why the per-residue split is taken as masks over this one group instead.
    gorge_sel = "not backbone and not name H* and resid " + " ".join(str(r) for r in r4.LINING)
    gorge = by_index(gorge_sel)
    # Residue numbers come from the receptor PDB: the trajectory's own resids are renumbered
    # per moleculetype and five of them are duplicated.
    gorge_resids = rec.select_atoms(gorge_sel).resids
    mols = r4.receptor_molecules(rec) + [u.select_atoms("resname UNK").indices]
    whole = u.atoms[np.concatenate(mols)]
    try:
        has_bonds = len(whole.bonds) > 0
    except Exception:
        has_bonds = False       # a .gro carries none; trjconv -pbc mol already made them whole

    dt_ps = float(u.trajectory.dt)
    skip = int(round(r4.EQUIL_NS * 1000 / dt_ps))
    u.trajectory[0]
    r4.unwrap_mols(whole, has_bonds)
    r4.make_whole(u, mols, u.trajectory.ts.dimensions)
    prev = whole.positions.copy()

    # Pass one: unwrap forward and keep the lab-frame coordinates. The superposition cannot be
    # done inside this loop, because the reference is the frame at EQUIL_NS -- the same
    # reference analyse() uses. Taking frame 0 instead would make the plotted RMSD and the
    # JSON's rmsd_mean two different quantities drawn on one axis, which is worse than not
    # plotting the JSON at all.
    t, depth, cas_on, pas_on = [], [], [], []
    ca_xyz, lig_lab, gorge_lab = [], [], []
    pas_sep = np.nan
    for i, ts in enumerate(u.trajectory):
        box = ts.dimensions
        r4.nojump_into(whole, prev, box)
        r4.make_whole(u, mols, box)
        prev = whole.positions.copy()
        lx = lig.positions
        t.append(i*dt_ps/1000.0)
        cas_on.append(distance_array(lx, site["cas"].positions, box=box).min() <= r4.CONTACT)
        pas_on.append(distance_array(lx, site["pas"].positions, box=box).min() <= r4.CONTACT)
        cas_c = site["cas"].center_of_geometry()
        ax = site["pas"].center_of_geometry() - cas_c
        sep = float(np.linalg.norm(ax)); ax = ax/sep
        if i == skip or not np.isfinite(pas_sep):
            pas_sep = sep
        depth.append(float(np.dot(lx.mean(0) - cas_c, ax)))
        ca_xyz.append(ca.positions.copy())
        lig_lab.append(lx.copy())
        gorge_lab.append(gorge.positions.copy())

    # Pass two: superpose every frame on the reference and carry the ligand with it.
    ca_xyz, lig_lab, gorge_lab = map(np.array, (ca_xyz, lig_lab, gorge_lab))
    ref_i = min(skip, len(ca_xyz) - 1)
    ca_ref_c = ca_xyz[ref_i].mean(0)
    ca_ref = ca_xyz[ref_i] - ca_ref_c
    ref_lig = lig_lab[ref_i] - ca_ref_c
    ca_rmsd, lig_rmsd, L, G = [], [], [], []
    for P, lx, gx in zip(ca_xyz, lig_lab, gorge_lab):
        Pc = P.mean(0)
        U = r4.fit_rotation(P - Pc, ca_ref)
        ca_rmsd.append(float(np.sqrt((((P - Pc) @ U - ca_ref)**2).sum(1).mean())))
        li = (lx - Pc) @ U
        lig_rmsd.append(float(np.sqrt(((li - ref_lig)**2).sum(1).mean())))
        L.append(li); G.append((gx - Pc) @ U)
    ca_rmsd, lig_rmsd = np.array(ca_rmsd), np.array(lig_rmsd)
    L = np.array(L); G = np.array(G)
    # Per-residue sidechain fluctuation. The gorge group is one flat selection, so each
    # residue is a column mask into it rather than a separate pass over the trajectory.
    gorge_res = []
    for r in r4.LINING:
        mask = gorge_resids == r
        if not mask.any():
            continue        # glycine: no sidechain heavy atom to fluctuate
        g = G[:, mask, :]
        gorge_res.append((r, float(np.sqrt(((g - g.mean(0))**2).sum(-1).mean(0)).mean())))

    return dict(
        run=run, t=np.array(t), dt_ps=dt_ps, pas_sep=pas_sep,
        ca_rmsd=ca_rmsd, lig_rmsd=lig_rmsd, depth=np.array(depth), equil_ns=r4.EQUIL_NS,
        cas_on=np.array(cas_on), pas_on=np.array(pas_on),
        rmsf_atom=np.sqrt(((L - L.mean(0))**2).sum(-1).mean(0)),
        lig_names=list(lig.names),
        rmsf_gorge_res=gorge_res,
    )


def panel(s: dict, ref: dict | None, out: Path):
    """One page per run. Equilibration is shaded rather than cut, because a run that has not
    settled is exactly what this is meant to reveal."""
    eq = r4.EQUIL_NS
    keep = s["t"] >= eq
    fig, ax = plt.subplots(3, 2, figsize=(13, 10))
    fig.suptitle(f"{s['run']} — rung-4 QC   (reduced archive, {s['dt_ps']:.0f} ps/frame)",
                 fontsize=13, fontweight="bold")

    def shade(a):
        a.axvspan(0, eq, color="0.9", zorder=0)
        a.set_xlabel("time (ns)")

    a = ax[0, 0]
    a.plot(s["t"], s["ca_rmsd"], lw=0.8, color="#333")
    a.axhline(CA_RMSD_HIGH, ls=":", color="orange", label=f"{CA_RMSD_HIGH:.0f} Å drifting")
    a.axhline(CA_RMSD_BROKEN, ls="--", color="crimson", label=f"{CA_RMSD_BROKEN:.0f} Å broken")
    a.set_ylabel("receptor Cα RMSD (Å)"); a.set_title("receptor stability")
    a.set_ylim(0, max(CA_RMSD_BROKEN*1.1, s["ca_rmsd"].max()*1.1)); a.legend(fontsize=7)
    shade(a)

    a = ax[0, 1]
    a.plot(s["t"], s["lig_rmsd"], lw=0.8, color="#1f77b4")
    a.axhline(s["lig_rmsd"][keep].mean(), ls="-", color="k", lw=0.8,
              label=f"archive mean {s['lig_rmsd'][keep].mean():.2f} Å")
    top = max(1.0, s["lig_rmsd"].max()*1.15)
    if ref and np.isfinite(ref.get("rmsd_mean", np.nan)):
        jm = ref["rmsd_mean"]
        if jm <= top:
            a.axhline(jm, ls="--", color="crimson", lw=0.8, label=f"JSON mean {jm:.2f} Å")
        else:
            # A JSON value off the top of the data's own scale is itself the finding. Letting
            # it set the axis would flatten the trace that shows why it is wrong.
            a.plot([], [], ls="--", color="crimson", label=f"JSON mean {jm:.1f} Å (off scale)")
    a.set_ylim(0, top)
    a.set_ylabel("ligand RMSD (Å)"); a.set_title("ligand, in the receptor's frame")
    a.legend(fontsize=7); shade(a)

    a = ax[1, 0]
    a.plot(s["t"], s["depth"], lw=0.8, color="#2ca02c")
    a.axhline(0, ls="--", color="k", lw=0.8, label="CAS (Trp86)")
    a.axhline(s["pas_sep"], ls="--", color="purple", lw=0.8, label="PAS (Trp286)")
    a.set_ylabel("depth along CAS→PAS axis (Å)"); a.set_title("position in the gorge")
    a.legend(fontsize=7); shade(a)

    a = ax[1, 1]
    n = len(s["rmsf_atom"])
    a.bar(np.arange(n), s["rmsf_atom"], color="#1f77b4", width=0.8)
    a.set_xticks(np.arange(n)); a.set_xticklabels(s["lig_names"], rotation=90, fontsize=5)
    a.set_ylabel("RMSF (Å)"); a.set_title("per-atom ligand fluctuation")
    a.axhline(s["rmsf_atom"].mean(), ls="--", color="k", lw=0.8,
              label=f"mean {s['rmsf_atom'].mean():.2f} Å")
    a.legend(fontsize=7)

    a = ax[2, 0]
    res = sorted(s["rmsf_gorge_res"], key=lambda x: x[0])
    a.bar(range(len(res)), [v for _, v in res], color="#9467bd")
    a.set_xticks(range(len(res))); a.set_xticklabels([str(r) for r, _ in res], fontsize=7)
    a.axhline(RMSF_GORGE_BROKEN, ls="--", color="crimson", lw=0.8,
              label=f"{RMSF_GORGE_BROKEN:.0f} Å non-physical")
    a.set_xlabel("gorge residue"); a.set_ylabel("sidechain RMSF (Å)")
    a.set_title("gorge sidechain fluctuation"); a.legend(fontsize=7)

    a = ax[2, 1]
    for row, (lab, v, col) in enumerate((("bridge", s["cas_on"] & s["pas_on"], "#d62728"),
                                         ("PAS", s["pas_on"], "purple"),
                                         ("CAS", s["cas_on"], "#2ca02c"))):
        a.fill_between(s["t"], row, row + 0.8, where=v, color=col, step="mid", lw=0)
        a.text(-0.01, row + 0.4, lab, ha="right", va="center", fontsize=8,
               transform=a.get_yaxis_transform())
    a.set_yticks([]); a.set_ylim(-0.1, 3.0)
    occ = (f"CAS {s['cas_on'][keep].mean():.3f}   PAS {s['pas_on'][keep].mean():.3f}   "
           f"bridge {(s['cas_on'] & s['pas_on'])[keep].mean():.3f}")
    if ref:
        occ += (f"\nJSON:  CAS {ref.get('occ_cas', float('nan')):.3f}   "
                f"PAS {ref.get('occ_pas', float('nan')):.3f}   "
                f"bridge {ref.get('occ_bridge', float('nan')):.3f}")
    a.set_title("contact occupancy\n" + occ, fontsize=8)
    shade(a)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110); plt.close(fig)


def check(d: dict) -> list[str]:
    """Threshold checks against one run's JSON. Order matters: a defect is reported as a
    defect, not as an unstable simulation."""
    f = []
    g = lambda k: d.get(k, float("nan"))
    if g("ca_rmsd_mean") > CA_RMSD_BROKEN:
        f.append(f"EXTRACTION: Ca RMSD {g('ca_rmsd_mean'):.1f} A, receptor not reassembled")
    if g("rmsf_gorge") > RMSF_GORGE_BROKEN:
        f.append(f"EXTRACTION: gorge RMSF {g('rmsf_gorge'):.1f} A is non-physical")
    if "ca_rmsd_mean" not in d:
        f.append("EXTRACTION: no ca_rmsd_mean, written by a pre-16b version of the script")
    if CA_RMSD_HIGH < g("ca_rmsd_mean") <= CA_RMSD_BROKEN:
        f.append(f"receptor drifting, Ca RMSD {g('ca_rmsd_mean'):.1f} A")
    if g("rmsd_mean") > LIG_RMSD_HIGH:
        f.append(f"ligand left its pose, RMSD {g('rmsd_mean'):.1f} A")
    if g("escaped") == 1:
        f.append("ligand left the gorge at some point")
    if g("ns_analysed") < NS_MIN:
        f.append(f"only {g('ns_analysed'):.1f} ns analysed")
    if g("occ_cas") + g("occ_pas") == 0:
        f.append("no contact with either CAS or PAS in any frame")
    return f


def fleet(records: list[dict]):
    """Distributions across runs, and the flag table. With 217 runs the per-run pages are not
    something anyone reads end to end; this is the page that says where to look."""
    import csv
    rows = []
    for d in records:
        fl = check(d)
        if fl:
            rows.append({"run": d.get("run", "?"),
                         "ca_rmsd_mean": d.get("ca_rmsd_mean"),
                         "rmsd_mean": d.get("rmsd_mean"),
                         "rmsf_gorge": d.get("rmsf_gorge"),
                         "escaped": d.get("escaped"),
                         "ns_analysed": d.get("ns_analysed"),
                         "flags": "; ".join(fl)})
    with FLAGS.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["run", "ca_rmsd_mean", "rmsd_mean", "rmsf_gorge",
                                           "escaped", "ns_analysed", "flags"])
        w.writeheader(); w.writerows(rows)

    panels = [("ca_rmsd_mean", "receptor Cα RMSD (Å)", CA_RMSD_BROKEN),
              ("rmsd_mean", "ligand RMSD (Å)", LIG_RMSD_HIGH),
              ("rmsf_lig", "ligand RMSF (Å)", None),
              ("rmsf_gorge", "gorge sidechain RMSF (Å)", RMSF_GORGE_BROKEN),
              ("occ_bridge", "CAS–PAS bridge occupancy", None),
              ("n_bridge_water", "bridging waters", None)]
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    fig.suptitle(f"rung-4 QC across {len(records)} completed runs", fontsize=13,
                 fontweight="bold")
    for a, (key, label, thr) in zip(axes.ravel(), panels):
        v = np.array([d.get(key, np.nan) for d in records], dtype=float)
        v = v[np.isfinite(v)]
        if len(v) == 0:
            a.text(0.5, 0.5, f"no {key}", ha="center", va="center"); a.set_axis_off(); continue
        a.hist(v, bins=min(30, max(5, len(v)//3)), color="#1f77b4", edgecolor="white")
        if thr is not None and v.max() > thr*0.5:
            a.axvline(thr, ls="--", color="crimson", lw=1)
        a.set_xlabel(label); a.set_ylabel("runs")
        a.set_title(f"median {np.median(v):.2f}   max {v.max():.2f}", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(SUMMARY, dpi=110); plt.close(fig)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=str, default=None, help="one run name")
    ap.add_argument("--fleet-only", action="store_true", help="skip the per-run panels")
    ap.add_argument("--new-only", action="store_true",
                    help="only render panels whose archive is newer than the existing page; "
                         "this is what makes the check cheap enough to run on every pull")
    ap.add_argument("--reduced", type=str, default=None)
    ap.add_argument("--jsons", type=str, default=None)
    a = ap.parse_args()
    if a.reduced: REDUCED = Path(a.reduced)
    if a.jsons:   JSONS = Path(a.jsons)

    records = []
    for f in sorted(JSONS.glob("*.json")):
        if f.name.endswith(".smoke.json"):
            continue
        records.append(json.loads(f.read_text()))
    print(f"JSON records: {len(records)}")

    if records:
        rows = fleet(records)
        print(f"fleet summary -> {SUMMARY.relative_to(ROOT)}")
        print(f"flagged runs : {len(rows)} / {len(records)} -> {FLAGS.relative_to(ROOT)}")
        for r in rows[:10]:
            print(f"  {r['run']}: {r['flags']}")

    if not a.fleet_only:
        runs = [a.only] if a.only else sorted(p.stem for p in REDUCED.glob("*.xtc"))
        if a.new_only:
            fresh = []
            for run in runs:
                png, xtc = OUTDIR/f"{run}.png", REDUCED/f"{run}.xtc"
                if not png.exists() or png.stat().st_mtime < xtc.stat().st_mtime:
                    fresh.append(run)
            skipped = len(runs) - len(fresh)
            runs = fresh
            if skipped:
                print(f"panels already current: {skipped}")
        by_run = {d.get("run"): d for d in records}
        ok = 0
        for run in runs:
            try:
                s = series_from_archive(run)
                panel(s, by_run.get(run), OUTDIR/f"{run}.png")
                ok += 1
                k = s["t"] >= s["equil_ns"]
                print(f"  {run}: Ca RMSD {s['ca_rmsd'][k].mean():.2f} A, "
                      f"ligand RMSD {s['lig_rmsd'][k].mean():.2f} A -> "
                      f"{(OUTDIR/f'{run}.png').relative_to(ROOT)}")
            except Exception as e:
                print(f"  {run}: FAILED {type(e).__name__}: {e}")
        print(f"per-run panels: {ok} / {len(runs)}")

    # Printed last and kept to one line, because the caller that matters is the periodic pull
    # monitor, which only shows the tail of this output. A defect has to be in that tail.
    bad = [d for d in records if any(f.startswith("EXTRACTION") for f in check(d))]
    if bad:
        print(f"QC: *** {len(bad)} / {len(records)} runs have an EXTRACTION defect: "
              + ", ".join(d.get("run", "?") for d in bad[:5])
              + (" ..." if len(bad) > 5 else ""))
        sys.exit(2)
    warn = [d for d in records if check(d)]
    print(f"QC: {len(records) - len(warn)} / {len(records)} clean"
          + (f", {len(warn)} flagged (simulation, not extraction)" if warn else ""))
