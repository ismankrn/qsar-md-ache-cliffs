"""Stage 16 — rung-4 dynamic descriptors from the production trajectories.

The point of rung 4 is to measure what a static pose cannot. Averaging rung-2 quantities over
frames would just be rung 2 with less noise, and stage 15 already showed that eight-fold
averaging buys nothing significant. So every descriptor here is either a *persistence*, a
*fluctuation*, or an *event count* -- quantities that require a trajectory to exist at all.

Grouped by what they test:

  the CAS-PAS bridging hypothesis, time-resolved
    occ_cas, occ_pas       fraction of frames contacting Trp86 / Trp286
    occ_bridge             fraction of frames contacting BOTH at once
    bridge_runs            number of separate bridging episodes per ns
    bridge_longest_ns      longest continuous bridging episode
      Rung 2 asked "does the docked pose bridge?" and found 92% did, which was why the
      criterion was replaced with axis projection. This asks how much of the time it does,
      which a single pose cannot answer.

  does the ligand stay where it was docked
    rmsd_mean, rmsd_max    ligand RMSD from the equilibrated pose, no realignment
    escaped                did it ever leave the gorge (COM beyond ESCAPE_A of start)
    depth_mean, depth_sd   position along the CAS->PAS axis, and how much it wanders

  flexibility, the thing a pose has none of
    rmsf_lig               mean per-atom fluctuation of the ligand
    rmsf_gorge             mean fluctuation of the 14 gorge-lining sidechains
    gate_chi1_sd           Tyr337 chi1 spread -- the swinging gate
    gate_flips             how often Tyr337 chi1 crosses between rotamer wells
    lig_dihedral_sd        mean spread of the ligand's own rotatable dihedrals

  contact durability, per residue
    occ_<resid>            fraction of frames within CONTACT of each lining residue

  water, which is why whole-system frames were kept
    n_bridge_water         mean count of waters H-bonding both ligand and protein at once

n_heavy is carried so that stage 10's size adjustment can be applied -- four rung-2/3
"findings" turned out to be size artefacts, and these descriptors are just as exposed.

The first EQUIL_NS of production is discarded on top of the 200 ps of restrained
equilibration. Frames are read at 10 ps for the cheap geometric terms and every
BURIAL_STRIDE-th frame for the ones needing all protein heavy atoms, which would otherwise
dominate the cost.

Usage:  python scripts/rung4_descriptors.py [--index I --total N] [--only RUN]
"""
from __future__ import annotations
import argparse, json, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np, pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array, capped_distance

ROOT = Path(__file__).resolve().parent.parent
RES, MD = ROOT/"results", ROOT/"md"
RUNS = MD/"runs"
RECEPTOR = ROOT/"docking"/"4ey7_rec_md.pdb"
OUT = RES/"16_rung4_per_run"

CONTACT = 4.5          # A, same cutoff as rung 2 so the two are comparable
EQUIL_NS = 5.0         # discard the first 5 ns of production
ESCAPE_A = 12.0        # A from the starting COM counts as leaving the gorge
BURIAL_STRIDE = 10     # whole-protein terms every 100 ps instead of every 10 ps
WATER_HB = 3.5         # A, heavy-atom donor/acceptor distance for a water bridge
SITES = {"cas": 86, "pas": 286, "gate": 337}
LINING = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]


def runs_of(csv):
    return pd.read_csv(csv).run.tolist()


def dihedral(p0, p1, p2, p3):
    """Dihedral in degrees for stacks of points, shape (n, 3)."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def circ_sd(deg):
    """Circular standard deviation, so a dihedral near +-180 is not reported as huge."""
    r = np.radians(deg)
    R = np.hypot(np.cos(r).mean(), np.sin(r).mean())
    return float(np.degrees(np.sqrt(max(-2*np.log(max(R, 1e-12)), 0.0))))


def analyse(run_dir: Path, traj: str = "prod", equil_ns: float = EQUIL_NS,
            min_frames: int = 100) -> dict:
    """traj/equil_ns are overridable so the extraction can be smoke-tested on the short
    equilibration trajectories before any production run exists. An extraction script used
    for the first time on the data it must not corrupt is a bad place to find a bug."""
    xtc = run_dir/f"{traj}.xtc"
    # Topology from the .tpr when it exists: a .gro carries no bonds, so lig.bonds comes back
    # empty and the ligand-dihedral descriptor silently becomes NaN. The .tpr has them.
    tpr = run_dir/f"{traj}.tpr"
    gro = run_dir/"npt.gro" if traj == "prod" else run_dir/"em.gro"
    topology = tpr if tpr.exists() else gro
    if not xtc.exists():
        raise FileNotFoundError("prod.xtc tidak ada")
    u = mda.Universe(str(topology), str(xtc))
    dt_ps = float(u.trajectory.dt)
    if dt_ps <= 0:
        raise RuntimeError("dt trajektori tidak terbaca")
    skip = int(round(equil_ns * 1000 / dt_ps))
    n_all = len(u.trajectory)
    if n_all - skip < min_frames:
        raise RuntimeError(f"terlalu sedikit frame terpakai: {n_all - skip} "
                           f"(total {n_all}, dibuang {skip})")

    # Residue selections MUST come from the receptor PDB, by atom index, not from resid
    # expressions on the trajectory. ParmEd numbers each moleculetype from 1, and the capped
    # receptor is three of them, so the .gro renumbers 535 residues into the range 1-257:
    # `resid 337` matches nothing, and a resid that does match is the wrong residue. Tyr337
    # appears in the trajectory as TYR74. Receptor atom order is preserved (verified: all
    # 8,173 residue names match, and the only 21 atom-name differences are the cap hydrogens
    # ParmEd renames HH31->H1 and CH3->C), so indices taken from the PDB apply directly.
    rec = mda.Universe(str(RECEPTOR))
    nr = rec.atoms.n_atoms
    if u.atoms.n_atoms < nr:
        raise RuntimeError(f"sistem lebih kecil dari reseptor: {u.atoms.n_atoms} < {nr}")
    if not (u.atoms.resnames[:nr] == rec.atoms.resnames).all():
        raise RuntimeError("nama residu reseptor tidak cocok dengan awal sistem -- "
                           "urutan atom berubah, indeks tidak bisa dipakai")

    def by_index(sel):
        ix = rec.select_atoms(sel).indices
        if len(ix) == 0:
            raise RuntimeError(f"seleksi reseptor kosong: {sel}")
        return u.atoms[ix]

    lig = u.select_atoms("resname UNK and not name H*")
    if lig.n_atoms == 0:
        raise RuntimeError("ligan (resname UNK) tidak ditemukan")
    prot_h = by_index("protein and not name H*")
    site = {k: by_index(f"resid {r} and not backbone and not name H*")
            for k, r in SITES.items()}
    lin = {r: by_index(f"resid {r} and not name H*") for r in LINING}
    gorge_sc = by_index("not backbone and not name H* and resid "
                        + " ".join(str(r) for r in LINING))
    gate = by_index(f"resid {SITES['gate']} and name N CA CB CG")
    if gate.n_atoms != 4:
        raise RuntimeError(f"atom chi1 Tyr337 ada {gate.n_atoms}, harus 4")
    water_o = u.select_atoms("resname HOH SOL WAT and name O OW")
    lig_pol = u.select_atoms("resname UNK and name N* O*")
    prot_pol = by_index("protein and name N* O*")

    # gorge axis, from the equilibrated frame
    u.trajectory[skip]
    cas_c, pas_c = site["cas"].center_of_geometry(), site["pas"].center_of_geometry()
    axis = pas_c - cas_c; axis /= np.linalg.norm(axis)
    ref_lig = lig.positions.copy()
    ref_com = lig.center_of_geometry()

    # ligand rotatable dihedrals, taken from bonded topology where available
    quads = []
    try:
        for b in lig.bonds:
            a1, a2 = b.atoms
            n1 = [x for x in a1.bonded_atoms if x.index != a2.index and x in lig]
            n2 = [x for x in a2.bonded_atoms if x.index != a1.index and x in lig]
            if n1 and n2:
                quads.append((n1[0].index, a1.index, a2.index, n2[0].index))
    except Exception:
        quads = []
    quads = quads[:40]

    occ = {k: 0 for k in ("cas", "pas", "bridge")}
    occ_res = {r: 0 for r in LINING}
    rmsd, depth, chi1, dih, nburial, nwat = [], [], [], [], [], []
    gcoords = []
    bridge_series = []
    coords = []
    nf = 0
    for ts in u.trajectory[skip:]:
        lx = lig.positions
        d_cas = distance_array(lx, site["cas"].positions).min()
        d_pas = distance_array(lx, site["pas"].positions).min()
        c, p = d_cas <= CONTACT, d_pas <= CONTACT
        occ["cas"] += c; occ["pas"] += p; occ["bridge"] += (c and p)
        bridge_series.append(bool(c and p))
        for r, sel in lin.items():
            occ_res[r] += distance_array(lx, sel.positions).min() <= CONTACT
        rmsd.append(float(np.sqrt(((lx - ref_lig)**2).sum(1).mean())))
        depth.append(float(np.dot(lig.center_of_geometry() - cas_c, axis)))
        g = gate.positions
        chi1.append(float(dihedral(g[0][None], g[1][None], g[2][None], g[3][None])[0]))
        if quads:
            q = np.array(quads)
            pos = u.atoms.positions
            dih.append(dihedral(pos[q[:, 0]], pos[q[:, 1]], pos[q[:, 2]], pos[q[:, 3]]))
        coords.append(lx.copy())
        gcoords.append(gorge_sc.positions.copy())
        if nf % BURIAL_STRIDE == 0:
            dm = distance_array(lx, prot_h.positions)
            nburial.append(float((dm.min(axis=1) <= CONTACT).mean()))
            if water_o.n_atoms and lig_pol.n_atoms:
                near_l = set(capped_distance(lig_pol.positions, water_o.positions,
                                             WATER_HB, return_distances=False)[:, 1])
                near_p = set(capped_distance(prot_pol.positions, water_o.positions,
                                             WATER_HB, return_distances=False)[:, 1])
                nwat.append(len(near_l & near_p))
        nf += 1

    C = np.array(coords)
    rmsf_lig = float(np.sqrt(((C - C.mean(0))**2).sum(-1).mean(0)).mean())
    G = np.array(gcoords)
    rmsf_gorge = float(np.sqrt(((G - G.mean(0))**2).sum(-1).mean(0)).mean())
    bs = np.array(bridge_series)
    runs_ = int(((bs[1:].astype(int) - bs[:-1].astype(int)) == 1).sum() + (1 if bs[0] else 0))
    longest = 0; cur = 0
    for b in bs:
        cur = cur + 1 if b else 0
        longest = max(longest, cur)
    ns_used = nf * dt_ps / 1000.0
    chi1 = np.array(chi1)
    # Wrap the DIFFERENCE, not the angle. Wrapping chi1 itself leaves 179 -> -179 looking
    # like a 358 degree jump, which reported 18 flips/ns for a dihedral whose spread was
    # 11.5 degrees -- an arithmetic impossibility that the smoke test exposed.
    dchi = np.diff(chi1)
    dchi = (dchi + 180.0) % 360.0 - 180.0
    flips = int((np.abs(dchi) > 60).sum())

    d = {
        "run": run_dir.name, "n_frames": nf, "ns_analysed": ns_used, "dt_ps": dt_ps,
        "occ_cas": occ["cas"]/nf, "occ_pas": occ["pas"]/nf, "occ_bridge": occ["bridge"]/nf,
        "bridge_runs_per_ns": runs_/ns_used if ns_used else np.nan,
        "bridge_longest_ns": longest*dt_ps/1000.0,
        "rmsd_mean": float(np.mean(rmsd)), "rmsd_max": float(np.max(rmsd)),
        "escaped": int(np.max([np.linalg.norm(c.mean(0) - ref_com) for c in C]) > ESCAPE_A),
        "depth_mean": float(np.mean(depth)), "depth_sd": float(np.std(depth)),
        "rmsf_lig": rmsf_lig, "rmsf_gorge": rmsf_gorge,
        "gate_chi1_sd": circ_sd(chi1), "gate_flips_per_ns": flips/ns_used if ns_used else np.nan,
        "lig_dihedral_sd": float(np.mean([circ_sd(np.array(dih)[:, j])
                                          for j in range(np.array(dih).shape[1])]))
                           if dih else np.nan,
        "frac_lig_buried_mean": float(np.mean(nburial)) if nburial else np.nan,
        "n_bridge_water": float(np.mean(nwat)) if nwat else np.nan,
        "n_heavy": lig.n_atoms,
    }
    for r in LINING:
        d[f"occ_{r}"] = occ_res[r]/nf
    return d


ap = argparse.ArgumentParser()
ap.add_argument("--index", type=int, default=0)
ap.add_argument("--total", type=int, default=1)
ap.add_argument("--only", type=str, default=None)
ap.add_argument("--traj", type=str, default="prod", help="prod | npt | nvt (untuk uji)")
ap.add_argument("--equil-ns", type=float, default=EQUIL_NS)
ap.add_argument("--min-frames", type=int, default=100)
ap.add_argument("--smoke", action="store_true",
                help="uji pada trajektori pendek: tulis ke JSON bersufiks, jangan timpa hasil")
a = ap.parse_args()

allruns = runs_of(RES/"14_md_runs.csv")
sel = [a.only] if a.only else allruns[a.index::a.total]
OUT.mkdir(parents=True, exist_ok=True)
print(f"tugas {a.index}/{a.total}: {len(sel)} run", flush=True)

done, skipped, failed = 0, 0, []
t0 = time.time()
for i, name in enumerate(sel, 1):
    f = OUT/(f"{name}.smoke.json" if a.smoke else f"{name}.json")
    if f.exists():
        skipped += 1; continue
    try:
        d = analyse(RUNS/name, traj=a.traj, equil_ns=a.equil_ns,
                    min_frames=a.min_frames)
        f.write_text(json.dumps(d, indent=2))
        done += 1
        if done <= 3 or done % 20 == 0:
            print(f"  {name}: {d['n_frames']} frame, occ_bridge {d['occ_bridge']:.3f}, "
                  f"rmsd {d['rmsd_mean']:.2f} A, {time.time()-t0:.0f} s", flush=True)
    except Exception as e:
        failed.append((name, f"{type(e).__name__}: {e}"))
        print(f"  {name}: GAGAL {type(e).__name__}: {e}", flush=True)

print(f"\nselesai {done}, dilewati {skipped}, gagal {len(failed)}", flush=True)
sys.exit(1 if failed else 0)
