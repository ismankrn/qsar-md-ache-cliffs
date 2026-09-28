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
    rmsd_mean, rmsd_max    ligand RMSD from the equilibrated pose, measured in the
                           receptor's own frame: the receptor is reassembled across the
                           periodic boundary and superposed on the reference by Ca before
                           the ligand is compared. Measuring this in lab coordinates
                           instead reported 45 A for a ligand that contacted Trp86 in
                           99.9% of frames -- it was the protein tumbling, not the ligand.
    escaped                did it ever leave the gorge (COM beyond ESCAPE_A of start)
    depth_mean, depth_sd   position along the CAS->PAS axis, and how much it wanders.
                           The axis is recomputed each frame, which makes it invariant to
                           the protein's rotation without needing any superposition.
    ca_rmsd_mean           receptor Ca RMSD after superposition -- a QC number, not a
                           hypothesis. It is the one quantity that catches a receptor left
                           broken across the box: physical values are 1-3 A.

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


def receptor_molecules(rec):
    """Index arrays of the receptor's separate molecules, found by walking peptide bonds.

    The capped receptor is three moleculetypes, and GROMACS wraps each one into the box
    independently. Which atoms belong together therefore has to be known before any
    coordinate can be trusted. Resid ranges are not usable for this -- the caps reuse
    numbers, so five resids are duplicated -- but atom order is preserved, so index
    ranges taken from the PDB apply to the trajectory directly.
    """
    res = rec.residues
    cuts = [0]
    for i in range(len(res) - 1):
        try:
            c = res[i].atoms.select_atoms("name C").positions[0]
            n = res[i + 1].atoms.select_atoms("name N").positions[0]
            gap = float(np.linalg.norm(c - n))
        except IndexError:
            gap = 99.0
        if gap > 2.0:
            cuts.append(i + 1)
    cuts.append(len(res))
    mols = [res[a:b].atoms.indices for a, b in zip(cuts[:-1], cuts[1:])]
    covered = sum(len(m) for m in mols)
    if covered != rec.atoms.n_atoms:
        raise RuntimeError(f"pemisahan molekul reseptor tidak lengkap: {covered} "
                           f"dari {rec.atoms.n_atoms} atom")
    return mols


def make_whole(u, mols, box):
    """Put every molecule into the first molecule's periodic image, in place.

    This assumes each molecule is already internally whole. In the raw mdrun output it is
    NOT: the first receptor molecule came back spanning 48 A more than it does in the PDB,
    i.e. cut in half by the box face. Bond-based unwrapping has to run first; see unwrap_mols.

    mdrun writes each molecule wrapped into the unit cell on its own, so the ligand can sit
    a whole box vector from the pocket it never left -- 25% of frames in the first run, which
    is what made the unfitted RMSD read 45 A while the ligand contacted Trp86 in 99.9% of
    frames. Molecule centres are tens of A apart against a ~99 A box, so rounding the
    separation to whole box vectors is unambiguous.
    """
    L = box[:3]
    pos = u.atoms.positions
    ref_c = pos[mols[0]].mean(0)
    moved = False
    for ix in mols[1:]:
        shift = np.round((pos[ix].mean(0) - ref_c) / L) * L
        if shift.any():
            pos[ix] -= shift
            moved = True
    if moved:
        u.atoms.positions = pos


def unwrap_mols(grp, has_bonds):
    """Make each molecule internally whole again, following its bonds.

    A molecule straddling a box face is written with its two halves a box vector apart.
    Only the bond graph says which atoms belong to which half, so this needs a topology
    that carries bonds -- the .tpr does, a .gro does not. On a reduced archive already
    passed through `trjconv -pbc mol` the molecules are whole and this is a no-op, which
    is why the reduced-archive validation could not have caught the split.

    Called ONCE, on the reference frame. It costs 397 ms per frame against 4 ms to read
    one, so calling it per frame would add 30 min to every run -- 108 GPU-hours over the
    design, for a result that nojump_into reproduces to the last digit.
    """
    if has_bonds:
        grp.unwrap(compound="fragments", inplace=True)


def nojump_into(grp, prev, box):
    """Carry a whole molecule forward across the periodic boundary, per atom, in place.

    Once the reference frame is whole, every later frame differs from the one before it by
    sub-angstrom thermal motion plus, occasionally, a whole box vector where mdrun rewrapped
    something. Rounding the per-atom displacement to box vectors removes exactly the latter
    and nothing else. 10 ps between frames against a ~99 A box leaves four orders of
    magnitude of margin, so this cannot silently mistake motion for a jump.
    """
    pos = grp.positions
    pos -= np.round((pos - prev) / box[:3]) * box[:3]
    grp.positions = pos
    return pos


def fit_rotation(P, Q):
    """Rotation U with P @ U superposed on Q. Both inputs must already be centred."""
    V, _, Wt = np.linalg.svd(P.T @ Q)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(V @ Wt))])
    return V @ D @ Wt


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

    ca = by_index("protein and name CA")
    mols = receptor_molecules(rec) + [u.select_atoms("resname UNK").indices]

    # The atoms whose mutual geometry is measured: receptor molecules plus the ligand.
    # Waters are left wrapped and handled by minimum image instead -- unwrapping 30,000 of
    # them every frame would cost far more than the one distance call that needs them.
    whole = u.atoms[np.concatenate(mols)]
    try:
        has_bonds = len(whole.bonds) > 0
    except Exception:
        has_bonds = False

    # Reference frame for everything that compares against a fixed structure.
    u.trajectory[skip]
    unwrap_mols(whole, has_bonds)
    make_whole(u, mols, u.trajectory.ts.dimensions)
    # Verification, not decoration: this is the check that failed on the raw trajectory and
    # exposed that per-molecule image placement alone was not enough.
    for i, ix in enumerate(mols[:-1]):
        grew = (np.ptp(u.atoms[ix].positions, axis=0)
                - np.ptp(rec.atoms[ix].positions, axis=0))
        if (grew > 15.0).any():
            raise RuntimeError(f"molekul reseptor {i} masih terbentang {grew.round(1)} A lebih "
                               f"dari di PDB setelah unwrap"
                               + ("" if has_bonds else
                                  " -- topologi tanpa ikatan, unwrap berbasis ikatan tidak "
                                  "mungkin; pakai .tpr atau trajektori yang sudah -pbc mol"))
    prev_whole = whole.positions.copy()
    ca_ref = ca.positions.copy(); ca_ref_c = ca_ref.mean(0); ca_ref = ca_ref - ca_ref_c
    ref_lig = lig.positions.copy() - ca_ref_c     # ligand in the receptor's own frame
    ref_com = ref_lig.mean(0)

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
    ca_rmsd = []
    for ts in u.trajectory[skip:]:
        box = ts.dimensions
        nojump_into(whole, prev_whole, box)
        make_whole(u, mols, box)
        prev_whole = whole.positions.copy()
        lx = lig.positions
        # box= on every distance: waters and the three receptor molecules are wrapped
        # independently, so a neighbour can be an image away without minimum-image.
        d_cas = distance_array(lx, site["cas"].positions, box=box).min()
        d_pas = distance_array(lx, site["pas"].positions, box=box).min()
        c, p = d_cas <= CONTACT, d_pas <= CONTACT
        occ["cas"] += c; occ["pas"] += p; occ["bridge"] += (c and p)
        bridge_series.append(bool(c and p))
        for r, sel in lin.items():
            occ_res[r] += distance_array(lx, sel.positions, box=box).min() <= CONTACT

        # Depth needs no superposition: taking the axis from the same frame it measures
        # makes it invariant to the protein's tumbling by construction. The old code froze
        # the axis in lab coordinates at frame 0, so depth_sd measured rotation.
        cas_c = site["cas"].center_of_geometry()
        axis = site["pas"].center_of_geometry() - cas_c
        axis /= np.linalg.norm(axis)
        depth.append(float(np.dot(lx.mean(0) - cas_c, axis)))

        # RMSD and RMSF do need it -- they compare against a frozen reference, so the
        # receptor is rotated back onto the reference orientation and the ligand carried
        # along with it. Without this the protein's free rotation over 45 ns dominates.
        P = ca.positions; Pc = P.mean(0)
        U = fit_rotation(P - Pc, ca_ref)
        ca_rmsd.append(float(np.sqrt((((P - Pc) @ U - ca_ref)**2).sum(1).mean())))
        li = (lx - Pc) @ U
        rmsd.append(float(np.sqrt(((li - ref_lig)**2).sum(1).mean())))
        coords.append(li)
        gcoords.append((gorge_sc.positions - Pc) @ U)

        g = gate.positions
        chi1.append(float(dihedral(g[0][None], g[1][None], g[2][None], g[3][None])[0]))
        if quads:
            q = np.array(quads)
            pos = u.atoms.positions
            dih.append(dihedral(pos[q[:, 0]], pos[q[:, 1]], pos[q[:, 2]], pos[q[:, 3]]))
        if nf % BURIAL_STRIDE == 0:
            dm = distance_array(lx, prot_h.positions, box=box)
            nburial.append(float((dm.min(axis=1) <= CONTACT).mean()))
            if water_o.n_atoms and lig_pol.n_atoms:
                near_l = set(capped_distance(lig_pol.positions, water_o.positions,
                                             WATER_HB, box=box,
                                             return_distances=False)[:, 1])
                near_p = set(capped_distance(prot_pol.positions, water_o.positions,
                                             WATER_HB, box=box,
                                             return_distances=False)[:, 1])
                nwat.append(len(near_l & near_p))
        nf += 1

    # A split molecule, or a fit that silently failed, shows up here and nowhere else:
    # a correctly reassembled AChE fits to 1-3 A, a broken one to tens of A.
    ca_rmsd = np.array(ca_rmsd)
    if ca_rmsd.max() > 8.0:
        raise RuntimeError(f"RMSD Ca setelah superposisi mencapai {ca_rmsd.max():.1f} A -- "
                           f"reseptor tidak utuh atau superposisi gagal")

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
        "ca_rmsd_mean": float(ca_rmsd.mean()),
    }
    for r in LINING:
        d[f"occ_{r}"] = occ_res[r]/nf
    return d

# Guarded so the geometry helpers above can be imported by the QC plot script. A QC plot
# that reimplements the unwrapping and the fit would be able to agree with a broken
# extraction, which is the one thing it must not do.
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, default=0)
    ap.add_argument("--total", type=int, default=1)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--traj", type=str, default="prod", help="prod | npt | nvt (untuk uji)")
    ap.add_argument("--equil-ns", type=float, default=EQUIL_NS)
    ap.add_argument("--min-frames", type=int, default=100)
    ap.add_argument("--runs-dir", type=str, default=None,
                    help="direktori berisi direktori run (default: md/runs di repositori)")
    ap.add_argument("--receptor", type=str, default=None,
                    help="PDB reseptor untuk pemetaan indeks atom")
    ap.add_argument("--out", type=str, default=None, help="direktori keluaran JSON")
    ap.add_argument("--manifest", type=str, default=None,
                    help="CSV manifest; tidak diperlukan kalau --only dipakai")
    ap.add_argument("--smoke", action="store_true",
                    help="uji pada trajektori pendek: tulis ke JSON bersufiks, jangan timpa hasil")
    a = ap.parse_args()

    # Paths are overridable so this runs standalone on a rented GPU instance, where the
    # repository layout does not exist -- only the tpr, the receptor and a work directory.
    if a.runs_dir: RUNS = Path(a.runs_dir)
    if a.receptor: RECEPTOR = Path(a.receptor)
    if a.out:      OUT = Path(a.out)

    # The manifest is only needed to enumerate runs; with --only there is nothing to enumerate,
    # and requiring it would make the script depend on a file the instance has no reason to hold.
    if a.only:
        sel = [a.only]
    else:
        manifest = Path(a.manifest) if a.manifest else RES/"14_md_runs.csv"
        if not manifest.exists():
            sys.exit(f"FATAL: manifest tidak ada: {manifest} (pakai --only atau --manifest)")
        sel = runs_of(manifest)[a.index::a.total]
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
