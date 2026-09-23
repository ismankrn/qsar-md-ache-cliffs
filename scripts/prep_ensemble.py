#!/usr/bin/env python
"""Stage 11b — prepare the ensemble receptors, all in the 4EY7 frame.

Each member is superposed onto 4EY7 before preparation, so every pose that
comes out of ensemble docking lands in one common frame. That matters: the
rung 2 descriptors are defined against 4EY7's Trp86 and Trp286 positions and
the gorge axis drawn between them, so poses from different receptors have to
be directly comparable without further transformation. The same docking box
is then used throughout.
"""
from __future__ import annotations
import subprocess, sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
from pdbfixer import PDBFixer
from openmm.app import PDBFile

ROOT = Path(__file__).resolve().parent.parent
STRUCT, DOCK = ROOT/"structures", ROOT/"docking"
REF, PH, BOX = "4ey7", 7.4, 30.0
ENS = (DOCK/"ensemble.txt").read_text().split()
MK = "/opt/miniforge3/envs/qsar-md/bin/mk_prepare_receptor.py"
if not Path(MK).exists():
    MK = str(Path(sys.executable).parent/"mk_prepare_receptor.py")

ref = mda.Universe(str(STRUCT/f"{REF}.pdb"))
def ca_one(u, ch="A"):
    sel = u.select_atoms(f"protein and chainID {ch} and name CA")
    seen, idx = set(), []
    for i, a in enumerate(sel):
        if a.resid not in seen:
            seen.add(a.resid); idx.append(i)
    return sel[idx]
ref_ca = ca_one(ref)
center = np.load(DOCK/f"{REF}_box_center.npy")
print(f"ensemble: {' '.join(p.upper() for p in ENS)}")
print(f"shared box: {np.round(center,2)}, {BOX} A\n")

ok = []
for pdb in ENS:
    out_pdbqt = DOCK/f"ens_{pdb}_rec.pdbqt"
    if out_pdbqt.exists():
        print(f"  {pdb.upper()}: already prepared"); ok.append(pdb); continue
    try:
        u = mda.Universe(str(STRUCT/f"{pdb}.pdb"))
        m_ca = ca_one(u)
        common = np.intersect1d(ref_ca.residues.resids, m_ca.residues.resids)
        rp = ref_ca.positions[np.isin(ref_ca.resids, common)]
        mp = m_ca.positions[np.isin(m_ca.resids, common)]
        n = min(len(rp), len(mp)); rp, mp = rp[:n], mp[:n]
        rc, mc = rp.mean(0), mp.mean(0)
        R, rms = mda_align.rotation_matrix(mp - mc, rp - rc)

        _all = u.select_atoms("protein and chainID A")
        first, last = int(_all.residues.resids.min()), int(_all.residues.resids.max())
        # the C-terminal residue is disordered and PDBFixer places its OXT
        # ambiguously between C and CB, which makes a 5-valent carbon
        rec = u.select_atoms(f"protein and chainID A and resid {first}:{last-1}")
        moved = rec.copy()
        moved.atoms.positions = (rec.positions - mc) @ R.T + rc
        raw = DOCK/f"ens_{pdb}_raw.pdb"
        moved.atoms.write(str(raw))

        fx = PDBFixer(filename=str(raw))
        fx.findMissingResidues()
        chains = list(fx.topology.chains())
        for k in list(fx.missingResidues):
            ch = chains[k[0]]
            if k[1] == 0 or k[1] == len(list(ch.residues())):
                del fx.missingResidues[k]
        fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
        fx.removeHeterogens(keepWater=False)
        fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(PH)
        clean = DOCK/f"ens_{pdb}_clean.pdb"
        with open(clean, "w") as fh:
            PDBFile.writeFile(fx.topology, fx.positions, fh, keepIds=True)

        r = subprocess.run(
            [sys.executable, MK, "--read_pdb", str(clean),
             "-o", str(DOCK/f"ens_{pdb}_rec"), "-p",
             "--box_size", str(BOX), str(BOX), str(BOX),
             "--box_center", *[f"{c:.3f}" for c in center]],
            capture_output=True, text=True)
        if r.returncode == 0 and out_pdbqt.exists():
            print(f"  {pdb.upper()}: superposed (CA RMSD {rms:.2f} A), prepared")
            ok.append(pdb)
        else:
            print(f"  {pdb.upper()}: meeko failed -- {r.stderr.strip()[-90:]}")
    except Exception as e:
        print(f"  {pdb.upper()}: {type(e).__name__}: {str(e)[:80]}")

(DOCK/"ensemble_ready.txt").write_text("\n".join(ok))
print(f"\nready: {len(ok)} of {len(ENS)} receptors")
print(f"written: docking/ensemble_ready.txt")
