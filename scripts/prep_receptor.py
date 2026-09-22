#!/usr/bin/env python
"""Siapkan reseptor AChE untuk docking.

Struktur 4EY7 (hAChE + donepezil, 2,35 A) dipilih karena donepezil adalah
pengikat DUAL-SITE: benzil di CAS (Trp86), indanon di PAS (Trp286). Artinya
konformasi gorge pada struktur ini sudah terbukti mampu menampung ligan yang
merentang kedua situs - persis yang dibutuhkan bis-ligan pada dataset ini.
Memakai struktur apo atau kompleks ligan CAS-saja berisiko gorge-nya terlalu
sempit, sehingga bis-ligan gagal terdocking bukan karena kimianya.

Langkah: ambil chain A -> buang heteroatom -> PDBFixer (H pada pH 7,4)
         -> PDBQT lewat Meeko.
"""
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import MDAnalysis as mda
from pdbfixer import PDBFixer
from openmm.app import PDBFile

ROOT = Path(__file__).resolve().parent.parent
STRUCT = ROOT / "structures"
DOCK = ROOT / "docking"
PH = 7.4
BOX_SIZE = 30.0          # A; ligan terpanjang berdiameter 31,6 A

ap = argparse.ArgumentParser()
ap.add_argument("--pdb", default="4ey7")
ap.add_argument("--chain", default="A")
ap.add_argument("--ligand-resname", default="E20")
a = ap.parse_args()
DOCK.mkdir(exist_ok=True)

src = STRUCT / f"{a.pdb}.pdb"
u = mda.Universe(str(src))

# ---- kotak docking dari geometri gorge -----------------------------
cas = u.select_atoms(f"chainID {a.chain} and resid 86 and not backbone")
pas = u.select_atoms(f"chainID {a.chain} and resid 286 and not backbone")
if cas.n_atoms == 0 or pas.n_atoms == 0:
    sys.exit("Trp86/Trp286 tidak ditemukan - periksa penomoran residu")
center = (cas.center_of_geometry() + pas.center_of_geometry()) / 2
print(f"jarak CAS-PAS : {np.linalg.norm(cas.center_of_geometry()-pas.center_of_geometry()):.2f} A")
print(f"pusat kotak   : {np.round(center,3)}")
print(f"ukuran kotak  : {BOX_SIZE} A")

# ---- ligan kristal, untuk validasi redocking -----------------------
lig = u.select_atoms(f"resname {a.ligand_resname} and chainID {a.chain}")
if lig.n_atoms:
    lig.write(str(DOCK / f"{a.pdb}_ligand_ref.pdb"))
    print(f"ligan referensi: {lig.n_atoms} atom -> {a.pdb}_ligand_ref.pdb")

# ---- reseptor: protein chain A saja, tanpa heteroatom --------------
# Air, NAG (glikosilasi), NO3 dibuang. Air gorge dihapus untuk docking;
# pada MD nanti sistem disolvasi ulang sehingga air kembali secara eksplisit.
# Pangkas terminus-C yang tidak teratur (resid 542, B-factor ~61, ~32 A dari
# situs aktif). PDBFixer menaruh OXT nyaris sama jauh dari C dan CB, sehingga
# perseptor ikatan berbasis jarak menyimpulkan karbon bervalensi 5. Residu ini
# tidak relevan secara kimia; dipangkas agar reseptor docking dan sistem MD
# nanti memakai rantai yang sama persis.
_all = u.select_atoms(f"protein and chainID {a.chain}")
first, last = int(_all.residues.resids.min()), int(_all.residues.resids.max())
rec = u.select_atoms(f"protein and chainID {a.chain} and resid {first}:{last-1}")
print(f"rantai dipangkas: resid {first}-{last-1} (residu {last} dibuang)")
raw = DOCK / f"{a.pdb}_rec_raw.pdb"
rec.write(str(raw))
print(f"\nreseptor mentah: {rec.n_atoms} atom, {rec.n_residues} residu")

fixer = PDBFixer(filename=str(raw))
fixer.findMissingResidues()
# jangan sambung celah di ujung rantai - hanya isi celah internal
keys = sorted(fixer.missingResidues.keys())
chains = list(fixer.topology.chains())
for k in list(keys):
    ch = chains[k[0]]
    if k[1] == 0 or k[1] == len(list(ch.residues())):
        del fixer.missingResidues[k]
if fixer.missingResidues:
    print(f"  celah internal diisi: {len(fixer.missingResidues)}")
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(PH)
clean = DOCK / f"{a.pdb}_rec_clean.pdb"
with open(clean, "w") as fh:
    PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
n_at = sum(1 for _ in fixer.topology.atoms())
print(f"reseptor bersih (H pada pH {PH}): {n_at} atom -> {clean.name}")

# ---- PDBQT lewat Meeko ---------------------------------------------
base = str(DOCK / f"{a.pdb}_rec")
cmd = [sys.executable, "/opt/anaconda3/envs/qsar-md/bin/mk_prepare_receptor.py",
       "--read_pdb", str(clean), "-o", base, "-p",
       "--box_size", str(BOX_SIZE), str(BOX_SIZE), str(BOX_SIZE),
       "--box_center", *[f"{c:.3f}" for c in center]]
r = subprocess.run(cmd, capture_output=True, text=True)
print(f"\nmk_prepare_receptor: {'OK' if r.returncode==0 else 'GAGAL'}")
if r.returncode != 0:
    print(r.stdout[-1500:]); print(r.stderr[-1500:])
else:
    for f in sorted(DOCK.glob(f"{a.pdb}_rec*.pdbqt")):
        print(f"  {f.name}  ({f.stat().st_size:,} B)")

# file kotak per-struktur: versi lama menimpa satu file bersama,
# sehingga menyiapkan reseptor kedua diam-diam merusak kotak yang pertama
np.save(DOCK / f"{a.pdb}_box_center.npy", center)
np.save(DOCK / "box_center.npy", center) if a.pdb == "4ey7" else None
with open(DOCK / f"{a.pdb}_box.txt", "w") as fh:
    fh.write(f"center {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n"
             f"size {BOX_SIZE} {BOX_SIZE} {BOX_SIZE}\n")
print(f"\nparameter kotak disimpan: docking/{a.pdb}_box.txt")
