#!/usr/bin/env python
"""Stage 11a — choose a receptor ensemble that spans the gorge's flexibility.

Docking into 4EY7 alone failed for ligands whose binding mode differs from
donepezil's: huperzine A landed 3.10 A from its crystal pose, and self-docking
it into its own structure gave 0.24 A. The cause was measured -- the Tyr337
side chain sits 2.53 A apart between 4EY5 and 4EY7 and blocks the true pose.

An ensemble is only useful if its members actually differ where it matters, so
structures are chosen by clustering on the side-chain positions of the
gorge-lining residues after superposition, not by resolution or by date.
Picking the N best-resolution structures would likely return N copies of the
same conformation.
"""
from __future__ import annotations
import json, sys, urllib.request, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
STRUCT = ROOT/"structures"
STRUCT.mkdir(exist_ok=True)
REF = "4ey7"
N_PICK = int(sys.argv[1]) if len(sys.argv) > 1 else 8
GORGE = [72, 74, 86, 121, 124, 286, 295, 297, 337, 338, 341]

def rcsb_search():
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity.pdbx_description",
                "operator": "contains_phrase", "value": "acetylcholinesterase"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entity_source_organism.scientific_name",
                "operator": "exact_match", "value": "Homo sapiens"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less", "value": 2.6}}]},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": 60},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                      "direction": "asc"}]}}
    r = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",
        data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    return [x["identifier"] for x in json.load(urllib.request.urlopen(r))["result_set"]]

ids = rcsb_search()
print(f"hAChE structures below 2.6 A: {len(ids)}", flush=True)

def fetch(p):
    f = STRUCT/f"{p.lower()}.pdb"
    if not f.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{p.upper()}.pdb", f)
        except Exception:
            return None
    return f

ref_u = mda.Universe(str(STRUCT/f"{REF}.pdb"))
def ca_one(u, chain="A"):
    sel = u.select_atoms(f"protein and chainID {chain} and name CA")
    seen, idx = set(), []
    for i, a in enumerate(sel):
        if a.resid not in seen:
            seen.add(a.resid); idx.append(i)
    return sel[idx]
ref_ca = ca_one(ref_u)

# side-chain coordinates of the gorge residues, in the 4EY7 frame
def gorge_vector(pdb_id):
    f = fetch(pdb_id)
    if f is None:
        return None
    try:
        u = mda.Universe(str(f))
        if "A" not in set(u.atoms.chainIDs):
            return None
        m_ca = ca_one(u)
        common = np.intersect1d(ref_ca.residues.resids, m_ca.residues.resids)
        if len(common) < 400:
            return None
        rp = ref_ca.positions[np.isin(ref_ca.resids, common)]
        mp = m_ca.positions[np.isin(m_ca.resids, common)]
        n = min(len(rp), len(mp)); rp, mp = rp[:n], mp[:n]
        rc, mc = rp.mean(0), mp.mean(0)
        R, ca_rmsd = mda_align.rotation_matrix(mp-mc, rp-rc)
        if ca_rmsd > 1.5:
            return None
        vec = []
        for res in GORGE:
            s = u.select_atoms(f"chainID A and resid {res} and not backbone "
                               f"and not name H*")
            t = ref_u.select_atoms(f"chainID A and resid {res} and not backbone "
                                   f"and not name H*")
            if s.n_atoms == 0 or s.n_atoms != t.n_atoms:
                return None
            vec.append(((s.positions - mc) @ R.T + rc).ravel())
        return np.concatenate(vec), float(ca_rmsd)
    except Exception:
        return None

vecs, keep = {}, []
for p in ids:
    r = gorge_vector(p)
    if r is not None:
        vecs[p], _ = r
        keep.append(p)
    print(f"  {p} {'ok' if r is not None else 'skipped'}", flush=True)
print(f"\nusable structures: {len(keep)}")

V = np.array([vecs[p] for p in keep])
n = len(keep)
D = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        d = np.sqrt(((V[i]-V[j])**2).reshape(-1, 3).sum(1).mean())
        D[i, j] = D[j, i] = d
print(f"gorge side-chain RMSD between structures: "
      f"median {np.median(D[np.triu_indices(n,1)]):.2f}, "
      f"max {D.max():.2f} A")

Z = linkage(squareform(D, checks=False), method="average")
labels = fcluster(Z, N_PICK, criterion="maxclust")
picked = []
for c in sorted(set(labels)):
    members = [i for i in range(n) if labels[i] == c]
    # the member closest to its cluster's centre represents it
    rep = min(members, key=lambda i: D[i, members].mean())
    picked.append(keep[rep])
if REF.upper() not in picked:
    picked[0] = REF.upper()      # keep the reference in the ensemble
print(f"\nensemble of {len(picked)}: {' '.join(picked)}")

sub = [keep.index(p) for p in picked]
print(f"spread within the ensemble: "
      f"median {np.median(D[np.ix_(sub,sub)][np.triu_indices(len(sub),1)]):.2f}, "
      f"max {D[np.ix_(sub,sub)].max():.2f} A")
(ROOT/"docking"/"ensemble.txt").write_text("\n".join(p.lower() for p in picked))
print(f"\nwritten: docking/ensemble.txt")
