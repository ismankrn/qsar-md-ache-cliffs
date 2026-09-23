#!/usr/bin/env python
"""Rung 2, refined — a bridging measure that actually measures bridging.

The first attempt used "nearest ligand atom within 4.5 A of each tryptophan
side chain". That criterion is nearly vacuous: the nearest atoms of Trp86 and
Trp286 are only 11.37 A apart, so a ligand spanning ~2.4 A satisfies it, and
92% of the set passed. The 15.2 A figure is centroid to centroid, not the gap
a ligand has to cross.

Replaced by a projection onto the gorge axis. The axis runs from the Trp86
side-chain centroid to the Trp286 side-chain centroid; every ligand heavy atom
is projected onto it, and the span covered is measured directly. A ligand that
engages both sites covers most of the axis; one that sits at the floor covers
little. This is what "spanning the gorge" means geometrically.

The hypothesis is also tested where it was always aimed: at homologation
cliffs, where the partners differ only in linker length. Testing it on
substituent cliffs, where the mechanism is something else entirely, only
dilutes it. Both strata are reported.
"""
from __future__ import annotations
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import MDAnalysis as mda
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import TanimotoSimilarity
from meeko import PDBQTMolecule, RDKitMolCreate
from scipy.stats import mannwhitneyu, spearmanr

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parent.parent
RES, POSES = ROOT/"results", ROOT/"docking"/"poses"
TARGET = "CHEMBL220"

u = mda.Universe(str(ROOT/"docking"/"4ey7_rec_clean.pdb"))
cas_c = u.select_atoms("resid 86 and not backbone and not name H*").center_of_geometry()
pas_c = u.select_atoms("resid 286 and not backbone and not name H*").center_of_geometry()
axis = pas_c - cas_c
L = np.linalg.norm(axis)
unit = axis / L
print(f"gorge axis Trp86 -> Trp286 : {L:.2f} A\n")


def axis_span(xyz):
    """How much of the CAS->PAS axis the ligand covers, as a fraction."""
    t = (xyz - cas_c) @ unit          # projection, 0 at CAS, L at PAS
    lo, hi = np.clip([t.min(), t.max()], 0, L)
    return {
        "axis_span_frac": float(max(hi - lo, 0) / L),
        "axis_reach_pas": float(np.clip(t.max(), 0, L) / L),   # how far up it gets
        "axis_foot_cas":  float(np.clip(t.min(), 0, L) / L),   # how low it sits
        "d_cas_centroid": float(np.linalg.norm(xyz - cas_c, axis=1).min()),
        "d_pas_centroid": float(np.linalg.norm(xyz - pas_c, axis=1).min()),
    }


score = pd.read_csv(RES/"07_docking_scores.csv")
rows = []
for cid in sorted(score.molecule_chembl_id.unique()):
    per = []
    for seed in (1, 2, 3):
        f = POSES/f"{cid}_seed{seed}.pdbqt"
        if not f.exists():
            continue
        try:
            m = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(
                PDBQTMolecule.from_file(str(f), skip_typing=True))[0])
            per.append(axis_span(np.array(m.GetConformer(0).GetPositions())))
        except Exception:
            pass
    if per:
        df = pd.DataFrame(per)
        rows.append({"molecule_chembl_id": cid,
                     **{c: df[c].median() for c in df.columns}})
sp = pd.DataFrame(rows).set_index("molecule_chembl_id")
sp.reset_index().to_csv(RES/"08_rung2_bridging.csv", index=False)
print(f"axis descriptors: {len(sp)} ligands")
print(f"  axis_span_frac: median {sp.axis_span_frac.median():.2f}, "
      f"range {sp.axis_span_frac.min():.2f}-{sp.axis_span_frac.max():.2f}\n")

# ---- classify each cliff pair: homologation or substituent -----------
lig = pd.concat([
    pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")[["molecule_chembl_id","smiles_std"]],
    pd.read_csv(RES/f"04_validation_ligands_{TARGET}.csv")[["molecule_chembl_id","smiles_std"]],
]).drop_duplicates("molecule_chembl_id").set_index("molecule_chembl_id")
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fp = {c: gen.GetFingerprint(Chem.MolFromSmiles(s))
      for c, s in lig.smiles_std.items() if Chem.MolFromSmiles(s)}

pairs = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
rec = []
for _, r in pairs.iterrows():
    a, b = r.ligand_a, r.ligand_b
    if a not in sp.index or b not in sp.index or a not in fp or b not in fp:
        continue
    # homologation cliff: binary ECFP4 cannot tell them apart, because Morgan
    # fingerprints are blind to chain length
    homolog = TanimotoSimilarity(fp[a], fp[b]) >= 0.9999
    rec.append({"group": r.group, "dp": r.delta,
                "kind": "homologation" if homolog else "substituent",
                "d_span": abs(sp.loc[a,"axis_span_frac"] - sp.loc[b,"axis_span_frac"]),
                "d_reach": abs(sp.loc[a,"axis_reach_pas"] - sp.loc[b,"axis_reach_pas"]),
                "d_pas": abs(sp.loc[a,"d_pas_centroid"] - sp.loc[b,"d_pas_centroid"])})
x = pd.DataFrame(rec)

print("DIFFERENCE WITHIN A PAIR, by stratum")
print(f"{'stratum':>26} {'n':>4} {'d_span':>8} {'d_reach':>8} {'d_pas':>8}")
print("-"*60)
strata = [("cliff, homologation", (x.group=="cliff")&(x.kind=="homologation")),
          ("cliff, substituent",  (x.group=="cliff")&(x.kind=="substituent")),
          ("control, homologation",(x.group=="control")&(x.kind=="homologation")),
          ("control, substituent", (x.group=="control")&(x.kind=="substituent"))]
for lab, mask in strata:
    s = x[mask]
    if len(s):
        print(f"{lab:>26} {len(s):>4} {s.d_span.mean():>8.3f} "
              f"{s.d_reach.mean():>8.3f} {s.d_pas.mean():>8.2f}")

print("\nTHE TEST — homologation cliffs vs their matched controls")
ch = x[(x.group=="cliff") & (x.kind=="homologation")]
co = x[(x.group=="control")]
for col, name in [("d_span","axis span"), ("d_reach","reach toward PAS"),
                  ("d_pas","distance to PAS")]:
    if len(ch) > 2 and len(co) > 2:
        pv = mannwhitneyu(ch[col].dropna(), co[col].dropna()).pvalue
        print(f"  {name:20} cliff {ch[col].mean():.3f}  control {co[col].mean():.3f}"
              f"   p = {pv:.4f}{'  *' if pv < 0.05 else ''}")
if len(ch) > 2:
    print(f"\n  within homologation cliffs, rho(|dp|, d_reach) = "
          f"{spearmanr(ch.dp, ch.d_reach).statistic:+.3f}  (n={len(ch)})")
x.to_csv(RES/"08_rung2_bridging_pairs.csv", index=False)
print(f"\nwritten: 08_rung2_bridging.csv, 08_rung2_bridging_pairs.csv")
