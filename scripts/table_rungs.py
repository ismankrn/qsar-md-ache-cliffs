#!/usr/bin/env python
"""Table — the four descriptor sets, what each one needs, and what each one costs.

Section 2.6 describes the rungs one paragraph at a time, which is the right order for
reading and the wrong one for comparing. This puts the input each rung requires, what
its descriptors capture, how many features it carries and what it cost to produce
beside each other, so that the ordering by cost which gives the section its title is
visible rather than asserted.

Both the feature counts and the costs are measured rather than typed.

The feature counts come from the ladder output, i.e. the number of columns each rung
was actually scored with, which is not the number of columns in its descriptor file.
The model drops the per-seed spread columns, so counting the file gives 36 for rung 2
where 34 were scored.

The costs come from the docking logs, which record seconds per run. The molecular
dynamics figure is derived from the throughput measured on this campaign's hardware,
which is stated as wall-clock on one GPU rather than in core-hours, since a
GPU-resident mdrun is not billed in cores and converting it would invent precision.

Writes results/32_rungs.tex.

Usage:
    python scripts/table_rungs.py
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, Descriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"

MD_NS = 25.0            # production length per ligand, Section 2.7
MD_NS_PER_DAY = 33.61   # measured on this campaign's nodes at 16 CPU cores per run


def feature_counts() -> dict:
    """What each rung was scored with, from the ladder output."""
    d = pd.read_csv(RES/"20_ladder.csv")
    a = d[d.framing == "A"].set_index("rung")
    n = {k: int(a.loc[v, "n_features"]) for k, v in
         (("r1_bit", "Rung 1 2D"), ("r2", "Rung 2 pose"),
          ("r2e", "Rung 2 ensemble"), ("r3", "Rung 3 MM/GBSA"), ("r4", "Rung 4 MD"))}
    # Framing C doubles every feature, since a pair row is the sum and the absolute
    # difference concatenated, so the count-based arm's own width is half of its entry.
    c = d[(d.framing == "C") & (d.rung == "Rung 1 2D")].iloc[0]
    assert "count" in c["arm"], c["arm"]
    n["r1_cnt"] = int(c["n_features"]) // 2
    n["r4_filt"] = int(d[(d.framing == "B") & (d.rung == "Rung 4 MD")].iloc[0]["n_features"])
    return n


def costs() -> dict:
    """Median seconds per ligand, from the runs that were actually timed."""
    dock = pd.read_csv(RES/"07_docking_scores.csv")
    ens = pd.read_csv(RES/"11_ensemble_docking.csv")
    return {
        "dock_s": float(dock.groupby("molecule_chembl_id").seconds.sum().median()),
        "dock_seeds": int(dock.groupby("molecule_chembl_id").seed.nunique().median()),
        "ens_s": float(ens.groupby("molecule_chembl_id").seconds.sum().median()),
        "ens_rec": int(ens.receptor.nunique()),
        "md_h": MD_NS/MD_NS_PER_DAY*24.0,
        "fp_ms": rung1_ms(),
    }


def rung1_ms() -> float:
    """Median milliseconds to build both rung-1 feature sets for one ligand.

    Timed here rather than asserted as "under a second". It is the only cost in the
    table small enough that an assertion would have passed unchallenged, which is the
    reason to measure it.
    """
    lig = pd.read_csv(RES/"02_ligand_set_CHEMBL220.csv")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ts = []
    for s in lig.smiles_std:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        t0 = time.perf_counter()
        gen.GetFingerprintAsNumPy(m)
        gen.GetCountFingerprintAsNumPy(m)
        Descriptors.CalcMolDescriptors(m)
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))*1000.0


def main() -> None:
    n, c = feature_counts(), costs()
    for k, v in {**n, **c}.items():
        print(f"    {k:12s} {v}")

    rows = [
        ("1", "Molecular graph",
         "ECFP4 substructure fingerprint, with physicochemical descriptors in the "
         "second feature set",
         f"{n['r1_bit']:,} or {n['r1_cnt']:,}".replace(",", "{,}"),
         f"{c['fp_ms']:.0f} ms"),
        ("2", "One docked pose in 4EY7",
         "Distances to the two sites and to the gate, bridging along the axis between "
         "them, cation--$\\pi$ contacts, per-residue contacts, ligand shape and size",
         str(n["r2"]), f"{c['dock_s']:.0f} s, {c['dock_seeds']} seeds"),
        ("2e", f"Docked poses in {c['ens_rec']} receptor conformations",
         "As rung 2 without the axis descriptors, reduced to the median over the "
         "ensemble",
         str(n["r2e"]), f"{c['ens_s']:.0f} s, {c['ens_rec']} receptors"),
        ("3", "The rung-2 pose",
         "MM/GBSA total with its molecular-mechanics, polar and non-polar components, "
         "buried surface area, per-residue decomposition",
         str(n["r3"]), "rung 2, plus one energy evaluation"),
        ("4", f"A {MD_NS:.0f} ns trajectory",
         "Contact occupancy at each end of the gorge, position along the axis, "
         "behaviour of the Tyr337 gate, fluctuation of ligand and residues",
         f"{n['r4']}, or {n['r4_filt']} filtered",
         f"$\\approx${c['md_h']:.0f} h on one GPU"),
    ]
    col = r">{\raggedright\arraybackslash}p"
    tex = [r"\begin{tabular}{l" + col + "{2.6cm}" + col + "{5.2cm}" + col + "{1.7cm}"
           + col + r"{2.4cm}}", r"\toprule",
           r"Rung & Input required & What the descriptors capture & Features "
           r"& Cost per ligand \\", r"\midrule"]
    tex += [" & ".join(r) + r" \\[2pt]" for r in rows]
    tex += [r"\bottomrule", r"\end{tabular}"]
    (RES/"32_rungs.tex").write_text("\n".join(tex) + "\n")
    print("  -> results/32_rungs.tex")


if __name__ == "__main__":
    main()
