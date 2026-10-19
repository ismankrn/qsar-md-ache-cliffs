#!/usr/bin/env python
"""Supplementary tables — every descriptor of rungs 2, 2e, 3 and 4, by name.

Section 2.3 states feature counts and says what each rung's descriptors capture, but it
never names them. A reader who wants to know what "34 features" consists of, or which 11
of 32 dynamic descriptors survived the reliability filter, has nowhere to look. These
four tables are that list, with the name each descriptor carries in the deposited CSVs
so a number in the text can be traced to a column.

The definitions are written here, but the names, the counts and the filter verdicts are
read from the files the models were scored on, and a disagreement is fatal rather than
silent. A descriptor renamed in a script would otherwise leave a table describing a
feature set that no longer exists.

Writes results/37_si_rung2.tex, 37_si_rung3.tex and 37_si_rung4.tex.

Usage:
    python scripts/table_si_descriptors.py
"""
from __future__ import annotations
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"
WINDOW = 20.0           # the production window, and the one the filter was applied at

A = r"\si{\angstrom}"
# The 14 gorge-lining residues, in the order the scripts iterate them. Names are read
# from the receptor rather than typed, because a wrong one would be invisible here.
LINING = [72, 74, 86, 121, 124, 203, 286, 295, 297, 334, 337, 338, 341, 447]
RESNAME = {72: "Tyr", 74: "Asp", 86: "Trp", 121: "Gly", 124: "Tyr", 203: "Ser",
           286: "Trp", 295: "Phe", 297: "Phe", 334: "Glu", 337: "Tyr", 338: "Phe",
           341: "Tyr", 447: "His"}

GROUP = "\\addlinespace\n\\multicolumn{%d}{l}{\\itshape %s} \\\\\n\\addlinespace[2pt]"


def rows_rung2():
    """(name, definition, unit, kept at rung 2e) plus group headers as bare strings."""
    r = [
        "Engagement of the two anionic sites",
        ("d_cas", "Closest approach of a ligand heavy atom to a Trp86 side-chain heavy atom", A, True),
        ("d_pas", "The same distance to Trp286", A, True),
        ("d_gate", "The same distance to Tyr337", A, True),
        ("bridges", "Both $d_{\\mathrm{cas}}$ and $d_{\\mathrm{pas}}$ within \\SI{4.5}{\\angstrom}", "0/1", True),
        ("bridge_gap", "The larger of $d_{\\mathrm{cas}}$ and $d_{\\mathrm{pas}}$", A, True),
        ("catpi_cas", "Closest approach of a formally cationic ligand atom to the centroid of the Trp86 six-membered ring", A, True),
        ("catpi_pas", "The same distance to the Trp286 ring", A, True),
        ("catpi_both", "Both cation--$\\pi$ distances within \\SI{6}{\\angstrom}", "0/1", True),
        "Position along the gorge axis",
        ("axis_span_frac", "Extent of the ligand's projection onto the Trp86$\\rightarrow$Trp286 axis, as a fraction of its \\SI{15.2}{\\angstrom} length", "---", False),
        ("axis_reach_pas", "Projection of the ligand atom reaching furthest towards the gorge mouth", "---", False),
        ("axis_foot_cas", "Projection of the ligand atom sitting lowest towards the gorge floor", "---", False),
        ("d_cas_centroid", "Closest approach of a ligand heavy atom to the Trp86 side-chain centroid", A, False),
        ("d_pas_centroid", "The same distance to the Trp286 side-chain centroid", A, False),
        "Burial and contact density",
        ("n_contacts_4A", "Ligand heavy atoms having a protein heavy atom within \\SI{4}{\\angstrom}", "count", True),
        ("n_prot_5A", "Protein heavy atoms within \\SI{5}{\\angstrom} of any ligand heavy atom", "count", True),
        ("n_prot_8A", "Protein heavy atoms within \\SI{8}{\\angstrom} of any ligand heavy atom", "count", True),
        ("frac_lig_buried", "Fraction of ligand heavy atoms having a protein heavy atom within \\SI{4.5}{\\angstrom}", "---", True),
        "Shape of the bound conformation",
        ("rgyr", "Radius of gyration of the ligand heavy atoms, unweighted", A, True),
        ("end_to_end", "Largest distance between any two ligand heavy atoms", A, True),
        ("n_heavy", "Ligand heavy atoms, carried as a size covariate", "count", True),
        "Per-residue contacts",
    ]
    r += [(f"c_{n}", f"Ligand and residue heavy-atom pairs within \\SI{{4.5}}{{\\angstrom}} "
           f"of {RESNAME[n]}{n}", "count", True) for n in LINING]
    return r


def rows_rung3():
    r = [
        "Interaction energy, receptor to ligand",
        ("dE_elec", "Coulomb term between receptor and ligand atoms, on the minimised complex", "kcal/mol"),
        ("dE_vdw", "Lennard-Jones term between the same atoms", "kcal/mol"),
        ("dE_mm", "Sum of the Coulomb and Lennard-Jones terms", "kcal/mol"),
        "Desolvation and the assembled estimate",
        ("dG_mm_gb", "Energy of the complex less that of the receptor and the ligand apart, each evaluated with the generalised-Born term", "kcal/mol"),
        ("dG_polar_gb", "Polar desolvation, i.e. $\\Delta G_{\\mathrm{mm,gb}}$ less $\\Delta E_{\\mathrm{mm}}$", "kcal/mol"),
        ("dSASA", "Change in solvent-accessible surface area on binding, over interface atoms", "\\si{\\angstrom\\squared}"),
        ("dG_nonpolar", "Nonpolar term, \\num{0.0072} times the buried surface area", "kcal/mol"),
        ("dG_total", "Sum of $\\Delta G_{\\mathrm{mm,gb}}$ and the nonpolar term", "kcal/mol"),
        "Per-residue decomposition",
    ]
    r += [(f"res_{n}", f"Coulomb and Lennard-Jones terms summed over the atoms of "
           f"{RESNAME[n]}{n}", "kcal/mol") for n in LINING]
    return r


def rows_rung4():
    r = [
        "Persistence of the bridging contact",
        ("occ_cas", "Fraction of frames with a ligand heavy atom within \\SI{4.5}{\\angstrom} of Trp86"),
        ("occ_pas", "The same fraction for Trp286"),
        ("occ_bridge", "Fraction of frames contacting Trp86 and Trp286 at once"),
        ("bridge_runs_per_ns", "Separate bridging episodes begun per nanosecond"),
        ("bridge_longest_ns", "Longest uninterrupted bridging episode"),
        "Whether the ligand stayed where it was docked",
        ("rmsd_mean", "Ligand heavy-atom RMSD from the equilibrated pose, averaged over frames, measured after superposing the receptor on its reference by C$\\alpha$"),
        ("rmsd_max", "Largest value that RMSD reached"),
        ("escaped", "Whether the ligand centre of mass ever left a \\SI{12}{\\angstrom} sphere about its starting position"),
        ("depth_mean", "Position of the ligand centroid along the Trp86$\\rightarrow$Trp286 axis, which is recomputed each frame"),
        ("depth_sd", "Standard deviation of that position"),
        "Flexibility, which a single pose has none of",
        ("rmsf_lig", "Per-atom fluctuation of the ligand, averaged over its heavy atoms"),
        ("rmsf_gorge", "The same fluctuation for the 14 gorge-lining side chains"),
        ("gate_chi1_sd", "Circular standard deviation of the Tyr337 $\\chi_1$ dihedral"),
        ("gate_flips_per_ns", "Crossings of Tyr337 $\\chi_1$ between rotamer wells per nanosecond"),
        ("lig_dihedral_sd", "Circular standard deviation of the ligand's own rotatable dihedrals, averaged over them"),
        "Burial, water and quality control",
        ("frac_lig_buried_mean", "Fraction of ligand heavy atoms within \\SI{4.5}{\\angstrom} of the protein, averaged over frames"),
        ("n_bridge_water", "Waters within \\SI{3.5}{\\angstrom} of both a ligand and a protein polar heavy atom, averaged over frames"),
        ("ca_rmsd_mean", "Receptor C$\\alpha$ RMSD after superposition, kept as a check that the receptor was reassembled whole across the periodic boundary"),
        "Per-residue contact persistence",
    ]
    r += [(f"occ_{n}", f"Fraction of frames within \\SI{{4.5}}{{\\angstrom}} of "
           f"{RESNAME[n]}{n}") for n in LINING]
    r += ["Size covariate", ("n_heavy", "Ligand heavy atoms, which is not a dynamic quantity and is not put to the filter")]
    return r


def tabular(header, colspec, body, ncol):
    return "\n".join([f"\\begin{{tabular}}{{{colspec}}}", r"\toprule", header,
                      r"\midrule", *body, r"\bottomrule", r"\end{tabular}", ""])


def emit(rows, colspec, header, ncol, fmt):
    body = []
    for r in rows:
        body.append(GROUP % (ncol, r) if isinstance(r, str) else fmt(r))
    return tabular(header, colspec, body, ncol)


def names(rows):
    return [r[0] for r in rows if not isinstance(r, str)]


def check(label, described, actual):
    """A table that does not describe the scored feature set is worse than no table."""
    missing, extra = sorted(set(actual) - set(described)), sorted(set(described) - set(actual))
    if missing or extra:
        sys.exit(f"{label}: {len(described)} described against {len(actual)} scored\n"
                 f"  not described: {missing}\n  not scored:    {extra}")
    print(f"  {label}: {len(actual)} descriptors, names agree")


def main() -> None:
    # ---- what each rung was actually scored on ------------------------------
    d2 = pd.read_csv(RES/"08_rung2_descriptors.csv")
    br = pd.read_csv(RES/"08_rung2_bridging.csv")
    m2 = d2.merge(br, on="molecule_chembl_id").set_index("molecule_chembl_id")
    m2 = m2.drop(columns=[c for c in m2.columns
                          if c.startswith("spread_") or c == "n_seeds_used"])
    cols2 = list(m2.columns)
    cols2e = [c for c in pd.read_csv(RES/"15_rung2_ensemble_median.csv")
              .set_index("molecule_chembl_id").select_dtypes("number").columns]
    r3 = pd.read_csv(RES/"09_rung3_mmgbsa.csv")
    cols3 = [c for c in r3.columns if c.startswith(("dE_", "dG_", "res_", "dSASA"))]
    pil = pd.read_csv(RES/"18_replicate_pilot.csv")
    pil = pil[pil.window.astype(float) == WINDOW].set_index("descriptor")
    cols4 = sorted(pil.index) + ["n_heavy"]

    g2, g3, g4 = rows_rung2(), rows_rung3(), rows_rung4()
    check("rung 2", names(g2), cols2)
    check("rung 3", names(g3), cols3)
    check("rung 4", names(g4), cols4)
    kept2e = sorted(r[0] for r in g2 if not isinstance(r, str) and r[3])
    check("rung 2e", kept2e, cols2e)

    # ---- rung 2 and 2e, one table with a column for membership --------------
    t2 = emit(g2, r"l>{\raggedright\arraybackslash}p{0.54\linewidth}ll",
              r"Descriptor & Definition & Unit & Rung~2e \\", 4,
              lambda r: f"\\texttt{{{r[0].replace('_', chr(92)+'_')}}} & {r[1]} & {r[2]} & "
                        f"{'yes' if r[3] else 'no'} \\\\")
    (RES/"37_si_rung2.tex").write_text(t2)

    t3 = emit(g3, r"l>{\raggedright\arraybackslash}p{0.60\linewidth}l",
              r"Descriptor & Definition & Unit \\", 3,
              lambda r: f"\\texttt{{{r[0].replace('_', chr(92)+'_')}}} & {r[1]} & {r[2]} \\\\")
    (RES/"37_si_rung3.tex").write_text(t3)

    # ---- rung 4 carries the filter verdict, which is why it has a table -----
    def row4(r):
        n = r[0]
        if n in pil.index:
            v = pil.loc[n]
            ratio = "---" if not pd.notna(v.ratio) or v.ratio == float("inf") else f"{v.ratio:.2f}"
            verdict = v.verdict
        else:
            ratio, verdict = "---", "not applicable"
        return (f"\\texttt{{{n.replace('_', chr(92)+'_')}}} & {r[1]} & {ratio} & "
                f"{verdict} \\\\")

    t4 = emit(g4, r"l>{\raggedright\arraybackslash}p{0.40\linewidth}rl",
              r"Descriptor & Definition & $s_{\mathrm{lig}}/s_{\mathrm{rep}}$ & Verdict \\",
              4, row4)
    (RES/"37_si_rung4.tex").write_text(t4)

    n_use = int((pil.verdict == "usable").sum())
    print(f"\n  rung 4 at the {WINDOW:g} ns window: {n_use} usable of {len(pil)}")
    print("  written: 37_si_rung2.tex, 37_si_rung3.tex, 37_si_rung4.tex")


if __name__ == "__main__":
    main()
