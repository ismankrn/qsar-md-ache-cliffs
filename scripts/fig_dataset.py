#!/usr/bin/env python
"""Figure 1 and Table 1 — the fingerprint that saturates, and the size confound.

Three panels, each answering a question a reader should ask before any model appears.
The target screen that once opened this figure was dropped, since the manuscript states
the criteria hAChE met and does not argue the choice against the other candidates.

(a) What binary ECFP4 does to the pairs in the set. Every cliff pair is constitutional,
    differing by 1 to 10 heavy atoms. Binary ECFP4 nevertheless reports Tanimoto exactly
    1.000 for 49 of the 97 -- molecules that differ by a median of four heavy atoms,
    declared identical. Count-based Tanimoto separates the same pairs at a median of
    0.878. Folding to bits discards the multiplicity that distinguishes them.

(b) The plane of the cliff definition, with similarity measured by the binary
    fingerprint that performed the selection. The empty band between the two groups is
    imposed by the criterion rather than by the chemistry.

(c) The size confound. The two members of a cliff pair differ by 4.4 heavy atoms on
    average against 1.8 for a control pair, so any size-correlated descriptor separates
    the groups for a reason unrelated to the cliff.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
TARGET = "CHEMBL220"
CLIFF_DELTA = 2.0


def esc(s) -> str:
    return str(s).replace("_", r"\_").replace("&", r"\&")


def main() -> None:
    feas = pd.read_csv(RES/"01_feasibility.csv")
    au = pd.read_csv(RES/"02_cliff_similarity_audit.csv")
    pr = pd.read_csv(RES/f"02_pairs_{TARGET}.csv")
    lg = pd.read_csv(RES/f"02_ligand_set_{TARGET}.csv")
    sp = pd.read_csv(RES/f"03_splits_ligands_{TARGET}.csv")
    spp = pd.read_csv(RES/f"03_splits_pairs_{TARGET}.csv")
    heavy = (pd.read_csv(RES/"08_rung2_descriptors.csv")
               .drop_duplicates("molecule_chembl_id")
               .set_index("molecule_chembl_id")["n_heavy"])

    pr = pr.assign(d_heavy=[abs(heavy.get(a, np.nan) - heavy.get(b, np.nan))
                            for a, b in zip(pr.ligand_a, pr.ligand_b)])

    FIGS.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1))

    # ---- (b) binary vs count Tanimoto ----
    ax = axes[0]
    s = ax.scatter(au.t_bit, au.t_count, c=au.d_heavy, cmap="viridis", s=26,
                   edgecolors="0.3", linewidths=0.4)
    n_sat = int((au.t_bit >= 0.9995).sum())
    ax.axvline(1.0, color="0.5", lw=0.8, ls=":")
    ax.annotate(f"{n_sat} of {len(au)} cliff pairs sit here:\n"
                f"binary Tanimoto = 1.000, a median of\n"
                f"{au[au.t_bit>=0.9995].d_heavy.median():.0f} heavy atoms apart",
                xy=(0.9995, au[au.t_bit >= 0.9995].t_count.quantile(0.25)),
                xytext=(0.902, 0.955), fontsize=8.2, ha="left", va="top",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.4",
                                connectionstyle="arc3,rad=-0.2"))
    plt.colorbar(s, ax=ax, label="heavy-atom difference", fraction=0.046, pad=0.03)
    ax.set_xlabel("Tanimoto, binary ECFP4", fontsize=8.5)
    ax.set_ylabel("Tanimoto, count ECFP4", fontsize=8.5)
    ax.text(0.0, 1.02, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    # ---- (c) the cliff definition ----
    ax = axes[1]
    for grp, col, lab in (("control", "#4393c3", "control"), ("cliff", "#b2182b", "cliff")):
        q = pr[pr.group == grp]
        ax.scatter(q.tanimoto, q.delta.abs(), s=22, color=col, alpha=0.75,
                   edgecolors="none", label=f"{lab} (n={len(q)})")
    ax.axhline(CLIFF_DELTA, color="0.35", lw=1.0, ls="--")
    ax.text(0.99, CLIFF_DELTA, f"|$\\Delta$pChEMBL| = {CLIFF_DELTA:g} ",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=8.2, va="bottom",
            color="0.35")
    # Verified identical to t_bit in the audit file: the pair-selection criterion is the
    # binary fingerprint, which panel (a) shows saturates. So the measure that chose these
    # pairs is the one that cannot tell molecules eight heavy atoms apart from identical.
    ax.set_xlabel("Tanimoto similarity, binary ECFP4 (the selection criterion)",
                  fontsize=8.5)
    ax.set_ylabel(r"$|\Delta$pChEMBL$|$", fontsize=8.5)
    ax.text(0.0, 1.02, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
    ax.legend(fontsize=8.2, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.98))

    # ---- (d) the size confound ----
    ax = axes[2]
    data = [pr[pr.group == g].d_heavy.dropna() for g in ("cliff", "control")]
    parts = ax.violinplot(data, positions=[0, 1], showmedians=True, widths=0.7)
    for pc, col in zip(parts["bodies"], ("#b2182b", "#4393c3")):
        pc.set_facecolor(col); pc.set_alpha(0.55)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        parts[key].set_color("0.3"); parts[key].set_linewidth(1.0)
    for i, d in enumerate(data):
        ax.text(i, d.max()+0.5, f"mean {d.mean():.1f}", ha="center", fontsize=8.2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["cliff", "control"], fontsize=8.5)
    ax.set_ylabel("Heavy-atom difference within a pair", fontsize=8.5)
    ax.text(0.0, 1.02, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    for ax in axes.ravel():
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=(0, 0.01, 1, 0.97))
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig1_dataset.{ext}", dpi=300)
    print("  -> figures/fig1_dataset.pdf / .png")

    # ---------------- Table 1 ----------------
    m = lg.merge(sp[["molecule_chembl_id", "split"]], on="molecule_chembl_id", how="left")
    rows = []
    for split in ("dev", "test"):
        ml = m[m.split == split]
        mp = spp[spp.split == split]
        rows.append(dict(split=split, ligands=len(ml),
                         cliff_ligands=int((ml.group == "cliff").sum()),
                         pairs=len(mp), cliff_pairs=int((mp.group == "cliff").sum()),
                         control_pairs=int((mp.group != "cliff").sum())))
    tot = dict(split="total", ligands=len(m), cliff_ligands=int((m.group == "cliff").sum()),
               pairs=len(spp), cliff_pairs=int((spp.group == "cliff").sum()),
               control_pairs=int((spp.group != "cliff").sum()))
    rows.append(tot)
    t1 = pd.DataFrame(rows)
    t1.to_csv(RES/"23_dataset_composition.csv", index=False)

    tex = [r"\begin{tabular}{lrrrrr}", r"\toprule",
           r"Split & Ligands & of which cliff & Pairs & Cliff pairs & Control pairs \\",
           r"\midrule"]
    for _, r in t1.iterrows():
        name = {"dev": "Development", "test": "Test (sealed)", "total": "Total"}[r.split]
        line = (f"{name} & {r.ligands} & {r.cliff_ligands} & {r.pairs} & "
                f"{r.cliff_pairs} & {r.control_pairs}" + r" \\")
        if r.split == "total":
            tex.append(r"\midrule")
        tex.append(line)
    tex += [r"\bottomrule", r"\end{tabular}"]
    (RES/"23_dataset_composition.tex").write_text("\n".join(tex) + "\n")
    print("  -> results/23_dataset_composition.csv / .tex")
    print(t1.to_string(index=False))

    print(f"\n  median |dpChEMBL|: cliff {pr[pr.group=='cliff'].delta.abs().median():.2f}, "
          f"control {pr[pr.group!='cliff'].delta.abs().median():.2f}")
    print(f"  mean heavy-atom difference: cliff {data[0].mean():.1f}, "
          f"control {data[1].mean():.1f}")
    print(f"  binary ECFP4 saturation: {n_sat} of {len(au)} cliff pairs at Tanimoto 1.000; "
          f"all {len(au)} are constitutional")


if __name__ == "__main__":
    main()
