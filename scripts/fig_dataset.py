#!/usr/bin/env python
"""Figure 1 and Table 1 — target selection, the fingerprint that saturates, and the size
confound.

Four panels, each answering a question a reader should ask before any model appears.

(a) Why this target. Six ChEMBL targets were screened. Activity cliffs split into two kinds:
    constitutional pairs, which differ in composition, and stereochemical pairs, which are
    stereoisomers. Only constitutional cliffs are a fair test of a 2D descriptor, because no
    stereo-blind fingerprint can resolve a stereoisomer pair by construction while any 3D
    method resolves it trivially -- a benchmark dominated by them is rigged in favour of
    structure. AChE is the only target screened where constitutional cliffs dominate, 100
    against 10; BACE1 inverts that at 5 against 101.

(b) What binary ECFP4 does to the pairs that survive. Every cliff pair in the final set is
    constitutional, differing by 1 to 10 heavy atoms. Binary ECFP4 nevertheless reports
    Tanimoto exactly 1.000 for 49 of the 97 -- molecules that differ by a median of four
    heavy atoms, declared identical. Count-based Tanimoto separates the same pairs at a
    median of 0.878. Folding to bits discards the multiplicity that distinguishes them.

(c) The cliff definition on the selected set, and (d) the size confound that follows from it:
    cliff partners differ by more heavy atoms than control partners, so any size-correlated
    descriptor separates the groups for the wrong reason. Stage 10 found 0 of 19 descriptors
    survive a size adjustment.

Also writes Table 1 (dataset composition) as a LaTeX fragment.

Usage:
    python scripts/fig_dataset.py
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
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2))

    # ---- (a) target screen ----
    ax = axes[0, 0]
    f = feas.sort_values("cliffs_constitutional", ascending=True)
    yy = np.arange(len(f))
    ax.barh(yy - 0.19, f.cliffs_constitutional, height=0.36, color="#2166ac",
            label="constitutional")
    ax.barh(yy + 0.19, f.cliffs_stereochemical, height=0.36, color="#d6604d",
            label="stereochemical")
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{n}\n{t}" for n, t in zip(f.pref_name.str.slice(0, 26),
                                                    f.target_chembl_id)], fontsize=6.5)
    for i, (v, ok) in enumerate(zip(f.cliffs_constitutional, f.verdict)):
        ax.text(max(f.cliffs_stereochemical.max(), f.cliffs_constitutional.max())*1.02, i,
                ok, fontsize=7, va="center",
                color="#1b7837" if ok == "PASS" else "0.55")
    ax.set_xlabel("Activity-cliff pairs found", fontsize=8.5)
    ax.text(0.0, 1.02, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
    ax.legend(fontsize=7, frameon=False, loc="lower right")

    # ---- (b) binary vs count Tanimoto ----
    ax = axes[0, 1]
    s = ax.scatter(au.t_bit, au.t_count, c=au.d_heavy, cmap="viridis", s=26,
                   edgecolors="0.3", linewidths=0.4)
    n_sat = int((au.t_bit >= 0.9995).sum())
    ax.axvline(1.0, color="0.5", lw=0.8, ls=":")
    ax.annotate(f"{n_sat} of {len(au)} cliff pairs sit here:\n"
                f"binary Tanimoto = 1.000, a median of\n"
                f"{au[au.t_bit>=0.9995].d_heavy.median():.0f} heavy atoms apart",
                xy=(0.9995, au[au.t_bit >= 0.9995].t_count.quantile(0.25)),
                xytext=(0.902, 0.955), fontsize=7, ha="left", va="top",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.4",
                                connectionstyle="arc3,rad=-0.2"))
    plt.colorbar(s, ax=ax, label="heavy-atom difference", fraction=0.046, pad=0.03)
    ax.set_xlabel("Tanimoto, binary ECFP4", fontsize=8.5)
    ax.set_ylabel("Tanimoto, count ECFP4", fontsize=8.5)
    ax.text(0.0, 1.02, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    # ---- (c) the cliff definition ----
    ax = axes[1, 0]
    for grp, col, lab in (("control", "#4393c3", "control"), ("cliff", "#b2182b", "cliff")):
        q = pr[pr.group == grp]
        ax.scatter(q.tanimoto, q.delta.abs(), s=22, color=col, alpha=0.75,
                   edgecolors="none", label=f"{lab} (n={len(q)})")
    ax.axhline(CLIFF_DELTA, color="0.35", lw=1.0, ls="--")
    ax.text(0.99, CLIFF_DELTA, f"|$\\Delta$pChEMBL| = {CLIFF_DELTA:g} ",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=7, va="bottom",
            color="0.35")
    # Verified identical to t_bit in the audit file: the pair-selection criterion is the
    # binary fingerprint, which panel (b) shows saturates. So the measure that chose these
    # pairs is the one that cannot tell molecules eight heavy atoms apart from identical.
    ax.set_xlabel("Tanimoto similarity, binary ECFP4 (the selection criterion)",
                  fontsize=8.5)
    ax.set_ylabel(r"$|\Delta$pChEMBL$|$", fontsize=8.5)
    ax.text(0.0, 1.02, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.98))

    # ---- (d) the size confound ----
    ax = axes[1, 1]
    data = [pr[pr.group == g].d_heavy.dropna() for g in ("cliff", "control")]
    parts = ax.violinplot(data, positions=[0, 1], showmedians=True, widths=0.7)
    for pc, col in zip(parts["bodies"], ("#b2182b", "#4393c3")):
        pc.set_facecolor(col); pc.set_alpha(0.55)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        parts[key].set_color("0.3"); parts[key].set_linewidth(1.0)
    for i, d in enumerate(data):
        ax.text(i, d.max()+0.5, f"mean {d.mean():.1f}", ha="center", fontsize=7.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["cliff", "control"], fontsize=8.5)
    ax.set_ylabel("Heavy-atom difference within a pair", fontsize=8.5)
    ax.text(0.0, 1.02, "d", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    for ax in axes.ravel():
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=(0, 0.01, 1, 0.96))
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
