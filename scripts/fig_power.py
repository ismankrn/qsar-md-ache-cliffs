#!/usr/bin/env python
"""Figure 3 — what this study could have detected, and what caps it.

A null result means nothing without the effect size the design could have seen. Three panels,
all from stage 17, all reading CSVs the stage-17 scripts now write so the figure and the text
cite the same numbers rather than two independent calculations that drift apart.

(a) Power against a rung-1 baseline of 0.200, from `power_proper.py`. This is the corrected
    analysis: the baseline is re-simulated inside every replicate, so rung 1's own sampling
    uncertainty is carried. Two earlier versions compared against a fixed 0.200 and reported
    power that was too high; `power_and_ceiling.py` still contains one of them, and its power
    table is deliberately not used here.

(b) Width of the 95% interval on the gap, which is flat near 0.30 across the whole range.
    That is the portable sentence: this study separates two R² values about 0.30 apart and
    not less, independent of where on the scale they sit.

(c) The measurement-noise ceiling against where the rungs actually land. Even on the
    pessimistic p90 spread the ceiling is 0.62, far above every rung, so what the ladder
    fails to explain is real signal and not assay noise. That matters: it is the difference
    between "nothing to find" and "not found by these descriptors".

Usage:
    python scripts/fig_power.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
TARGET_POWER = 0.80


def need(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"missing {path.relative_to(ROOT)} — run the stage-17 script that writes it")
    return pd.read_csv(path)


def main() -> None:
    pw = need(RES/"17_power_curve.csv").sort_values("candidate_r2")
    nc = need(RES/"17_noise_ceiling.csv")
    ladder = pd.read_csv(RES/"20_ladder.csv") if (RES/"20_ladder.csv").exists() else None
    base = float(pw.baseline_r2.iloc[0])

    x = pw.candidate_r2.to_numpy(float)
    power = pw.power.to_numpy(float)
    # Linear interpolation on the simulated grid, not a fitted curve: six points cannot
    # support a parametric form, and drawing one would imply precision the design lacks.
    hit = np.nan
    if power.max() >= TARGET_POWER and power.min() < TARGET_POWER:
        hit = float(np.interp(TARGET_POWER, power, x))

    FIGS.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9))

    # ---- (a) power ----
    ax = axes[0]
    ax.plot(x, 100*power, "o-", color="#2166ac", ms=5, lw=1.3)
    ax.axhline(100*TARGET_POWER, color="0.6", lw=0.9, ls=":")
    ax.axvline(base, color="#b2182b", lw=1.0, ls="--")
    ax.text(base, 4, f" rung 1 = {base:.3f}", color="#b2182b", fontsize=7,
            rotation=90, va="bottom")
    if np.isfinite(hit):
        ax.annotate(f"{TARGET_POWER:.0%} power\nat $R^2$ ≈ {hit:.2f}", xy=(hit, 100*TARGET_POWER),
                    xytext=(hit-0.02, 46), fontsize=7.5, ha="right",
                    arrowprops=dict(arrowstyle="->", lw=0.8, color="0.4"))
    if ladder is not None:
        r4 = ladder[(ladder.framing == "B") & ladder.key.eq("r4")].value
        if len(r4):
            ax.axvline(float(r4.iloc[0]), color="#1b7837", lw=1.0, ls="-.")
            ax.text(float(r4.iloc[0]), 96, " rung 4 observed", color="#1b7837", fontsize=7,
                    rotation=90, va="top")
    ax.set_xlabel("Candidate descriptor set, true $R^2$", fontsize=8.5)
    ax.set_ylabel("Power to call it better than rung 1 (%)", fontsize=8.5)
    ax.text(0.0, 1.02, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
    ax.set_ylim(-3, 103)

    # ---- (b) interval width ----
    ax = axes[1]
    ax.plot(x, pw.ci_width.to_numpy(float), "s-", color="#762a83", ms=5, lw=1.3)
    ax.axhline(float(pw.ci_width.mean()), color="0.6", lw=0.9, ls=":")
    ax.text(0.98, float(pw.ci_width.mean()), f"mean {pw.ci_width.mean():.2f} ",
            transform=ax.get_yaxis_transform(), fontsize=7.5, color="0.4",
            va="bottom", ha="right")
    ax.set_xlabel("Candidate descriptor set, true $R^2$", fontsize=8.5)
    ax.set_ylabel("95% interval width on the gap", fontsize=8.5)
    ax.text(0.0, 1.02, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
    ax.set_ylim(0, max(0.45, float(pw.ci_width.max())*1.25))

    # ---- (c) ceiling vs where the rungs land ----
    ax = axes[2]
    pess = nc[nc.sigma_kind == "p90"].iloc[0]
    ax.axhspan(float(pess.ceiling_r2), 1.0, color="0.92", zorder=0)
    ax.text(0.5, (float(pess.ceiling_r2)+1.0)/2, "unreachable:\nassay noise", ha="center",
            va="center", fontsize=7.5, color="0.45")
    for _, r in nc.iterrows():
        ax.axhline(float(r.ceiling_r2), color="0.55", lw=0.8, ls=":")
        # The line sits at the CEILING, not at sigma. Labelling it with sigma alone made
        # the smallest spread produce the highest line, which reads as an error. The
        # median and mean ceilings are 0.04 apart, so their labels are pushed clear of
        # each other rather than left to overlap.
        dy = {"median": 0.035, "mean": -0.040}.get(r.sigma_kind, 0.0)
        ax.text(1.02, float(r.ceiling_r2) + dy,
                f"σ {r.sigma_kind} = {r.sigma:.2f}  →  ceiling {r.ceiling_r2:.3f}",
                fontsize=5.8, color="0.45", va="center",
                transform=ax.get_yaxis_transform())
    if ladder is not None:
        b = ladder[ladder.framing == "B"]
        xs = np.arange(len(b))
        ax.bar(xs, b.value.to_numpy(float), color="#2166ac", width=0.6)
        ax.set_xticks(xs)
        # "Rung 2 pose" and "Rung 2 ensemble" both start "Rung 2", so keep the qualifier
        ax.set_xticklabels([s.replace("Rung ", "R").replace(" ", "\n", 1) for s in b.rung],
                           fontsize=6.8)
        ax.set_xlabel("Rung (framing B)", fontsize=8.5)
    ax.set_ylabel("$R^2$", fontsize=8.5)
    ax.set_ylim(-0.1, 1.0)
    ax.text(0.0, 1.02, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")

    for ax in axes:
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    n = int(pw.n_pairs.iloc[0])
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig3_power.{ext}", dpi=300)
    print(f"  -> figures/fig3_power.pdf / .png")
    if np.isfinite(hit):
        print(f"  {TARGET_POWER:.0%} power needs R2 ~ {hit:.3f}; rung 1 baseline {base:.3f}")
    print(f"  interval width mean {pw.ci_width.mean():.3f}")
    print(f"  pessimistic ceiling {float(pess.ceiling_r2):.3f} (sigma p90 = {pess.sigma:.3f})")


if __name__ == "__main__":
    main()
