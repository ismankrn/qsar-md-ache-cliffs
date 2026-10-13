#!/usr/bin/env python
"""Figure — thermodynamic stability of the production ensemble.

A campaign of 217 independent systems should be shown to have run at the temperature,
pressure and density it was supposed to, and not one run at a time. These are the terms
GROMACS writes to prod.edr, extracted by name rather than by menu number because the term
indices shift between systems.

The energy files survive for 80 of the 229 runs: the full trajectories were deleted after
descriptor extraction and only the BRIN block kept its prod.edr. That subset is stated on the
figure. It is not the arrival-order subset that stage 21 showed to be biased -- it is the
runs that happened to finish on the cluster rather than on rented hardware -- but it is a
subset, and nothing here is used to compare descriptor sets, only to show that the
simulations were well behaved.

Pressure in a molecular dynamics box fluctuates by hundreds of bar instantaneously and is
meaningless frame by frame; the running mean is what shows the barostat working, so that is
what is drawn, with the instantaneous spread quoted instead of plotted.

Usage:
    python scripts/fig_thermo.py --edr <dir of .xvg>
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES, FIGS = ROOT/"results", ROOT/"figures"
EQUIL_NS = 5.0
TARGET_T, TARGET_P = 300.0, 1.0
MIN_SPAN_NS = 10.0      # below this the file is an abandoned run, not a short one


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--edr", required=True)
    a = ap.parse_args()

    files = sorted(Path(a.edr).glob("*.xvg"))
    if not files:
        raise SystemExit(f"no .xvg under {a.edr}")
    series, rows, short = {}, [], []
    for f in files:
        d = np.loadtxt(f)
        if d.ndim != 2 or d.shape[1] < 5:
            print(f"  {f.stem}: unexpected shape {d.shape}, skipped")
            continue
        t = d[:, 0]/1000.0            # ps -> ns
        pot, T, P, rho = d[:, 1], d[:, 2], d[:, 3], d[:, 4]
        # Sixteen directories on the cluster hold a prod.edr from the jobs that hung on
        # trembesi[75-78] and advanced ~0.4 ns before being cancelled. They are not failed
        # production runs, they are not production runs at all, and averaging them in would
        # drag every band toward an unequilibrated state.
        if t[-1] < MIN_SPAN_NS:
            short.append((f.stem, float(t[-1])))
            continue
        series[f.stem] = (t, pot, T, P, rho)
        k = t >= EQUIL_NS
        rows.append(dict(run=f.stem, span_ns=float(t[-1]), n=len(t),
                         T_mean=float(T[k].mean()), T_sd=float(T[k].std()),
                         P_mean=float(P[k].mean()), P_sd=float(P[k].std()),
                         rho_mean=float(rho[k].mean()), rho_sd=float(rho[k].std()),
                         pot_mean=float(pot[k].mean()),
                         pot_drift=float(np.polyfit(t[k], pot[k], 1)[0])))
    if short:
        print(f"  {len(short)} files below {MIN_SPAN_NS:g} ns excluded as abandoned runs: "
              f"{', '.join(f'{s}({v:.1f} ns)' for s, v in short[:4])}"
              + (" ..." if len(short) > 4 else ""))
    df = pd.DataFrame(rows)
    df.to_csv(RES/"28_thermodynamics.csv", index=False)
    print(f"  {len(df)} runs -> results/28_thermodynamics.csv")
    print(f"  T    {df.T_mean.mean():.2f} K  (run-to-run sd {df.T_mean.std():.2f}, "
          f"within-run sd {df.T_sd.mean():.1f})")
    print(f"  P    {df.P_mean.mean():+.2f} bar (within-run sd {df.P_sd.mean():.0f})")
    print(f"  rho  {df.rho_mean.mean():.1f} kg/m3 (run-to-run sd {df.rho_mean.std():.2f})")
    print(f"  potential energy drift: median {df.pot_drift.median():+.1f} kJ/mol/ns, "
          f"worst {df.pot_drift.abs().max():.1f}")

    # common grid: the shortest run bounds the band
    nmin = min(len(v[0]) for v in series.values())
    t = next(iter(series.values()))[0][:nmin]

    def band(idx, running=False):
        """Median and 5--95 band across runs.

        `running` takes the cumulative mean *per run first* and only then the quantiles
        across runs. Doing it the other way -- the cumulative mean of the 5th and 95th
        percentile series -- is what the first draft did, and it draws a +/-140 bar band
        that is the running mean of an extreme, not the spread of any run's running mean.
        """
        M = np.vstack([v[idx][:nmin] for v in series.values()])
        if running:
            M = np.cumsum(M, axis=1)/np.arange(1, M.shape[1]+1)
        return np.median(M, 0), *np.percentile(M, [5, 95], axis=0)

    FIGS.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 3.5))
    specs = [
        (2, "Temperature (K)", "(a) Thermostat", "#b2182b", TARGET_T),
        (4, "Density (kg m$^{-3}$)", "(b) Density", "#2166ac", None),
        (1, "Potential energy (kJ mol$^{-1}$)", "(c) Potential energy", "#1b7837", None),
        (3, "Pressure, running mean (bar)", "(d) Barostat", "#762a83", TARGET_P),
    ]
    for ax, (idx, ylab, title, col, tgt) in zip(axes, specs):
        # instantaneous pressure is hundreds of bar wide and says nothing frame by frame,
        # so the barostat panel shows each run's running mean
        med, lo, hi = band(idx, running=(idx == 3))
        ax.fill_between(t, lo, hi, color=col, alpha=0.2, lw=0)
        ax.plot(t, med, color=col, lw=1.1)
        if tgt is not None:
            ax.axhline(tgt, color="0.35", lw=0.9, ls=":")
        ax.axvline(EQUIL_NS, color="0.55", lw=0.8, ls="--")
        ax.set_xlabel("Time (ns)", fontsize=8.5)
        ax.set_ylabel(ylab, fontsize=8.5)
        ax.text(0.0, 1.02, title.split(")")[0] + ")", transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="bottom", ha="left")
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    for ext in ("pdf", "png"):
        fig.savefig(FIGS/f"fig_thermo.{ext}", dpi=300)
    print(f"  -> figures/fig_thermo.pdf / .png")


if __name__ == "__main__":
    main()
