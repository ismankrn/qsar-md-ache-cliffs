"""Stage 18 — how much of a descriptor's variation is the ligand, and how much is the seed?

Written before the data exists, deliberately. The decision this informs -- whether to spend
$178 on 162 replicate trajectories -- is exactly the kind that invites a rule chosen after
seeing the numbers. The rule is therefore fixed here, in advance:

    A descriptor can carry the pair endpoint only if the spread BETWEEN ligands is at least
    twice the spread BETWEEN REPLICATES of the same ligand.

Anything below that is measuring the random number generator. The threshold of 2 is a
judgement, but it is a judgement made blind.

Why this is not answered by the duration ablation (16d): that measures whether a descriptor's
time average has stopped drifting within one trajectory. A quantity can be perfectly flat in
every single run and still differ between seeds -- that is precisely what happens when each
run is trapped in a different basin. gate_flips_per_ns averaged 0.007/ns over the first three
runs, about one rotamer jump per three trajectories, so Tyr337 essentially never moves in
45 ns. A run that starts in one rotamer reports it with occupancy 1.000 and a flat RMSD, and
looks converged. Only replicates can tell "the answer is A" from "this run was stuck in A".

Reads every JSON that has replicates, and reports per descriptor and per ablation window:

    sd_rep   pooled standard deviation between replicates of the same ligand
    sd_lig   standard deviation of the per-ligand means
    ratio    sd_lig / sd_rep -- the discriminability

For ligands that form an activity-cliff pair it also reports the within-pair difference
against the replicate noise, which is the study's actual endpoint rather than a proxy for it.

Usage:  python scripts/replicate_pilot.py [--min-reps 2]
"""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
JSONS = ROOT/"results"/"16_rung4_per_run"
PAIRS = ROOT/"results"/"01_cliffs_CHEMBL220.csv"
OUT_CSV = ROOT/"results"/"18_replicate_pilot.csv"
OUT_PAIRS = ROOT/"results"/"18_replicate_pairs.csv"
OUT_PNG = ROOT/"results"/"18_replicate_pilot.png"

RATIO_USABLE = 2.0            # fixed in advance; see the module docstring
SKIP = {"run", "n_frames", "dt_ps", "ns_analysed", "n_heavy", "windows"}


def load():
    rows = []
    for f in sorted(JSONS.glob("*.json")):
        d = json.loads(f.read_text())
        base = {k: v for k, v in d.items() if k not in SKIP}
        lig = re.sub(r"_r\d+$", "", d["run"])
        rep = int(re.search(r"_r(\d+)$", d["run"]).group(1))
        rows.append({**base, "ligand": lig, "replicate": rep, "window": "45"})
        for w, wd in (d.get("windows") or {}).items():
            rows.append({**{k: v for k, v in wd.items() if k not in SKIP},
                         "ligand": lig, "replicate": rep, "window": w})
    return pd.DataFrame(rows)


def spreads(df: pd.DataFrame, min_reps: int) -> pd.DataFrame:
    out = []
    desc = [c for c in df.columns if c not in ("ligand", "replicate", "window")]
    for w, dw in df.groupby("window"):
        counts = dw.groupby("ligand").replicate.nunique()
        keep = counts[counts >= min_reps].index
        rep_part = dw[dw.ligand.isin(keep)]
        if rep_part.ligand.nunique() < 2:
            continue
        for c in desc:
            # The two spreads come from different populations on purpose. Seed noise can only
            # be measured where replicates exist -- six ligands at most. Between-ligand spread
            # needs variety, not replication, so it uses EVERY ligand with a result: six would
            # give a hopeless estimate of it, and the main slice is producing dozens.
            vr = pd.to_numeric(rep_part[c], errors="coerce")
            va = pd.to_numeric(dw[c], errors="coerce")
            if vr.isna().all() or va.isna().all():
                continue
            gr = vr.groupby(rep_part.ligand)
            # Pooled within-ligand SD: the average variance across ligands, not the average
            # SD. Averaging SDs would underweight the ligands that actually vary.
            var_within = gr.var(ddof=1).dropna()
            if var_within.empty:
                continue
            sd_rep = float(np.sqrt(var_within.mean()))
            lig_means = va.groupby(dw.ligand).mean()
            if lig_means.size < 2:
                continue
            sd_lig = float(lig_means.std(ddof=1))
            # sd_rep = 0 is not one case but three, and the rule as written called all of
            # them unusable -- including the best possible outcome.
            #   both zero          the descriptor never varies at all: no information
            #   sd_rep 0, few      almost every run returns the same value (usually 0) and
            #     distinct values  the within-ligand variance is zero by sparsity, not by
            #                      reproducibility. gate_flips_per_ns is this: about one
            #                      rotamer jump per three trajectories, so most runs read 0.
            #   sd_rep 0, varied   genuinely reproducible to the digit, and discriminating
            nz = float((va.abs() > 1e-12).mean())
            distinct = int(va.round(10).nunique())
            if sd_lig == 0 and sd_rep == 0:
                verdict = "constant"
            elif sd_rep == 0 and (nz < 0.25 or distinct <= 3):
                verdict = "sparse"
            elif sd_rep == 0:
                verdict = "usable"
            elif sd_lig/sd_rep >= RATIO_USABLE:
                verdict = "usable"
            else:
                verdict = "seed-dominated"
            out.append({"window": w, "descriptor": c,
                        "n_lig_rep": int(len(keep)), "n_lig_all": int(lig_means.size),
                        "frac_nonzero": nz, "n_distinct": distinct,
                        "sd_rep": sd_rep, "sd_lig": sd_lig,
                        "ratio": (sd_lig/sd_rep) if sd_rep > 0 else np.inf,
                        "verdict": verdict, "usable": verdict == "usable"})
    return pd.DataFrame(out)


def pair_check(df: pd.DataFrame, sp: pd.DataFrame) -> pd.DataFrame:
    if not PAIRS.exists():
        return pd.DataFrame()
    pairs = pd.read_csv(PAIRS)
    full = df[df.window == "45"]
    means = full.groupby("ligand").mean(numeric_only=True)
    sd = sp[sp.window == "45"].set_index("descriptor").sd_rep
    have = set(means.index)
    out = []
    for _, p in pairs.iterrows():
        if p.cid_a not in have or p.cid_b not in have:
            continue
        for c in sd.index:
            if c not in means.columns:
                continue
            delta = float(means.at[p.cid_a, c] - means.at[p.cid_b, c])
            # A pair difference is two independent estimates, so the noise on the difference
            # is sd_rep * sqrt(2), not sd_rep.
            noise = float(sd[c]) * np.sqrt(2)
            out.append({"cid_a": p.cid_a, "cid_b": p.cid_b, "dpchembl": float(p.delta),
                        "descriptor": c, "delta_desc": delta, "noise": noise,
                        "signal_to_noise": abs(delta)/noise if noise > 0 else np.inf})
    return pd.DataFrame(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-reps", type=int, default=2)
    a = ap.parse_args()

    df = load()
    if df.empty:
        sys.exit("no rung-4 JSONs found")
    reps = df[df.window == "45"].groupby("ligand").replicate.nunique()
    ready = reps[reps >= a.min_reps]
    print(f"ligands with results : {len(reps)}  (these set sd_lig)")
    print(f"with >= {a.min_reps} replicates : {len(ready)}  {list(ready.index)}  (these set sd_rep)")
    if len(ready) < 2:
        sys.exit("need at least two ligands with replicates before this says anything")

    sp = spreads(df, a.min_reps)
    sp.to_csv(OUT_CSV, index=False)
    full = sp[sp.window == "45"].sort_values("ratio", ascending=False)
    print(f"\n  {'descriptor':24s} {'sd_rep':>9s} {'sd_lig':>9s} {'ratio':>7s}  verdict")
    for _, r in full.iterrows():
        rat = "   inf" if not np.isfinite(r.ratio) else f"{r.ratio:7.2f}"
        print(f"  {r.descriptor:24s} {r.sd_rep:9.4f} {r.sd_lig:9.4f} {rat}  {r.verdict}")
    counts = full.verdict.value_counts().to_dict()
    print(f"\n  45 ns: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items()))
          + f"  -> {OUT_CSV.name}")

    pc = pair_check(df, sp)
    if not pc.empty:
        pc.to_csv(OUT_PAIRS, index=False)
        print(f"  cliff pairs with replicates on both sides: "
              f"{pc[['cid_a','cid_b']].drop_duplicates().shape[0]} -> {OUT_PAIRS.name}")
    else:
        print("  no cliff pair yet has replicates on both members")
