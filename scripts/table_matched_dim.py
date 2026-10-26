"""Stage 29 — the matched-dimensionality control, as a table.

The headline comparison gives rung 1 thousands of features and rung 4 twelve, so the
obvious objection is that the fingerprint wins on width rather than on content. The
control cuts the fingerprint to the width of rung 4 with a univariate filter fitted
inside each training fold, and the result is easier to read as a grid than as prose:
three framings, three numbers each, and the question of whether the ordering survives.

Every value is read from the benchmark tables that rung4_model.py writes, not retyped,
so the table cannot drift from the numbers the models actually produced. The selector
sits inside the pipeline, so the k features are chosen from training data alone; doing
it on the full matrix first would let the held-out fold pick the features it is then
asked to predict from, which inflates the control and is invisible in the score.

Outputs results/29_matched_dim.csv and .tex.
"""
import re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT/"results"

# framing -> (file, metric column, the arm label the generator wrote for each row)
SRC = {
    "A": ("19_rung4_framingA_w20.csv", "r2"),
    "B": ("19_rung4_benchmark_w20.csv", "r2_B"),
    "C": ("19_rung4_framingC_w20.csv", "auc_roc"),
}

def best(df, metric, pred):
    sub = df[df.descriptors.map(pred)]
    if sub.empty:
        raise SystemExit(f"arm tidak ditemukan untuk {metric}")
    return float(sub[metric].max()), int(sub.loc[sub[metric].idxmax(), "n_features"])

rows = []
for fr, (fname, metric) in SRC.items():
    d = pd.read_csv(RES/fname)
    # rung 1 as it would be used: the better of its two feature sets
    full, _ = best(d, metric, lambda s: s.startswith("rung1") and "k=" not in s)
    # rung 1 cut to the width of rung 4, selected inside each training fold.
    # The width is read from the arm's own label and not from n_features, which records
    # the matrix the selector was given rather than the k it kept.
    cut, _ = best(d, metric, lambda s: s.startswith("rung1") and "k=" in s)
    label = d[d.descriptors.str.contains("k=", na=False)].descriptors.iloc[0]
    width = int(re.search(r"k=(\d+)", label).group(1))
    # rung 4 as reported in the ladder: the better of filtered and unfiltered
    r4, _ = best(d, metric, lambda s: s.startswith("rung4"))
    rows.append({"framing": fr, "width": width, "rung1_full": full,
                 "rung1_cut": cut, "rung4": r4,
                 "ordering_reversed": bool(cut < r4 <= full)})

df = pd.DataFrame(rows)
df.to_csv(RES/"29_matched_dim.csv", index=False)

# The resolvable framing is A, which is where the only development interval excludes
# zero. Naming it in the table is the whole point: the two framings whose ordering the
# control reverses are the two the design cannot resolve.
RESOLVED = {"A": "yes", "B": "no", "C": "no"}
fmt = lambda fr, v: f"${v:+.3f}$" if fr in "AB" else f"${v:.3f}$"

lines = [r"\begin{tabular}{ccrrrcc}", r"\hline",
         r"Framing & Width & Rung 1, full & Rung 1, cut & Rung 4 & Ordering "
         r"& Resolvable \\", r"\hline"]
for r in df.itertuples():
    lines.append(
        f"{r.framing} & {r.width} & {fmt(r.framing, r.rung1_full)} & "
        f"{fmt(r.framing, r.rung1_cut)} & {fmt(r.framing, r.rung4)} & "
        f"{'reversed' if r.ordering_reversed else 'unchanged'} & "
        f"{RESOLVED[r.framing]} \\\\")
lines += [r"\hline", r"\end{tabular}"]
(RES/"29_matched_dim.tex").write_text("\n".join(lines) + "\n")

print(df.to_string(index=False))
print(f"\n-> results/29_matched_dim.csv")
print(f"-> results/29_matched_dim.tex")
