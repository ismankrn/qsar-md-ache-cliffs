#!/usr/bin/env python
"""Table — the learning algorithms and the settings they were given.

Built by introspecting the pipelines in rung4_model.py rather than by transcribing them,
so the table cannot drift from the code that produced the results. For each family it
reports the scikit-learn estimator used on the regression framings and on the
classification framing, every setting that differs from the library default, and whether
the features reach the estimator standardised.

Usage:  python scripts/table_algorithms.py
"""
from __future__ import annotations
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/"scripts"))
import rung4_model as M

RES = ROOT/"results"
SKIP = {"n_jobs", "random_state", "verbose", "copy_X", "copy"}
# Reported even when it equals the library default, because the family is named for it.
ALWAYS = {"n_neighbors"}
PRETTY = {"alpha": "alpha", "n_components": "components", "n_neighbors": "neighbours",
          "C": "C", "gamma": "gamma", "n_estimators": "trees", "max_iter": "max.\\ iter",
          "probability": "probability"}


def fmt(v):
    if isinstance(v, float) and v == int(v):
        return str(int(v))
    if v is True:
        return "on"
    return str(v)


def dissect(factory):
    """(estimator class name, {setting: value}, standardised?) for one pipeline."""
    pipe = factory()
    est = pipe.named_steps["est"]
    base = type(est)()
    out = {}
    for k, v in est.get_params(deep=False).items():
        if k in SKIP:
            continue
        if k in ALWAYS or v != base.get_params(deep=False).get(k):
            out[PRETTY.get(k, k)] = fmt(v)
    return type(est).__name__, out, "sc" in pipe.named_steps


def settings(r, c):
    """One settings cell, naming a side only where the two differ."""
    rs = r[1] if r else {}
    cs = c[1] if c else {}
    if not r or not c:
        d = rs or cs
        return ", ".join(f"{k} $=$ {v}" for k, v in d.items()) or "library defaults"
    shared = {k: v for k, v in rs.items() if cs.get(k) == v}
    parts = [f"{k} $=$ {v}" for k, v in shared.items()]
    parts += [f"{k} $=$ {v} (reg.)" for k, v in rs.items() if k not in shared]
    parts += [f"{k} $=$ {v} (clf.)" for k, v in cs.items() if k not in shared]
    return ", ".join(parts) or "library defaults"


def main() -> None:
    order = ["Ridge (linear)", "Logistic regression", "PLS", "k-NN (k=5)",
             "SVM (RBF)", "Random Forest", "Gradient Boosting"]
    rows = []
    for fam in order:
        r = dissect(M.REG[fam]) if fam in M.REG else None
        c = dissect(M.CLF[fam]) if fam in M.CLF else None
        name = fam.replace("k-NN (k=5)", "$k$-nearest neighbours")

        # One estimator cell rather than two columns: where the two sides differ only by
        # the suffix, the shared stem is factored out, which is what makes the table fit.
        if r and c:
            rn, cn = r[0], c[0]
            stem = rn[:len(rn) - len("Regressor")] if rn.endswith("Regressor") else ""
            est = (f"{stem}" + r"\allowbreak" + r"\{Regressor, Classifier\}"
                   if stem and cn == stem + "Classifier" else f"{rn}, {cn}")
            where = "A, B, C"
        elif r:
            est, where = r[0], "A, B"
        else:
            est, where = c[0], "C"

        rows.append(" & ".join([name, est, where, settings(r, c),
                                "yes" if (r or c)[2] else "no"]) + r" \\")

    out = (
        r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{2.9cm}"
        r">{\raggedright\arraybackslash}p{4.5cm}"
        r"c"
        r">{\raggedright\arraybackslash}p{4.3cm}c@{}}" "\n"
        r"\toprule" "\n"
        r"Family & scikit-learn estimator & Framings & Settings & Scaled \\" "\n"
        r"\midrule" "\n" + "\n".join(rows) + "\n" r"\bottomrule" "\n"
        r"\end{tabular}" "\n")
    (RES/"39_algorithms.tex").write_text(out)
    for row in rows:
        print("   ", " | ".join(x.strip() for x in row.rstrip(" \\\\").split("&")))
    print(f"\n  seed {M.SEED}")
    print("  -> results/39_algorithms.tex")


if __name__ == "__main__":
    main()
