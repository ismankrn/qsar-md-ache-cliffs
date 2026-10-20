# Submitting to J. Chem. Inf. Model.

Everything below is built from this repository. Figures come from `scripts/fig_*.py`,
tables from `scripts/table_*.py`, and every number in the text traces to a file under
`results/` — see `results/INDEX.md` for which stage produced what.

## Build

```bash
cd manuscript
latexmk -pdf si.tex            # first: main.tex reads si.aux for the Table S numbers
latexmk -pdf main.tex
latexmk -pdf cover_letter.tex
```

## What to upload

| file | what it is |
|---|---|
| `manuscript/main.pdf` | the Article, double spaced with line numbers |
| `manuscript/si.pdf` | Supporting Information, Tables S1–S6 |
| `manuscript/cover_letter.pdf` | cover letter |
| `figures/fig_toc.pdf` | table-of-contents graphic, 3.25 × 1.75 in |

## Before submitting

**Rebuild `fig_thermo` on the cluster.** It is the one figure that cannot be
regenerated here, because it reads the GROMACS energy files rather than anything under
`results/`. Its script already carries the font size ACS requires; the committed PDF
still has the old one, which falls below the 4.5 pt lettering floor once published.

```bash
python scripts/fig_thermo.py --edr <directory of .xvg files>
```

**Turn off the drafting switch** in `manuscript/main.tex` if the journal does not want
line numbers and double spacing at submission. JCIM accepts both, so this is optional.

```latex
\draftingfalse   % was \draftingtrue
% and comment out \doublespacing and \linenumbers
```

## Conformance notes

Checked against the author guidelines (coden `jcisd8`) and the QSAR/QSPR editorial they
bind submissions to, Jorgensen, *J. Chem. Inf. Model.* **2006**, *46*, 937.

- Title is 9 words, under the preferred 12.
- Abstract is 250 words and states what is novel, as requirement 1 asks.
- All 271 molecular structures are in Table S1, as requirement 4 asks, and are also in
  the deposit.
- Every reported number is a pooled out-of-fold estimate on data not in the training
  fold, which is guideline 1. The sealed test split was never opened; Methods says why.
- Supporting Information is described in the nonsentence form with the file extension.
- A Data Availability Statement and a Funding statement are present, and the competing
  interest sentence is the one ACS prescribes.
- Figure lettering is at least 4.5 pt once scaled to the 7 in double-column maximum,
  for every figure except `fig_thermo` until it is rebuilt.

Guideline 7 discourages three things the paper touches, and the cover letter answers
each: an extensively modelled dataset, a high ratio of descriptors to data points, and
new descriptors without evidence of superiority.
