# Manuscript

    latexmk -pdf main.tex     # build
    latexmk -c                # clean aux files

Figures resolve from `../figures` via `\graphicspath`; generated tables are
`\input` from `../results`. Nothing here should contain a number typed by hand:
if a value appears in the text, a CSV in `../results` should produce it, and
`../results/INDEX.md` says which stage wrote which file.

`\draftingtrue` in `main.tex` turns on line numbers and makes `\TODO` red. Set it
false before submission and both disappear — but the `\TODO`s must be resolved
first, not hidden.

## Artefacts already wired in

| float | file | source script |
|---|---|---|
| Figure 2 (ladder) | `../figures/fig2_ladder.pdf` | `scripts/fig_ladder.py` |
| Figure 3 (power) | `../figures/fig3_power.pdf` | `scripts/fig_power.py` |
| Figure 5 (duration) | `../figures/fig5_duration.pdf` | `scripts/fig_duration.py` |
| Figure 1 (hAChE structure) | `../figures/Figure1_hAChE_structure.pdf` | `figures/fig1a.pml`, `fig1b.pml` + `scripts/compose_figures.py` |
| Figure 2 (redocking) | `../figures/Figure2_redocking_validation.pdf` | `figures/fig2*.pml` + `scripts/compose_figures.py` |
| Figure 3 (dataset) | `../figures/fig1_dataset.pdf` | `scripts/fig_dataset.py` |
| Table 1 (composition) | `../results/23_dataset_composition.tex` | `scripts/fig_dataset.py` |
| Figure 4 (reliability) | `../figures/fig4_reliability.pdf` | `scripts/fig_reliability.py` |
| Table 2 (benchmark) | `../results/22_benchmark_main.tex` | `scripts/table_benchmark.py` |
| Tables S1, S2 | `../results/22_benchmark_si_{reg,clf}.tex` | `scripts/table_benchmark.py` |
| Table 3 (subset) | `../results/21_subset_inflation.tex` | `scripts/table_subset_inflation.py` |

## Not yet built

- Figure 2 intervals — rerun `fig_ladder.py --n-boot 300`
- `refs.bib` — every entry is an unverified stub

## A note on the two PyMOL figures

Figures 1 and 2 are rendered from the `.pml` scripts in `figures/` and composed by
`scripts/compose_figures.py`. **PyMOL is not installed in this environment**, so the
panels cannot be regenerated here — the committed PDFs and PNGs are the artefacts.
Keep them. Their file names predate the `figN_` convention used by the newer
scripts; the LaTeX numbering is by order of appearance, not by file name.
