# Manuscript

    latexmk -pdf main.tex     # build
    latexmk -c                # clean aux files

Figures resolve from `../figures` via `\graphicspath`; generated tables are
`\input` from `../results`. Nothing here should contain a number typed by hand:
if a value appears in the text, a CSV in `../results` should produce it, and
`../results/INDEX.md` says which stage wrote which file.

Results and Discussion are one section: each interpretation sits with the number it
rests on, the cross-cutting argument is gathered in "What this study tested, and what
it did not", and Limitations closes the section. `sections/04_discussion.tex` is gone.

`\draftingtrue` in `main.tex` turns on line numbers and makes `\TODO` red. Set it
false before submission and both disappear — but the `\TODO`s must be resolved
first, not hidden.

## Style rules

Six rules came out of drafting Sections 2.1 and 2.2 against the author's own papers in
`writing-references/`. Run `python scripts/check_manuscript.py` before committing a
drafted section; it exits non-zero if an enforced rule fails.

1. A drafted section carries no `\TODO` and no `\NOTE`. Guidance that is worth the
   reader's attention belongs in the prose; the rest is dropped.
2. Paragraphs are not left short.
3. No colons in prose or in captions. Headings are exempt, since a heading is a title
   rather than a sentence, but the drafted ones avoid them anyway.
4. The words "rung" and "ladder" do not appear before the subsection that defines them.
   Descriptor sets are named by what they are, as in "the pose-geometry descriptors".
5. A caption identifies what each panel shows. The numbers and the interpretation go in
   the prose, so a figure caption over roughly 120 words is a sign that something
   belongs elsewhere. A table caption is held to 25 words, i.e. a title and nothing
   more, because the column headings already say what each column holds and anything
   past that is prose in the wrong place.
6. Every figure panel is referenced explicitly in the prose, as `Figure~\ref{...}(a)`.
   All three figures in Sections 2.1 and 2.2 were first drafted referencing the figure
   only as a whole. Every panel's content was discussed, so nothing was missing from the
   argument, and a reader still had no way to tell which panel to look at.

7. Procedure in Methods is written in the past tense. The reference papers are split
   on this, with the 2018 BCSJ paper fully past-tense passive and the 2019 Kanazawa
   paper present-tense, and the author settled it in favour of the past. Three kinds
   of present tense stay. A pointer to the paper's own apparatus, as in "is shown in
   Figure~2"; a timeless fact about the protein or the chemistry, as in "Trp86 is
   located at the catalytic anionic site"; and a general claim that is not a step
   taken here, as in "no stereo-blind fingerprint can resolve such a pair".

Voice is matched to those papers: past-tense passive for procedure, "by using X",
"with the PDB ID being", "as shown in Figure~N", "i.e." for enumeration, and a roadmap
sentence opening each section. Register is matched; the non-native constructions in the
sources are not reproduced.

## Artefacts already wired in

| float | file | source script |
|---|---|---|
| Figure 2 (ladder) | `../figures/fig2_ladder.pdf` | `scripts/fig_ladder.py` |
| Figure 6 (power) | `../figures/fig3_power.pdf` | `scripts/fig_power.py` |
| Figure 8 (duration) | `../figures/fig5_duration.pdf` | `scripts/fig_duration.py` |
| Figure 1 (hAChE structure) | `../figures/Figure1_hAChE_structure.pdf` | `figures/fig1a.pml`, `fig1b.pml` + `scripts/compose_figures.py` |
| Figure 2 (redocking) | `../figures/Figure2_redocking_validation.pdf` | `figures/fig2*.pml` + `scripts/compose_figures.py` |
| Table 2 (system composition) | `../results/27_system_composition.tex` | `scripts/table_system_composition.py` |
| Figure 4 (thermodynamics) | `../figures/fig_thermo.pdf` | `scripts/fig_thermo.py` |
| Figure 5 (MD stability) | `../figures/fig_md_stability.pdf` | `scripts/fig_md_stability.py` |
| Figure 3 (interactions) | `../figures/fig_interactions.pdf` | `scripts/fig_interactions.py` |
| Figure 4 (dataset) | `../figures/fig1_dataset.pdf` | `scripts/fig_dataset.py` |
| Table 1 (composition) | `../results/23_dataset_composition.tex` | `scripts/fig_dataset.py` |
| Figure 7 (reliability) | `../figures/fig4_reliability.pdf` | `scripts/fig_reliability.py` |
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

`prolif` was added to the environment for Figure 3 (additive: prolif, dill,
multiprocess; nothing downgraded). PLIP was not used because its visual output needs
PyMOL, which is not installed here.

## Data pulled from the cluster for these figures

`table_system_composition.py` and `fig_thermo.py` read files that live only on BRIN
and are not in this repository:

* 217 `md/systems/*/build.json` — system composition. Pull with
  `ssh mahameru 'cd qsar-md && tar cf - md/systems/*/build.json' > builds.tar`.
* 64 usable `md/runs/*/prod.edr` — thermodynamic terms, extracted on the cluster
  with `gmx energy` selecting Potential, Temperature, Pressure and Density **by
  name**, because the term numbers differ between systems. 16 further energy files
  cover 0.4 ns and are the cancelled trembesi[75-78] jobs, not production runs.

The extracted `.xvg` files are 15 MB and are not committed; regenerate them with the
recipe above if the figure needs rebuilding.
