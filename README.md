# Do dynamic descriptors explain activity cliffs that 2D-QSAR cannot?

An end-to-end QSAR → molecular dynamics study on human acetylcholinesterase
(hAChE, ChEMBL target CHEMBL220).

**Hypothesis.** Activity cliffs — near-identical molecules with wildly different
potency — are the canonical failure mode of QSAR. For hAChE, half of them are
*homologation cliffs*: pairs differing only in linker length, where the linker
determines whether a bis-cationic ligand can bridge the catalytic anionic site
(Trp86) and the peripheral anionic site (Trp286), 15.2 Å apart at either end of a
narrow ~20 Å gorge. That is a geometric and dynamic property. Descriptors derived
from molecular dynamics should capture it; 2D descriptors provably cannot.

## Status

| Stage | State |
|---|---|
| 01 Target screening (6 candidates) | done |
| 02 Ligand set curation | done — 271 ligands, 183 pairs |
| 03 Frozen train/test split | done — seed 20260902, test sealed |
| 04 Continuous-Δ validation set | done — 108 pairs, no gap |
| 05 Rung 1: 2D-QSAR baseline | done — 40 model runs |
| 06 Docking protocol validation | done — 0.75 ± 0.01 Å redocking |
| 07 Production docking (354 ligands) | **next** |
| 08–11 MD, descriptors, rungs 2–4 | not started |

The central hypothesis is **not yet tested**. What is established is that the
question is worth asking.

## Key results so far

**2D-QSAR cannot predict the cliffs, and it is not a matter of model choice or
data volume.** Across six algorithm families (Ridge, PLS, k-NN, SVM, Random
Forest, Gradient Boosting) and two descriptor sets, the best pairwise-Δp model
reaches R² = 0.200. Retraining on 5,016 compounds — 23× more data — raised
ordinary potency prediction to R² = 0.698, competitive with the published AChE
QSAR literature, while cliff recovery did not improve at all. The failure is
representational, not statistical.

**39 of 77 development cliff pairs are provably unpredictable.** Their binary
ECFP4 vectors are bit-identical, because Morgan fingerprints are blind to chain
length. The feature difference is exactly zero while true potency differs by
100× to 4,786×, so any fingerprint-based model must output Δp = 0.

**The docking protocol is validated and its limits are quantified.** Donepezil
redocks into 4EY7 at 0.75 ± 0.01 Å over five seeds. Cross-docking degrades for
small ligands with a different binding mode (huperzine A, 3.10 Å), and the cause
is identified: the Tyr337 side chain shifts 2.53 Å between structures, blocking
the true pose. Self-docking the same ligand into its own receptor gives 0.24 Å,
proving the limitation is receptor rigidity rather than the protocol.

## Layout

```
scripts/     pipeline, numbered by stage in results/INDEX.md
results/     outputs, prefixed 01_–06_ by the stage that produced them
docking/     receptors, ligand poses, validation runs
structures/  crystal structures 4EY7, 4EY5, 4EY6, 4M0F
figures/     publication figures (PNG 300 dpi + vector PDF)
data/        cached ChEMBL downloads
```

See [`results/INDEX.md`](results/INDEX.md) for the stage-by-stage file index.

## Reproducing

```bash
conda env create -f environment.yml
conda activate qsar-md
```

`environment.yml` is the portable spec; on a Linux/CUDA host conda-forge resolves
the GPU-enabled OpenMM build from the same file.
`environment-full-osx-arm64.yml` is the fully pinned local snapshot.

Two pinned choices matter for reproducibility: the ligand force field is
**openff-2.3.0 (Sage)**, and protonation at pH 7.4 uses this project's own
rule-based module (`scripts/protonation.py`, validated 16/16 on a panel including
six real AChE drugs) because dimorphite-dl 2.x returns chemically impossible
states.
