# Stage 33 — does the 2D baseline change if it sees the protonated ligand?

Rung 1 is computed on `smiles_std`, the standardised neutral structure. Docking,
MM/GBSA and the simulations all use `smiles_md`, the dominant protonation state at
pH 7.4. The ladder therefore does not describe one chemical species from bottom to top.

This is not cosmetic. 232 of the 271 ligands differ between the two columns, and ECFP4
scores the neutral and protonated form of the *same* molecule at a median Tanimoto of
0.706 — more different than many of the cliff pairs are from each other.

## Result

Rung 1 rescored on `smiles_md`, everything else held identical: the frozen folds, both
feature sets, the same six and five algorithm families, the same pipelines, the same
three framings. Each cell is the best of the two feature sets, as rung 1 is scored.

| framing | neutral | protonated | change |
|---|---:|---:|---:|
| A (per-ligand pChEMBL, R²)   | 0.386 | **0.403** | +0.018 |
| B (signed pair Δ, R²)        | 0.200 | **0.200** | +0.000 |
| C (cliff/control, ROC AUC)   | 0.805 | **0.816** | +0.010 |

The neutral column reproduces the published ladder exactly, which is the check that
this script replicates rung 1 faithfully.

**The protonated baseline is never worse.** The published rung 1 is, if anything,
slightly understated, so every gap between rung 1 and a structure-based rung is
slightly *under*stated too. The conclusion that no structure-based rung resolvably
beats 2D is unaffected and marginally strengthened.

## Why framing B does not move at all

Identical to six decimal places, which is too clean to be luck. Protonation flips a
median of 22 ECFP4 bits per ligand (max 40). But framing B trains on the *difference*
between the two members of a pair, and the members of a near-identical pair protonate
the same way, so the change cancels:

| | |
|---|---:|
| dev pairs whose ECFP4 difference vector is unchanged | 130 / 145 (89.7%) |
| differing bits in the remaining 15 pairs | median 4, max 8, of 2048 |

That residue is small enough that the fitted model is unchanged to six decimals.

This is worth stating in the paper. On the endpoint that defines an activity cliff, the
protonation state is invisible to a 2D fingerprint *by construction*, because what the
model sees is a difference between two molecules that protonate alike.

## What the manuscript must say

1. Rung 1 uses the neutral structures, the structure-based rungs the protonated ones.
2. This was tested, and the baseline does not fall.

Produced by `scripts/rung1_protonation.py` -> `results/33_rung1_protonation.csv`.
