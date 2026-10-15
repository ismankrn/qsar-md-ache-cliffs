# Stage 36 — the replicate pilot re-run on six ligands

The filter published in stage 18 was computed while two of the six replicate sets were
still running, and while only 30 of the 217 ligands had finished. It therefore rested
on four ligands for the replicate spread and 30 for the between-ligand spread. Both
sets are now complete, so the pilot was re-run unchanged.

The rule was not touched. A descriptor passes only when the spread between ligands is
at least twice the spread between replicates of the same ligand, which is the criterion
fixed before any replicate trajectory existed.

## What changed in the filter

| | first run | re-run |
|---|---:|---:|
| ligands with replicates | 4 | **6** |
| ligands for the between-ligand spread | 30 | **217** |
| seed-dominated, of 160 combinations | 96 | **100** |
| descriptors passing at 20 ns | 11 | **11** |

The count is unchanged and the membership is not. Nine of the eleven survivors are
shared. `frac_lig_buried_mean` and `n_bridge_water` drop out; `bridge_runs_per_ns` and
`occ_124` come in. Nine of the eleven are now contact occupancies, against eight
before, so the finding that only persistence survives in quantity is strengthened
rather than disturbed.

## What changed in the ladder

Only rung 4 filtered moves, and only on two framings.

| framing | first run | re-run | rung 1 | gap |
|---|---:|---:|---:|---:|
| A, $R^2$ | 0.056 | 0.043 | 0.386 | −0.317 |
| B, $R^2$ | 0.124 | **0.152** | 0.200 | −0.048 |
| C, ROC AUC | 0.793 | **0.796** | 0.805 | −0.010 |

Rung 4 remains below rung 1 on all three framings, so nothing in the conclusion turns
on this. The framing-B gap narrows from −0.077 to −0.048 and stays negative.

The superseded output is kept as `18_replicate_pilot_4lig_superseded.csv`.
