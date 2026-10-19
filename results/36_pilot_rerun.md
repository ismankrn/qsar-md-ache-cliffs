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

## Bootstrap intervals, recomputed

300 replicates over units, the algorithm re-selected inside every replicate.

| framing | gap, first run | gap, re-run | interval, re-run | excludes zero |
|---|---:|---:|---|---|
| A, $R^2$ | −0.384 | −0.375 | [−0.596, −0.184] | **yes** |
| B, $R^2$ | −0.111 | −0.093 | [−0.319, +0.127] | no |
| C, ROC AUC | −0.002 | +0.003 | [−0.102, +0.103] | no |

The structure is unchanged. Framing A is the only resolved comparison in the study and
it is resolved against rung 4, while B and C remain unresolved. The bootstrap levels
carry the known duplicate-resampling bias and are not performance estimates; the gap is
the quantity to read.

## Duration ablation, re-run

`fig_duration.py` reads the filter live, so re-running it was the whole fix. The fixed
series holds the 20 ns selection at every window and is the series that isolates length.

| window | A, first run | A, re-run | B, first run | B, re-run | C, first run | C, re-run |
|---:|---:|---:|---:|---:|---:|---:|
| 5 ns | 0.039 | **0.030** | 0.116 | **0.149** | 0.792 | **0.747** |
| 10 ns | 0.027 | **0.044** | 0.198 | **0.196** | 0.805 | **0.788** |
| 20 ns | 0.056 | **0.043** | 0.124 | **0.152** | 0.793 | **0.796** |

The 20 ns column now matches the ladder exactly, which is the check that matters: the
same descriptor set scored by the same code has to give the same number in both places,
and it does (A 0.043, B 0.152, C 0.796).

The reading is unchanged and is in places firmer. Framing A moved over a range of 0.029
across the three windows and now moves over 0.014. Framing B still peaks at 10 ns.
Framing C now rises monotonically with length, by 0.049 from 5 to 20 ns, where before it
moved without direction; that is well inside the bootstrap interval of roughly 0.10 on
this framing, so it does not become a dependence on length.

The per-window filter is also less uneven than before. It passes 13, 13 and 11
descriptors at 5, 10 and 20 ns, against 14, 13 and 11 before, but the membership at 5
and 10 ns still differs by one, so the adaptive series continues to confound length with
descriptor choice and the fixed series remains the one to read.

The 45-ligand subset that reaches 30 ns does not rise past 20 ns, as before.
