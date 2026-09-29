# Stage 18 — replicate pilot: is the descriptor the ligand, or the seed?

**Written before the trajectories were run.** The decision it informs — whether to spend $178
on the 162 deferred replicate trajectories — is the kind that invites a criterion chosen after
seeing the numbers, so the criterion is fixed below instead.

## The question

A replicate here is the same system simulated again from a different velocity seed. Same
ligand, same pose, same force field, same box; the system files are symlinked, so this is
provable rather than asserted. Each replicate is equilibrated independently from NVT onward,
which is what makes them independent samples rather than branches of one trajectory.

The objection to running them is reasonable: if 50 ns is long enough for the system to
converge, a different seed should give the same answer. The first three production runs show
that this is true for some descriptors and false for others, in the same 45 ns:

| | 5 ns | 45 ns | drift |
|---|---:|---:|---:|
| `occ_bridge` | 0.667 | 0.668 | 0.2% |
| `frac_lig_buried_mean` | 0.991 | 0.994 | 0.3% |
| `depth_mean` | 7.989 | 8.039 | 0.6% |
| `rmsf_gorge` | 0.626 | 0.714 | 12.4% |
| `n_bridge_water` | 2.007 | 1.693 | 18.6% |
| `gate_flips_per_ns` | 0.000 | 0.007 | — |

Convergence is a property of an observable, not of a system. A flat Cα RMSD says the system
stopped drifting; it says nothing about whether the time average of a slow quantity has
reached its ensemble value. `n_bridge_water` falls monotonically through every window and has
not levelled off at 45 ns, and a per-frame mean has no definitional excuse for that the way
an RMSF does.

**The duration ablation cannot answer the seed question, and this is the important part.**
It measures whether one trajectory's average has stopped moving. A quantity can be perfectly
flat in every individual run and still differ between seeds — which is exactly what happens
under basin trapping. `gate_flips_per_ns` averages 0.007/ns, about one Tyr337 rotamer jump
per three trajectories, so the gate essentially never moves in 45 ns. A run beginning in one
rotamer will report it at occupancy 1.000 with a flat RMSD and zero drift, and look converged.
Nothing inside a single trajectory distinguishes *"the answer is A"* from *"this run was stuck
in A"*. Only replicates do.

Production also uses `tcoupl = v-rescale` and `pcoupl = c-rescale`, both stochastic at every
step, so divergence continues throughout rather than being set once at t=0.

## Design

Six ligands, two extra replicates each: 12 runs of 50 ns, about 85 GPU-hours, **roughly $13**.

All six are cliff-pair ligands — the 81 that the full design gives three replicates — so these
12 runs are a literal subset of the deferred 162 and are not thrown away if the full set is
later funded. Their `_r1` trajectories come free from the instance already running, which is
why they were chosen: five of the six sit early in that slice.

| ligand | r1 source |
|---|---|
| CHEMBL118968 | done |
| CHEMBL1642973 | done |
| CHEMBL1643944 | running |
| CHEMBL303574 | slice position 10 |
| CHEMBL3642047 | slice position 17 |
| CHEMBL62664 | slice position 27 |

The last two are a real activity-cliff pair: **CHEMBL62664 / CHEMBL303574, ΔpChEMBL 2.28,
Tanimoto 1.00, constitutional.** With replicates on both sides it becomes possible to ask the
study's actual question rather than a proxy for it — is the within-pair descriptor difference
larger than the noise on that difference? — for one pair, at a cost of $2.

## The rule, fixed in advance

> A descriptor can carry the pair endpoint only if the spread between ligands is at least
> **twice** the spread between replicates of the same ligand.

Below that it is measuring the random number generator. The factor of 2 is a judgement; it is
a judgement made blind. For the pair endpoint specifically, a difference is two independent
estimates, so the noise on it is `sd_rep * sqrt(2)` and not `sd_rep`.

`scripts/replicate_pilot.py` implements this and reports per descriptor and per ablation
window, so it also shows whether a longer trajectory buys down the replicate noise or not.

Three outcomes, all of them decisive:

- **Most descriptors pass.** The 162 replicates are worth funding, and the pair endpoint
  stands on measured ground.
- **The slow descriptors fail and the fast ones pass.** Then `n_bridge_water`,
  `gate_flips_per_ns` and `bridge_runs_per_ns` are unmeasurable at 50 ns and should be dropped
  rather than reported — and dropping them is a finding, not a retreat.
- **Most fail.** Rung 4 at 50 ns cannot support a per-pair claim at all, and the honest
  conclusion is about the timescale rather than about the ligands.

## A reproducibility gap found while preparing this

`md/mdp/prod.mdp` sets no `ld-seed`, so grompp draws a random one per invocation and writes
it only into the `.tpr`. The manifest records `gen_seed` — the NVT velocity seed — and nothing
about the production thermostat's random stream. Confirmed by dumping the tprs: three
replicates of one ligand carry ld-seeds −8801800, −1208092685 and 1252786170.

This is good for the pilot, which wants independence, and bad for reproducibility, which the
paper will claim. The seeds are recoverable only from the tprs, which exist now and will not
exist forever, so they were harvested into `results/14d_prod_ld_seeds.csv` while they could be.
