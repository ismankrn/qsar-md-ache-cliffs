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

---

## Result (2026-10-02): 12 of 12 replicate runs complete

Six ligands, two extra replicates each, against 30 ligands for the between-ligand spread.
The rule fixed in advance — between-ligand spread at least twice the between-replicate
spread — gives, at the full 45 ns window:

**15 usable, 15 seed-dominated, 1 sparse, 1 constant.**

Half the descriptor set measures the random number generator more than it measures the
ligand. Without replicates all 30 would have gone into the model and 15 of them would
have been noise wearing the costume of signal. That alone is what the $13 bought.

### What survives, and what does not

| usable | ratio | | seed-dominated | ratio |
|---|---:|---|---|---:|
| `occ_86` | 465 | | `occ_121` | 1.91 |
| `occ_cas` | 201 | | `rmsd_max` | 1.72 |
| `occ_bridge` | 39.3 | | `rmsf_lig` | 1.46 |
| `occ_pas`, `occ_286` | 37.4 | | `bridge_longest_ns` | 1.30 |
| `occ_124` | 11.9 | | `rmsf_gorge` | 1.22 |
| `occ_72` | 5.9 | | `rmsd_mean` | 0.99 |
| `occ_337` | 5.0 | | `ca_rmsd_mean` | 0.83 |
| `lig_dihedral_sd` | 2.7 | | `bridge_runs_per_ns` | 0.71 |
| `frac_lig_buried_mean` | 2.6 | | `gate_flips_per_ns` | 0.36 |
| `depth_mean` | 2.6 | | | |
| `n_bridge_water` | 2.2 | | | |

**The CAS–PAS bridging family is the strongest in the set.** The hypothesis the study was
built on is carried by the descriptors that reproduce best, not by the ones that happen to
be noisy. The casualties are the mobility metrics — ligand RMSD and RMSF distinguish seeds
about as well as they distinguish ligands — and the rate descriptors, which was already
expected from `gate_flips_per_ns` averaging about one event per three trajectories.

`ca_rmsd_mean` is the positive control and it behaves: it is a property of the receptor,
it should not discriminate between ligands, and at 0.83 it lands squarely in the
seed-dominated half. A method that flags a known negative as negative is worth more than
one that only produces the answers it was pointed at.

### How much of this survives at 20 ns

The count of usable descriptors is close to flat with trajectory length:

| window | usable | seed-dominated |
|---:|---:|---:|
| 5 ns | 12 | 17 |
| 10 ns | 13 | 16 |
| 20 ns | 12 | 18 |
| 30 ns | 13 | 17 |
| 45 ns | **15** | 15 |

Twenty nanoseconds of analysed trajectory keeps 12 of the 15. The three bought by the
full run are `frac_lig_buried_mean` (0.89 at 5 ns, 2.64 at 45), `occ_124` (2.15 to 11.9)
and `gate_chi1_sd`. Everything else is already there by 10 ns, and some are *better*
short: `occ_337` falls from 22.4 to 5.0 as the window lengthens, and `occ_bridge` from
65.7 to 39.3.

More decisive than the counts: **shortening does not make the seed noise worse.**

| `sd_rep` | 5 ns | 10 ns | 20 ns | 30 ns | 45 ns |
|---|---:|---:|---:|---:|---:|
| `occ_bridge` | 0.006 | 0.006 | 0.006 | 0.006 | 0.010 |
| `n_bridge_water` | 0.662 | 0.557 | 0.547 | 0.572 | 0.584 |
| `depth_mean` | 1.041 | 1.030 | 0.948 | 0.923 | 0.920 |

Reproducibility of the descriptors that matter is flat or better at the shorter windows.
The case for 50 ns rests on three descriptors, not on the set as a whole.

### Limits of this reading

Six ligands set `sd_rep`, and two of them have two replicates rather than three because
their `_r1` is still queued on vast1. Ratios for the near-saturated occupancies
(`occ_86`, `occ_cas`, occupancy ~1.0) swing between 387 and 2376 across windows — the
denominator is near zero and the precision is illusory; read them as "reproducible", not
as a number. All six are cliff-pair ligands, which is the right population for a pair
endpoint but not a random sample of the 217.

### Robustness: restricting to the four ligands with all three replicates

Dropping to `--min-reps 3` removes the two ligands whose `_r1` is still queued, leaving
four to set `sd_rep`. The verdict count moves from 15 usable to 13, and the useful
question is which descriptors move with it.

**Stable** — same order of magnitude, same verdict, on either subset: the whole CAS–PAS
family (`occ_cas` 398 vs 201, `occ_86` 380 vs 465, `occ_bridge` 33 vs 39, `occ_pas` 31 vs
37), `frac_lig_buried_mean` (2.61 vs 2.64), `lig_dihedral_sd` (2.60 vs 2.72), `depth_mean`
(2.28 vs 2.59), and every descriptor already far below the line (`bridge_runs_per_ns`,
`occ_338`, `ca_rmsd_mean`).

**Unstable** — flips across the threshold depending on the subset:

| | 4 ligands | 6 ligands |
|---|---:|---:|
| `n_bridge_water` | 1.85 seed-dominated | 2.16 usable |
| `depth_sd` | 1.82 seed-dominated | 2.09 usable |
| `rmsf_lig` | 1.99 seed-dominated | 1.46 seed-dominated |

So the threshold of 2.0 is doing real work only for descriptors that sit near it, and for
those three the pilot does not settle the question — six ligands is not enough to place a
ratio to better than about ±0.3. `n_bridge_water` matters most of the three, being the
one genuinely novel physical quantity in the set, and its status should be treated as
open rather than decided.

What the pilot does settle is the two ends: the bridging occupancies are reproducible by
a wide margin, and the mobility and rate descriptors are not, on either subset.
