# Stage 29 — literature check on the three claims that needed one

`INDEX.md` has carried an OPEN AT DRAFTING note since the descriptors were designed: the
rung-4 set had no citations, and two novelty claims were unverified. This is what eight
targeted searches found on 2026-10-13. Where nothing was found, that is recorded as *not
found*, which is not the same as *does not exist*.

## 1. Has anyone used MD-derived descriptors on activity cliffs?

**No study found that tests MD-derived descriptors for activity-cliff prediction.** The
surrounding literature is well populated and the gap is specific, which is the useful form
of this answer.

**The closest precedent is structure-based, and its structure favours us.** Husby, Bottegoni
et al., *Structure-Based Predictions of Activity Cliffs*, J. Chem. Inf. Model. 2015,
55(5):1062--1076, doi:10.1021/ci500742b. Read carefully before citing the abstract, which
says activity cliffs "can be accurately predicted by advanced structure-based methods":

| arm | accuracy |
|---|---|
| template docking, using the X-ray pose of the **stronger** ligand | median 90% |
| ALiBERO pocket set, **selected using the cliff labels** | 71--100% |
| **cross-docking, the prospective setting** | **median 53%** |

The authors call the ALiBERO arm "a purely retrospective exercise" themselves. The benchmark
is 146 cliff pairs from 158 X-ray complexes across nine targets, so **both partners have a
crystal structure** -- a setting this study does not have and most prospective work does not
either. The task is pairwise ranking, not regression on a signed difference. And there is
**no ligand-based or 2D comparison anywhere in the paper**.

So the honest positioning is not "nobody has done structure-based cliff prediction" but:
the strongest reported result uses the answer's crystal pose, the realistic arm sits near
chance for a two-way ranking, and no 2D baseline was run. Our contribution is the
prospective setting, the harder endpoint, and the matched baseline.

Also relevant: Dablander et al., *Exploring QSAR models for activity-cliff prediction*,
J. Cheminform. 2023, doi:10.1186/s13321-023-00708-w -- nine QSAR models on three targets,
finding they frequently fail on cliffs, and that ECFP remains the strongest representation.

A title exists, *Pharmacophore Modeling of Targets Infested with Activity Cliffs via
Molecular Dynamics Simulation Coupled with QSAR* (KDR case study), which could not be
retrieved or confirmed. **Find it before submission**; it is the one item that could change
this answer.

## 2. Has subset-baseline inflation been named?

**Not found**, in two searches with different phrasings. Adjacent ideas exist -- coverage
bias in public datasets (Nat. Commun. 2024), sampling and publication bias in QSAR reviews --
but none describes the specific mechanism: an expensive descriptor set arriving incrementally
in an order set by computational cost, compared against a baseline's published full-set
figure.

Two searches is thin evidence. **Word the claim as "we are not aware of a named treatment"
rather than "the first"**, and keep the demonstration itself as the contribution, which does
not depend on priority.

## 3. Citations per descriptor group

| group | status |
|---|---|
| `occ_*` contact occupancies | **SIFt**: Deng Z, Chuaqui C, Singh J. *J. Med. Chem.* 2004, 47(2):337--344. The foundational residue-wise interaction fingerprint. |
| the same, over a trajectory | **Not novel.** ProLIF (Bouysset & Fiorucci, *J. Cheminform.* 2021, doi:10.1186/s13321-021-00548-6) encodes interaction fingerprints from MD trajectories, and a 2024 paper (doi:10.1186/s13321-024-00822-3) aggregates them into occupancies with a 30% presence threshold. **Cite both and drop any claim that trajectory-averaged occupancy is new.** |
| `depth_mean`, `depth_sd` (CAS--PAS axis projection) | no precedent found; our own construction, as recorded |
| `rmsf_*`, `n_bridge_water`, `gate_chi1_sd`, `gate_flips_per_ns` | **still unverified** -- these are standard MD observables and almost certainly have canonical citations, but they were not searched for here |

## What this changes in the manuscript

1. **Husby et al. must be cited and engaged with**, not just listed. The abstract's claim and
   the paper's cross-docking number point in opposite directions, and saying so is a
   contribution in itself.
2. **The trajectory-occupancy descriptors are not novel machinery.** The novelty is the
   application to cliffs, the matched baseline and the stated detection limit. Reword
   anywhere that implies otherwise.
3. **MoleculeACE carries a correction.** van Tilborg, Alenicheva & Grisoni, *J. Chem. Inf.
   Model.* 2022, 62(23):5938--5951, doi:10.1021/acs.jcim.2c01073, with a later correction
   for a bug in train/test splitting caused by mislabelled cliff pairs. If its numbers are
   quoted, quote the corrected ones and say which.
4. Remaining work: the KDR paper, and citations for the fluctuation and water descriptors.
