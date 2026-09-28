# results/ — index by pipeline stage

Files are prefixed with the stage that produced them. Stage numbers follow the
logical pipeline order, not the order in which the work happened.

| Stage | Produced by | Files |
|---|---|---|
| **01** ChEMBL screening & target selection | `chembl_screen.py` | `01_feasibility.csv`, `01_compounds_*.csv`, `01_cliffs_*.csv` |
| **02** Ligand set curation | `build_ligand_set.py` | `02_ligand_set_CHEMBL220.csv`, `02_pairs_CHEMBL220.csv`, `02_cliff_similarity_audit.csv` |
| **03** Frozen train/test split | `make_splits.py` | `03_splits_ligands_CHEMBL220.csv`, `03_splits_pairs_CHEMBL220.csv` |
| **04** Continuous-Δ validation set | `build_validation_set.py` | `04_validation_ligands_CHEMBL220.csv`, `04_validation_pairs_CHEMBL220.csv` |
| **05** Rung 1 — 2D-QSAR baseline | `rung1_2d_qsar.py`, `rung1_benchmark.py`, `rung1_robustness.py` | `05_rung1_benchmark.csv`, `05_rung1_robustness.csv`, `05_baseline_2d_predictions.csv` |
| **06** Docking protocol validation | `redock_replicas.py`, `box_experiment.py` | `06_redocking_validation.csv`, `06_box_experiment.csv` |
| **07** Production docking | `run_docking.py`, `analyse_docking.py` | `07_docking_scores.csv`, `07_docking_per_ligand.csv`, `07_docking_pair_analysis.csv` |
| **08** Rung 2 — pose geometry | `rung2_descriptors.py`, `rung2_bridging.py`, `rung2_model.py` | `08_rung2_descriptors.csv`, `08_rung2_bridging.csv`, `08_rung2_bridging_pairs.csv`, `08_rung2_benchmark.csv` |
| **09** Rung 3 — MM/GBSA | `rung3_mmgbsa.py`, `rung3_analyse.py` | `09_rung3_mmgbsa.csv`, `09_rung3_pair_analysis.csv` |
| **10** Size-adjusted tests | `size_adjusted_tests.py` | `10_size_adjusted_tests.csv` |
| **11** Ensemble docking (8 receptors) | `select_ensemble.py`, `prep_ensemble.py`, `run_ensemble_docking.py`, `analyse_ensemble.py` | `11_ensemble_docking.csv`, `11_ensemble_per_ligand.csv`, `11_ensemble_pair_analysis.csv` |
| **12** OpenFF→GROMACS export validation | `validate_gmx_export.py` | `12_gmx_export_validation.md` |
| **12b** Receptor chain breaks capped for MD | `prep_receptor_md.py` | `12b_receptor_chain_breaks.md` |
| **13** MD system build (217 systems) | `build_md_systems.py` | `13_md_ligands.csv`, `13_md_system_build.md`, `13b_pose_clearance.csv` |
| **14** MD run layout (379 trajectories) | `make_md_runs.py` | `14_md_runs.csv` |
| **15** Rung 2 on ensemble poses + bootstrap | `rung2_ensemble.py`, `rung2_ensemble_control.py`, `rung2_ensemble_model.py`, `rung2_ensemble_bootstrap.py` | `15_rung2_ensemble.md`, `15_rung2_ensemble_*.csv`, `15b_rung2_ensemble_median_fixed4ey7.csv` |
| **16** Rung 4 dynamic descriptors (running) | `rung4_descriptors.py` | `16_rung4_per_run/*.json` |
| **16b** PBC + superposition defect in stage 16, and its repair | `rung4_descriptors.py` | `16b_rung4_pbc_fix.md` |
| **16c** Rung 4 QC plots and threshold flags | `rung4_qc_plots.py` | `16c_qc_summary.png`, `16c_qc_flags.csv`, `16c_qc/<run>.png` (untracked) |
| **16d** Duration ablation: R² against simulated time | `rung4_descriptors.py` | `16d_duration_ablation.md`, `windows` key in each run JSON |
| **17** Power analysis and noise ceiling | `power_and_ceiling.py`, `power_proper.py` | `17_power_analysis.md` |

## Notes

- **Only CHEMBL220 (hAChE) is carried forward.** The other five targets in stage 01
  were screening candidates; their files are kept as the audit trail for target
  selection but are not used downstream.
- **`05_baseline_2d_predictions.csv` is superseded.** It came from `baseline_2d.py`,
  which computed its own folds instead of loading the frozen split and filtered
  zero-variance descriptors across the whole dataset. Use `05_rung1_benchmark.csv`
  and `05_rung1_robustness.csv` instead.
- **Stage 03 is frozen.** Seed 20260902. Every descriptor rung must be scored on
  these folds, or the comparison between rungs is meaningless.
- **The test split is sealed** and has not been scored by any model.

## Where the ladder stands

| Rung | Descriptors | Framing B, best R² |
|---|---|---:|
| 1 | 2D fingerprints + physicochemical | **+0.200** |
| — | raw docking score | cliff/control separation −0.03 kcal/mol |
| 2 | static pose geometry | +0.167 |
| 3 | single-snapshot MM/GBSA | −0.055 |

After adjusting for molecular size, **0 of 19 descriptors** separate cliff
pairs from control pairs. No static structure-based descriptor explains these
activity cliffs, and the CAS–PAS bridging hypothesis is unsupported by both
geometry and energetics.

`docking/poses/` is not versioned: 1,416 files, 16 MB, regenerable by
`run_docking.py` with its fixed seeds, or moved with `scripts/sync.sh`.

## The size confound

Cliff partners differ by 4.4 heavy atoms on average, controls by 1.8, so any
size-correlated descriptor separates the groups for the wrong reason. Four
descriptors that looked significant — `n_contacts_4A`, `dSASA`,
`dG_nonpolar`, `axis_span_frac` — are size artefacts.

Matching harder does not fix it. Across all 1,911 similar AChE pairs, potency
difference grows with size difference: pairs differing by 0–2 heavy atoms are
iso-potent 48% of the time, pairs differing by 6–9 only 18%. Size-matched
iso-potent controls are intrinsically rare. `rebuild_controls.py` records the
attempt — it cost 26 of 86 control pairs and still left a 2.9 vs 4.4 gap, so
it is not used. The confound is handled statistically instead, in
`size_adjusted_tests.py`.

Note on that script: trust the **stratified** column, not the residual one.
Regressing out a covariate that itself differs between groups can manufacture
differences in the residuals.

## Stages not yet run

| Stage | Description | Cost |
|---|---|---|
| 11 | Ensemble docking — tests whether wrong poses explain the failure | free, ~6 days |
| 12 | MD, 16 ligands — tests whether dynamics explains it | ~$43 |
