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

## Stages not yet run

| Stage | Description | Gate |
|---|---|---|
| 07 | Production docking, 354 ligands | free, ~3–6 h locally |
| 08 | System building and parameterisation | free |
| 09 | Pilot MD, 5 ligands | free, validates the MD pipeline |
| 10 | Production MD, 354 ligands | ~$482 on rented GPUs |
| 11 | Rungs 2–4: static 3D, MM/GBSA, dynamic descriptors | free |
