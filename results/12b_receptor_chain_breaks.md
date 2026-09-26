# Stage 12b — fabricated peptide bonds in the receptor, and how they were removed

## The defect

4EY7 is missing residues 259-264 and 495-497. `prep_receptor.py` calls
`findMissingResidues()` and has logic to fill internal gaps, but it never fires: the
intermediate file is written by MDAnalysis, which emits no SEQRES records, and PDBFixer
identifies missing residues by comparing coordinates against SEQRES. With no SEQRES it
reports nothing missing, prints nothing, and the gaps pass through silently.

`app.ForceField` then treats each chain as continuous and fabricates peptide bonds across
the gaps:

| bond | distance | energy |
|---|---:|---:|
| C258-N265 | 7.108 Å | 136,665 kJ/mol |
| C494-N498 | 4.050 Å | 30,220 kJ/mol |
| | **sum** | **166,885** |
| | measured total bond energy | 185,963 |

**90% of the system's entire bond energy was two bonds that do not exist.**

It surfaced from a single `grompp` line in stage 12 — a warning about excluded atoms
0.908 Å apart — which was nearly dismissed as noise.

## Earlier stages are unaffected

Docking and rung 2 use a rigid receptor. Rung 3 computes
`E_complex - E_protein - E_ligand` with the protein frozen, so the fabricated bonds appear
identically in two terms and cancel exactly. **Nothing needs to be redone.** Only MD is
sensitive: minimisation would drag the segments together with ~10^5 kJ/mol/nm of force and
deform the fold.

## Loop modelling was tried and abandoned

| attempt | outcome |
|---|---|
| Minimise, only modelled residues free | pinned against frozen flanking atoms, 6.75e6 kJ/mol |
| Free a 3-residue flank as well | 0 clashes, but **4 non-proline cis peptides** |
| Add omega restraints to trans | all trans, but 5 new clashes, bond energy 100,831 |
| Anneal at 500 K | diverged to NaN |

Non-proline cis peptides occur in ~0.03% of real residues; four in one six-residue loop is
an artefact. PDBFixer builds loops crudely and expects minimisation to rescue them, and
here it could not. Expected bond energy for this protein is ~1,000 kJ/mol; the best attempt
reached 100,831.

## What was done instead: cap the breaks

What has to go is the fabricated bond, not the absence of the loop. The chain is split
into three segments at the breaks, each given its own chain id, and every segment terminus
is capped — ACE on the N-side, NME on the C-side. Caps are placed by internal coordinates,
so the procedure is deterministic: no minimisation of the protein, no random seed.

Neutral caps rather than charged termini, which would have put six artificial formal
charges on the protein surface.

    segments: 4-258, 265-494, 498-541      3 chains, 3 ACE + 3 NME, 8,173 atoms

| check | before | after |
|---|---:|---:|
| Bond energy | 185,963 | **1,145** |
| Nonbonded | — | **−16,308** (negative: properly packed) |
| Total potential | — | 17,323 |
| Max force | 3.3e8 | 2.1e5 |
| Gorge displacement | — | **0.0000 Å** |
| Short contacts < 1.6 Å | 21 | 18 |

**The gorge is bit-identical**, so all 354 docked poses from stage 07 remain valid. The
remaining sub-1.6 Å contacts are geminal hydrogen pairs from hydrogen placement; there were
more of them in the original receptor, and only one involves a cap. They relax in the first
minimisation.

## Three traps worth remembering

1. **Writing an intermediate PDB through MDAnalysis drops SEQRES**, which silently disables
   PDBFixer's missing-residue detection. Feed PDBFixer the original file.
2. **`Topology.createStandardBonds()` creates no bonds to hydrogen**, and misses NME's
   N-CH3. Caps built by hand need their bonds built by hand.
3. **PDBFixer adds OXT to what it takes to be the chain terminus.** Capping that residue
   with NME then gives it both a terminal carboxylate and an amide link; the template
   matcher reports it as "matches CSER, but ... 1 C atom too many". Drop OXT before capping.

Produced by `scripts/prep_receptor_md.py` -> `docking/4ey7_rec_md.pdb`.
Use this file for MD. `4ey7_rec_clean.pdb` remains the docking receptor and must not be
regenerated, or the poses would no longer correspond to it.
