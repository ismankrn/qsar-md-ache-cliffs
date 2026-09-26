# Stage 12 — OpenFF/OpenMM → GROMACS export validation

Rung 4 is parameterised with OpenFF Sage 2.3.0 through OpenMM, matching rung 3, but
production runs GROMACS on BRIN. A faulty export would not fail loudly — it would yield
plausible-looking, wrong trajectories. This stage proves the export is faithful before
379 systems are built on it.

**System:** 4EY7 + CHEMBL104387 (seed 2 pose), TIP3P, 0.15 M NaCl, 91,381 atoms,
98.54 Å cubic box. Protein `amber14-all.xml`, ligand OpenFF 2.3.0 via
`SMIRNOFFTemplateGenerator`. Single-point, PME, 1.0 nm cutoffs, no constraints.

## Result — export is faithful

| term | relative difference |
|---|---:|
| Torsions | **4.2e-08** |
| Angles | 1.7e-05 |
| Bonds | 2.6e-05 |
| Nonbonded | 4.0e-05 |
| Total potential | 7.0e-05 |

Per-atom forces, which are what determine dynamics (a constant energy offset moves nothing):

```
rms deviation 0.2144 kJ/mol/nm  against rms magnitude 1468.48   ->  1.5e-04
component correlation  0.999999999294
```

Torsions matching to 4e-8 is the strongest evidence: torsions are the most
parameter-dense term (hundreds of types, each with phase, multiplicity, periodicity)
and the most likely to be corrupted by a bad export.

The residual is fully explained. Ligand forces deviate by 1.75 rms versus 0.16 for
protein and 0.21 for solvent, because ParmEd writes `b0` to 5 decimals in nm while
OpenFF Sage values are arbitrary reals. A 1.8e-5 nm rounding on a bond with
k ≈ 265,000 gives ~9 kJ/mol/nm — the order of the observed 2.95 maximum. Physically
this shifts equilibrium bond lengths by 1e-4 Å, three orders below thermal
fluctuation. Irrelevant to every descriptor.

## Confounds that had to be removed first

Each of these differs by default between the two engines and masquerades as an export
error. Four wrong hypotheses were tested and discarded before the real cause was found
by reading the generated topology rather than guessing.

| confound | default OpenMM | default GROMACS | resolution |
|---|---|---|---|
| Dispersion correction | on | off | off in both |
| vdW cutoff treatment | hard cutoff | potential-shift | raw in both (`vdw-modifier = none`) |
| PME grid and alpha | chosen independently | chosen independently | GROMACS forced to OpenMM's, via `ewald-rtol = erfc(alpha*rc)` |
| Coordinate precision | exact | `.gro` = 0.001 nm | ship a binary `.trr` |
| **Water rigidity** | `rigidWater=False` → flexible | takes `[ settles ]` → rigid | `define = -DFLEXIBLE` |

**The `[ settles ]` trap is the one to remember.** ParmEd writes water twice: flexible
bonds and angles under `#ifdef FLEXIBLE`, and a rigid `[ settles ]` block under `#else`.
GROMACS takes settles by default, so the comparison was between rigid and flexible water
— two different physical systems, not an export error. This alone accounted for a force
deviation of 22.90 kJ/mol/nm; defining `-DFLEXIBLE` collapsed it to 0.2144, a factor of 107.

Two hypotheses that were tested and **refuted**, recorded so they are not re-tried:
- *`.gro` rounding explains the residual.* It does not. A `.trr` round trip shifts
  coordinates by at most 1.1e-5 Å and the bond energy by 0.055 kJ/mol, not the 21.7 seen.
- *PME discretisation explains the residual.* It does not. Refining OpenMM's own grid
  from 79³ to 158³ moves its forces by rms 0.0428 kJ/mol/nm — three orders too small.

## Settings for production — the opposite of the validation

The validation deliberately used flexible water and no constraints to isolate parameters.
Production must **not**:

- leave `FLEXIBLE` undefined, so water uses `[ settles ]`
- `constraints = h-bonds`

Rigid water plus h-bond constraints is what permits the 2 fs timestep. Flexible water
would force ~0.5 fs and quadruple the cost for nothing.

Reproduce with `scripts/validate_gmx_export.py`; artefacts in `md/export_check/`.
Note that `addSolvent` places water and ions stochastically, so the OpenMM reference is
valid only for the system exported in the same run — verify by checksum before comparing.
