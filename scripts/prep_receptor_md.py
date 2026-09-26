"""Stage 12b — receptor preparation for MD: chain breaks capped, no loops modelled.

The problem this solves. 4EY7 is missing residues 259-264 and 495-497. `prep_receptor.py`
writes its intermediate file with MDAnalysis, which emits no SEQRES, so PDBFixer never
learns those residues are absent and the gaps pass through silently. `app.ForceField`
then treats each chain as continuous and fabricates peptide bonds across them:

    C258-N265   7.108 A  ->  136,665 kJ/mol
    C494-N498   4.050 A  ->   30,220 kJ/mol

That is 90% of the system's entire bond energy. Docking and rung 2 used a rigid receptor,
and rung 3 computes E_complex - E_protein - E_ligand with the protein frozen, so the fake
bonds cancel exactly -- none of those stages are affected. MD is not so forgiving:
minimisation would drag the segments together and deform the fold.

Why capping rather than loop modelling. PDBFixer builds loops crudely and expects
minimisation to rescue them; here it could not. Three attempts produced, in turn, a loop
pinned against frozen flanking atoms (6.75e6 kJ/mol), four non-proline cis peptides,
and -- with omega restrained to trans -- five fresh clashes and 100,831 kJ/mol of bond
energy against an expected ~19,600. Annealing at 500 K then diverged to NaN. What
actually has to go is the fabricated bond, not the absence of the loop. Capping removes
it without inventing a conformation that cannot be defended, at the cost of two surface
holes 21-49 A from the gorge, which is not where this study looks.

The chain is split into three segments at the breaks, each given its own chain id, and
every segment terminus is capped with ACE (N-side) or NME (C-side). Caps are placed by
internal coordinates, so the result is deterministic -- no minimisation, no random seed.
Neutral caps are used throughout rather than charged termini, which would introduce six
artificial formal charges into the protein surface.
"""
from pathlib import Path
import numpy as np
from pdbfixer import PDBFixer
from openmm import unit
from openmm.app import (PDBFile, ForceField, Modeller, Topology, Element,
                        Simulation, CutoffNonPeriodic)

ROOT = Path(__file__).resolve().parent.parent
STRUCT, DOCK = ROOT/"structures", ROOT/"docking"
PDB, CHAIN, PH = "4ey7", "A", 7.4
DROP_RESID = 542          # disordered C-terminus; see prep_receptor.py

# ---------------------------------------------------------------- stage the input
keep = []
for line in (STRUCT/f"{PDB}.pdb").read_text().splitlines():
    rec = line[:6]
    if rec == "SEQRES" and line[11] == CHAIN:
        keep.append(line)
    elif rec == "ATOM  " and len(line) > 21 and line[21] == CHAIN:
        try:
            if int(line[22:26]) == DROP_RESID:
                continue
        except ValueError:
            pass
        keep.append(line)
    elif rec in ("HEADER", "CRYST1", "SSBOND"):
        keep.append(line)
keep.append("END")
staged = DOCK/f"{PDB}_rec_seqres.pdb"
staged.write_text("\n".join(keep) + "\n")

# ---------------------------------------------------------------- repair, but no loops
fixer = PDBFixer(filename=str(staged))
fixer.findMissingResidues()
print(f"celah yang dilaporkan PDBFixer: {len(fixer.missingResidues)}")
fixer.missingResidues = {}                 # model nothing; the breaks get capped instead
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()                    # sidechain atoms only, now that loops are off

pos = np.array(fixer.positions.value_in_unit(unit.angstrom))
residues = list(fixer.topology.residues())
resids = [int(r.id) for r in residues]

# ---------------------------------------------------------------- find the breaks
breaks = [i for i in range(len(residues)-1) if resids[i+1] != resids[i] + 1]
print(f"putus rantai: {len(breaks)}")
for i in breaks:
    print(f"  antara {residues[i].name}{resids[i]} dan {residues[i+1].name}{resids[i+1]}"
          f"  ({resids[i+1]-resids[i]-1} residu hilang)")
segments = []
start = 0
for i in breaks:
    segments.append((start, i)); start = i + 1
segments.append((start, len(residues)-1))
print(f"segmen: " + ", ".join(f"{resids[a]}-{resids[b]}" for a, b in segments))

# ---------------------------------------------------------------- cap geometry
def place(a, b, c, bond, angle_deg, dih_deg):
    """Position an atom bonded to c, with angle(b,c,new) and dihedral(a,b,c,new)."""
    ang, dih = np.radians(angle_deg), np.radians(dih_deg)
    bc = c - b; bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc); n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-bond*np.cos(ang), bond*np.sin(ang)*np.cos(dih), bond*np.sin(ang)*np.sin(dih)])
    return c + d[0]*bc + d[1]*m + d[2]*n

def atom_xyz(res, name):
    for a in res.atoms():
        if a.name == name:
            return pos[a.index]
    raise KeyError(f"{name} tidak ada di {res.name}{res.id}")

def methyl_h(a, b, c):
    """Three tetrahedral hydrogens on methyl carbon c, which is bonded to b."""
    return [place(a, b, c, 1.090, 109.5, d) for d in (60.0, 180.0, 300.0)]

def ace_for(res):
    """ACE capping the N-terminus of res: CH3-C(=O)- bonded to its N.

    Modeller.addHydrogens leaves the cap methyls without hydrogens and template
    matching then rejects the half-built residue, so the cap hydrogens are placed
    here by the same internal-coordinate construction as the heavy atoms.
    """
    N, CA, C = atom_xyz(res,"N"), atom_xyz(res,"CA"), atom_xyz(res,"C")
    c   = place(C,  CA, N, 1.335, 121.7, 180.0)   # carbonyl C, extended phi
    o   = place(CA, N,  c, 1.229, 122.9,   0.0)   # O cis to CA  -> trans peptide
    ch3 = place(CA, N,  c, 1.508, 116.6, 180.0)
    h1, h2, h3 = methyl_h(o, c, ch3)
    return [("CH3","C",ch3), ("HH31","H",h1), ("HH32","H",h2), ("HH33","H",h3),
            ("C","C",c), ("O","O",o)]

def nme_for(res):
    """NME capping the C-terminus of res: -NH-CH3 bonded to its C."""
    N, CA, C, O = atom_xyz(res,"N"), atom_xyz(res,"CA"), atom_xyz(res,"C"), atom_xyz(res,"O")
    # Anchor the amide N on the carbonyl O, not on an assumed psi. Placing it at
    # dihedral N-CA-C-Nnew = 180 ignores where O actually sits, and whenever the
    # residue's real psi is near 180 the new N lands on top of O -- it came out at
    # 0.136 A from O/SER541. Referencing O instead fixes N in the peptide plane,
    # anti to CA about the C=O axis, at the standard O-C-N angle.
    n   = place(CA, O,  C, 1.335, 123.0, 180.0)
    h   = place(O,  C,  n, 1.010, 119.0, 180.0)   # amide H anti to the carbonyl O
    ch3 = place(CA, C,  n, 1.449, 121.9, 180.0)   # trans peptide
    h1, h2, h3 = methyl_h(C, n, ch3)
    return [("N","N",n), ("H","H",h), ("CH3","C",ch3),
            ("HH31","H",h1), ("HH32","H",h2), ("HH33","H",h3)]

# ---------------------------------------------------------------- rebuild topology
top = Topology()
newpos = []
CHAIN_IDS = "ABCDEFGH"
for seg_i, (a, b) in enumerate(segments):
    ch = top.addChain(CHAIN_IDS[seg_i])
    cap_n = ace_for(residues[a])
    r = top.addResidue("ACE", ch)
    for nm, el, xyz in cap_n:
        top.addAtom(nm, Element.getBySymbol(el), r); newpos.append(xyz)
    for ri in range(a, b+1):
        src = residues[ri]
        r = top.addResidue(src.name, ch, id=str(resids[ri]))
        for at in src.atoms():
            # PDBFixer adds OXT to what it takes to be the chain terminus. Every segment
            # end is capped with NME here, so a terminal carboxylate on the same residue
            # would be chemically impossible -- the template matcher reports it as
            # "matches CSER, but ... 1 C atom too many".
            if at.name == "OXT":
                continue
            top.addAtom(at.name, at.element, r); newpos.append(pos[at.index])
    cap_c = nme_for(residues[b])
    r = top.addResidue("NME", ch)
    for nm, el, xyz in cap_c:
        top.addAtom(nm, Element.getBySymbol(el), r); newpos.append(xyz)
top.createStandardBonds()

# createStandardBonds() works from heavy-atom templates and creates no bonds to
# hydrogen at all -- it assumes hydrogens come later from Modeller. It also misses
# NME's N-CH3. Since the cap atoms are built here, their bonds are too.
CAP_BONDS = {
    "ACE": [("CH3","HH31"), ("CH3","HH32"), ("CH3","HH33"), ("CH3","C"), ("C","O")],
    "NME": [("N","H"), ("N","CH3"), ("CH3","HH31"), ("CH3","HH32"), ("CH3","HH33")],
}
existing = {frozenset((b[0].index, b[1].index)) for b in top.bonds()}
added = 0
for r in top.residues():
    if r.name not in CAP_BONDS:
        continue
    byname = {a.name: a for a in r.atoms()}
    for n1, n2 in CAP_BONDS[r.name]:
        a1, a2 = byname.get(n1), byname.get(n2)
        if a1 is None or a2 is None:
            continue
        if frozenset((a1.index, a2.index)) not in existing:
            top.addBond(a1, a2); existing.add(frozenset((a1.index, a2.index))); added += 1
print(f"ikatan tutup yang ditambahkan manual: {added}")
# not "heavy atoms": this count includes the 21 cap hydrogens built above
# (3 ACE x 3 H + 3 NME x 4 H). The heavy-atom count is 4,140.
print(f"\nsetelah penutupan: {top.getNumAtoms()} atom, {top.getNumResidues()} residu, "
      f"{top.getNumChains()} rantai")

# ---------------------------------------------------------------- hydrogens
ff = ForceField("amber14-all.xml")
mod = Modeller(top, np.array(newpos)*unit.angstrom)
mod.addHydrogens(ff, pH=PH)
# The 494->498 gap is only 4.05 A across, and both caps were built pointing into it at
# an extended dihedral, so NME and ACE collided (HH32...C at 0.835 A, forces of 3e8
# kJ/mol/nm). The caps are free to rotate away from each other; relax them and nothing
# else. Only the cap atoms move, so the protein -- and the gorge -- stay exactly as the
# crystal left them.
import openmm
sys_cap = ff.createSystem(mod.topology, nonbondedMethod=CutoffNonPeriodic,
                          nonbondedCutoff=1.2*unit.nanometer)
cap_atoms = {a.index for a in mod.topology.atoms() if a.residue.name in ("ACE", "NME")}
for a in mod.topology.atoms():
    if a.index not in cap_atoms:
        sys_cap.setParticleMass(a.index, 0.0)
sim = Simulation(mod.topology, sys_cap, openmm.VerletIntegrator(1*unit.femtosecond),
                     openmm.Platform.getPlatformByName("CPU"))
sim.context.setPositions(mod.positions)
KJ = unit.kilojoule_per_mole
e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
sim.minimizeEnergy(tolerance=1.0*KJ/unit.nanometer, maxIterations=20000)
st = sim.context.getState(getEnergy=True, getPositions=True)
print(f"relaksasi tutup ({len(cap_atoms)} atom bebas): "
      f"{e0:,.0f} -> {st.getPotentialEnergy().value_in_unit(KJ):,.1f} kJ/mol")
mod.positions = st.getPositions()

out = DOCK/f"{PDB}_rec_md.pdb"
with open(out, "w") as fh:
    PDBFile.writeFile(mod.topology, mod.positions, fh, keepIds=True)
print(f"reseptor MD: {mod.topology.getNumAtoms()} atom -> {out.name}")

# ---------------------------------------------------------------- sanity
system = ff.createSystem(mod.topology)
import openmm
for i, f in enumerate(system.getForces()):
    f.setForceGroup(i)
gb = {i for i, f in enumerate(system.getForces()) if isinstance(f, openmm.HarmonicBondForce)}
ctx = openmm.Context(system, openmm.VerletIntegrator(1*unit.femtosecond),
                     openmm.Platform.getPlatformByName("CPU"))
ctx.setPositions(mod.positions)
e = ctx.getState(getEnergy=True, groups=gb).getPotentialEnergy().value_in_unit(
    unit.kilojoule_per_mole)
print(f"energi ikatan: {e:,.1f} kJ/mol   (sebelumnya 185,963 dengan ikatan fiktif)")
