#!/usr/bin/env python
"""Rule-based protonation at pH 7.4 for MD setup.

Written after dimorphite-dl 2.0.2 was found to return tacrine as neutral
(pKa 9.8) and to claim a quaternary ammonium could be uncharged. Explicit
rules are auditable and can be stated in a methods section; a black-box
predictor that is silently wrong cannot.

Rules applied (pH 7.4):
  deprotonated  carboxylic acid, tetrazole, sulfonic/phosphonic acid,
                acylsulfonamide
  protonated    aliphatic amine (1'/2'/3'), amidine, guanidine,
                amidine-like aromatic N (e.g. 4-aminoquinolines: tacrine)
  left neutral  aniline, amide, alcohol, phenol, pyridine, nitrile,
                sulfonamide NH
  preserved     pre-existing formal charges (quaternary ammonium etc.)

Imidazole (pKa ~7) is genuinely ambiguous at 7.4 and is flagged rather
than guessed.
"""
from __future__ import annotations
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# --- groups that lose a proton at 7.4 -------------------------------
# NOTE: each entry holds a LIST of SMARTS. Never comma-join them: a comma is
# also the OR operator *inside* [...] brackets, so splitting on it corrupts
# patterns like [NX3;H2,H1,H0].
ACIDS = [
    ("carboxylic_acid", ["[CX3](=O)[OX2H1]"]),
    ("tetrazole",       ["[nH]1nnnc1", "c1[nH]nnn1", "[nX3H1]1nnnc1"]),
    ("sulfonic_acid",   ["[SX4](=O)(=O)[OX2H1]"]),
    ("phosphonic_acid", ["[PX4](=O)([OX2H1])"]),
    ("acylsulfonamide", ["[CX3](=O)[NX3H1][SX4](=O)(=O)"]),
]
# --- groups that gain a proton at 7.4 -------------------------------
# (name, [smarts], index of the atom in the match that gets the proton)
BASES = [
    # aliphatic amine: exclude amide/sulfonamide N, aniline, imine, nitrile
    ("aliphatic_amine",
     ["[NX3;H2,H1,H0;!$(N~[!#6;!#1]);!$(N-[C,S,P]=[O,S,N]);!$(N-a);"
      "!$(N#*);!$(N=*);!$([N+])]"], 0),
    # amidine: the sp2 N is the basic site, and tertiary amidines (no N-H,
    # e.g. cyclic amidines, pKa ~11) are just as basic as primary ones
    ("amidine",
     ["[NX3;!$(N-[C,S]=[O,S]);!$(N-a)][CX3;!$(C-[O,S])]=[NX2;!$(N-[C,S]=[O,S])]"], 2),
    ("guanidine",
     ["[NX3;!$(N-[C,S]=[O,S])][CX3](=[NX2;!$(N-[C,S]=[O,S])])[NX3]"], 2),
    # 2-/4-aminopyridine-type vinylogous amidines are strongly basic
    # (4-aminopyridine pKa 9.1; tacrine 9.8), unlike plain pyridine (5.2).
    # The exocyclic N must be a real amine: an amide N donates its lone pair
    # to the carbonyl instead, leaving the ring N weakly basic (pKa ~4).
    ("aminoazine_N",
     ["[nX2;$(n:c-[NX3;H2,H1;!$(N-[C,S]=[O,S])]),"
      "$(n:c:c-[NX3;H2,H1;!$(N-[C,S]=[O,S])]),"
      "$(n:c:c:c-[NX3;H2,H1;!$(N-[C,S]=[O,S])])]"], 0),
]
AMBIGUOUS = [("imidazole", ["c1cnc[nH]1", "c1c[nH]cn1"])]


def _match(mol, smarts_list):
    """First non-empty match across a list of SMARTS alternatives."""
    for s in smarts_list:
        patt = Chem.MolFromSmarts(s)
        if patt is None:
            continue
        hits = mol.GetSubstructMatches(patt)
        if hits:
            return hits
    return ()


def protonate(smiles: str, ph: float = 7.4) -> dict:
    """Return dominant state at pH, its charge, applied rules, and ambiguity."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"smiles": smiles, "charge": 0, "rules": [], "ambiguous": True,
                "error": "unparseable"}
    mol = Chem.RWMol(mol)
    applied, ambiguous = [], False

    # acids -> anion
    for name, smarts in ACIDS:
        for hit in _match(mol, smarts):
            for idx in hit:
                a = mol.GetAtomWithIdx(idx)
                if a.GetSymbol() in ("O", "N") and a.GetTotalNumHs() > 0 \
                        and a.GetFormalCharge() == 0:
                    a.SetFormalCharge(-1)
                    a.SetNumExplicitHs(max(0, a.GetTotalNumHs() - 1))
                    a.SetNoImplicit(True)
                    applied.append(name)
                    break

    # bases -> cation
    for name, smarts, site in BASES:
        for hit in _match(mol, smarts):
            if site >= len(hit):
                continue
            idx = hit[site]
            a = mol.GetAtomWithIdx(idx)
            if a.GetSymbol().lower() != "n" or a.GetFormalCharge() != 0:
                continue
            if a.GetIsAromatic():
                if a.GetTotalNumHs() != 0 or a.GetDegree() != 2:
                    continue          # only pyridine-type ring N
                a.SetFormalCharge(1)
                a.SetNumExplicitHs(1)
            else:
                if a.GetTotalValence() >= 4:
                    continue          # already quaternary
                a.SetFormalCharge(1)
                a.SetNumExplicitHs(a.GetTotalNumHs() + 1)
            a.SetNoImplicit(True)
            applied.append(name)

    for name, smarts in AMBIGUOUS:
        if _match(mol, smarts):
            ambiguous = True
            applied.append(f"{name}(ambiguous)")

    try:
        out = mol.GetMol()
        Chem.SanitizeMol(out)
        return {"smiles": Chem.MolToSmiles(out),
                "charge": Chem.GetFormalCharge(out),
                "rules": sorted(set(applied)), "ambiguous": ambiguous,
                "error": None}
    except Exception as e:
        m0 = Chem.MolFromSmiles(smiles)
        return {"smiles": smiles, "charge": Chem.GetFormalCharge(m0),
                "rules": [], "ambiguous": True, "error": str(e)[:60]}


VALIDATION = [
    ("tacrine",       "Nc1c2ccccc2nc2c1CCCC2", 1),
    ("donepezil",     "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2", 1),
    ("galantamine",   "CN1CCC23C=CC(O)CC2Oc2c(OC)ccc(c23)C1", 1),
    ("rivastigmine",  "CCN(C)C(=O)Oc1cccc(C(C)N(C)C)c1", 1),
    ("edrophonium",   "CC[N+](C)(C)c1cccc(O)c1", 1),
    ("huperzine A",   "CC=C1C2CC3=C(C=CC(=O)N3)C1(N)CC(C)=C2", 1),
    ("aspirin",       "CC(=O)Oc1ccccc1C(=O)O", -1),
    ("phenol",        "Oc1ccccc1", 0),
    ("aniline",       "Nc1ccccc1", 0),
    ("acetamide",     "CC(N)=O", 0),
    ("pyridine",      "c1ccncc1", 0),
    ("benzenesulfonamide", "NS(=O)(=O)c1ccccc1", 0),
    # cyclic amidine, pKa ~11 - no N-H, still basic
    ("cyclic amidine",  "c1ccc2c(c1)CN1CCCCC1=N2", 1),
    # 2-acylaminopyridine: amide N is not a donor, ring N stays weak (pKa ~4)
    ("2-acetamidopyridine", "CC(=O)Nc1ccccn1", 0),
    ("2-aminopyridine",     "Nc1ccccn1", 1),
    ("DBU (amidine)",       "C1CCC2=NCCCN2CC1", 1),
]

if __name__ == "__main__":
    print(f"{'compound':22} {'expect':>7} {'got':>5}  {'':2} rules")
    print("-" * 72)
    n_ok = 0
    for name, smi, expect in VALIDATION:
        r = protonate(smi)
        ok = r["charge"] == expect
        n_ok += ok
        print(f"{name:22} {expect:>7} {r['charge']:>5}  {'OK' if ok else 'XX'} "
              f"{','.join(r['rules']) or '-'}")
    print(f"\n{n_ok}/{len(VALIDATION)} correct")
