#!/usr/bin/env python
"""Screen ChEMBL targets for activity-cliff study feasibility.

Answers one question per target: does it have enough well-measured compounds,
spanning enough potency, forming enough activity cliffs, to support the study?

  python chembl_screen.py --search "beta-secretase"      # find target IDs
  python chembl_screen.py --targets CHEMBL203 CHEMBL244  # screen them
  python chembl_screen.py --preset kinases

Raw downloads are cached under data/chembl_cache/, so re-runs are free.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, rdMMPA, SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.DataStructs import BulkTanimotoSimilarity

RDLogger.DisableLog("rdApp.*")

API = "https://www.ebi.ac.uk/chembl/api/data"
ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "chembl_cache"
RESULTS = ROOT / "results"

# ---- feasibility thresholds (the study's entry criteria) ------------
MIN_COMPOUNDS   = 80     # enough to train + hold out a cliff-aware test set
MIN_LOG_SPAN    = 3.0    # pChEMBL range, log units
MIN_CLIFF_PAIRS = 25     # strict cliffs needed to make claims about them
SIM_CUTOFF      = 0.90   # ECFP4 Tanimoto: "near-identical structure"
CLIFF_DELTA     = 2.0    # 100-fold potency change
SOFT_DELTA      = 1.0    # 10-fold, reported for context
REPLICATE_SD    = 1.0    # discard compounds whose repeat measurements disagree

PRESETS = {
    "kinases":   ["CHEMBL203", "CHEMBL279", "CHEMBL267", "CHEMBL2971"],
    "proteases": ["CHEMBL244", "CHEMBL4822"],
    "cns":       ["CHEMBL220", "CHEMBL4822"],
}


# ------------------------------------------------------------------ API
def _get(endpoint: str, params: dict, tolerant: bool = False) -> list[dict]:
    """Paginated ChEMBL REST fetch.

    ChEMBL 500s on deep offsets (beyond ~16k), so callers should filter
    server-side to keep result sets small. With tolerant=True a page that
    keeps failing ends the fetch with partial data instead of crashing.
    """
    out, offset, limit = [], 0, 1000
    while True:
        p = dict(params, limit=limit, offset=offset, format="json")
        data = None
        for attempt in range(6):
            try:
                r = requests.get(f"{API}/{endpoint}.json", params=p, timeout=120)
                r.raise_for_status()
                data = r.json()
                break
            except requests.RequestException:
                if attempt < 5:
                    time.sleep(min(45, 3 * 2 ** attempt))
        if data is None:
            if tolerant:
                print(f"      warning: page at offset {offset:,} failed; "
                      f"continuing with {len(out):,} records", file=sys.stderr)
                break
            raise RuntimeError(f"ChEMBL fetch failed at offset {offset}")
        key = [k for k in data if k != "page_meta"][0]
        batch = data[key]
        out.extend(batch)
        total = data["page_meta"]["total_count"]
        offset += limit
        print(f"      fetched {min(offset, total):,}/{total:,}", end="\r", file=sys.stderr)
        if offset >= total or not batch:
            break
        time.sleep(0.3)
    print(" " * 44, end="\r", file=sys.stderr)
    return out


def search_targets(query: str) -> pd.DataFrame:
    rows = _get("target", {"pref_name__icontains": query,
                           "target_type": "SINGLE PROTEIN"})
    return pd.DataFrame([{
        "target_chembl_id": t["target_chembl_id"],
        "pref_name": t["pref_name"],
        "organism": t.get("organism"),
    } for t in rows])


def fetch_activities(tid: str, refresh: bool = False) -> pd.DataFrame:
    """Fetch binding activities, filtered server-side to dodge deep paging."""
    cache = CACHE / f"{tid}.csv.gz"
    if cache.exists() and not refresh:
        return pd.read_csv(cache, low_memory=False)
    print(f"    downloading {tid} from ChEMBL...", file=sys.stderr)
    keep = ["molecule_chembl_id", "canonical_smiles", "standard_type",
            "standard_value", "standard_units", "standard_relation",
            "pchembl_value", "assay_chembl_id", "data_validity_comment",
            "potential_duplicate", "target_pref_name"]
    rows = []
    # split by assay type: each subset stays well under the paging limit,
    # and requiring pChEMBL server-side removes most junk before transfer
    for stype in ("IC50", "Ki", "Kd"):
        got = _get("activity", {
            "target_chembl_id": tid,
            "standard_type": stype,
            "assay_type": "B",
            "pchembl_value__isnull": "false",
            "standard_relation": "=",
        }, tolerant=True)
        if got:
            print(f"      {stype}: {len(got):,}", file=sys.stderr)
        rows.extend(got)
    df = pd.DataFrame([{k: a.get(k) for k in keep} for a in rows])
    CACHE.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False, compression="gzip")
    return df


# -------------------------------------------------------------- cleaning
_saltremover = SaltRemover.SaltRemover()
_uncharger = rdMolStandardize.Uncharger()


def standardize(smiles) -> str | None:
    # NaN from the CSV cache is a float and is truthy - guard on type
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = _saltremover.StripMol(mol, dontRemoveEverything=True)
        mol = rdMolStandardize.Cleanup(mol)
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if not frags:
            return None
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())  # keep parent
        mol = _uncharger.uncharge(mol)
        if mol.GetNumHeavyAtoms() < 6:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


def clean(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Quality filters. A cliff built on bad data is an artifact, not a finding."""
    audit = {"raw_records": len(df)}
    if df.empty:
        return df, audit

    d = df.copy()
    d = d[d["standard_relation"] == "="]                     # censored data out
    audit["exact_relation"] = len(d)
    d = d[d["standard_units"] == "nM"]
    d = d[d["pchembl_value"].notna()]
    audit["has_pchembl"] = len(d)
    d = d[d["data_validity_comment"].isna()]                 # ChEMBL's own flags
    audit["validity_ok"] = len(d)
    d = d[d["potential_duplicate"].fillna(0).astype(int) == 0]
    d["pchembl_value"] = d["pchembl_value"].astype(float)

    # one assay type only - mixing IC50 with Ki fabricates cliffs
    if len(d):
        best = d["standard_type"].value_counts().idxmax()
        d = d[d["standard_type"] == best]
        audit["assay_type"] = best
    audit["single_assay_type"] = len(d)

    d["smiles_std"] = d["canonical_smiles"].map(standardize)
    d = d[d["smiles_std"].notna()]
    audit["standardized"] = len(d)

    # collapse replicates; drop compounds whose repeats disagree badly
    g = d.groupby("smiles_std")["pchembl_value"]
    agg = pd.DataFrame({"pchembl": g.median(), "n_meas": g.size(),
                        "sd": g.std().fillna(0.0)})
    agg = agg[agg["sd"] <= REPLICATE_SD]
    audit["inconsistent_dropped"] = int((g.std().fillna(0) > REPLICATE_SD).sum())

    ids = d.groupby("smiles_std")["molecule_chembl_id"].first()
    out = agg.join(ids).reset_index()
    audit["final_compounds"] = len(out)
    return out, audit


# ----------------------------------------------------------- cliff finding
_fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def _flat(smiles) -> str | None:
    """Constitution only - stereochemistry stripped."""
    if not isinstance(smiles, str) or not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m, isomericSmiles=False) if m else None


def find_cliffs(df: pd.DataFrame) -> pd.DataFrame:
    """Pairs that are structurally near-identical but far apart in potency."""
    mols = [Chem.MolFromSmiles(s) for s in df["smiles_std"]]
    fps = [_fpgen.GetFingerprint(m) for m in mols]
    pot = df["pchembl"].to_numpy()
    sm = df["smiles_std"].tolist()
    cid = df["molecule_chembl_id"].tolist()
    flat = [_flat(x) for x in sm]

    rows = []
    for i in range(len(fps) - 1):
        sims = np.array(BulkTanimotoSimilarity(fps[i], fps[i + 1:]))
        d = np.abs(pot[i + 1:] - pot[i])
        hits = np.where((sims >= SIM_CUTOFF) & (d >= SOFT_DELTA))[0]
        for h in hits:
            j = i + 1 + h
            s = float(sims[h])
            rows.append({
                "cid_a": cid[i], "cid_b": cid[j],
                "smiles_a": sm[i], "smiles_b": sm[j],
                "pchembl_a": pot[i], "pchembl_b": pot[j],
                "delta": float(d[h]), "tanimoto": s,
                "sali": float(d[h] / (1 - s)) if s < 1 else np.inf,
                "strict": bool(d[h] >= CLIFF_DELTA),
                # ECFP4 ignores stereochemistry, so enantiomer pairs score
                # T=1.0; they are a different phenomenon from constitutional
                # cliffs and must be counted separately
                "cliff_class": ("stereochemical"
                                if flat[i] is not None and flat[i] == flat[j]
                                else "constitutional"),
            })
    return pd.DataFrame(rows).sort_values("sali", ascending=False) if rows else pd.DataFrame()


def mmp_cliffs(df: pd.DataFrame, cap: int = 600) -> int:
    """Matched molecular pairs (single cut) that are also cliffs - the strict
    medicinal-chemistry definition: identical except one substituent."""
    sub = df.head(cap)
    index: dict[str, list] = {}
    for _, r in sub.iterrows():
        mol = Chem.MolFromSmiles(r["smiles_std"])
        if mol is None:
            continue
        try:
            frags = rdMMPA.FragmentMol(mol, maxCuts=1, resultsAsMols=False)
        except Exception:
            continue
        for core, chains in frags:
            parts = [p for p in (chains or "").split(".") if p]
            if len(parts) != 2:
                continue
            scaffold, rgroup = sorted(parts, key=len, reverse=True)
            index.setdefault(scaffold, []).append((r["pchembl"], rgroup))

    n = 0
    for _, members in index.items():
        if len(members) < 2:
            continue
        for a in range(len(members) - 1):
            for b in range(a + 1, len(members)):
                if members[a][1] == members[b][1]:
                    continue
                if abs(members[a][0] - members[b][0]) >= CLIFF_DELTA:
                    n += 1
    return n


# ------------------------------------------------------------------ report
def screen(tid: str, refresh: bool = False) -> dict:
    print(f"\n>> {tid}")
    raw = fetch_activities(tid, refresh)
    if raw.empty:
        print("   no activity records")
        return {"target_chembl_id": tid, "verdict": "NO DATA"}

    name = raw["target_pref_name"].dropna().iloc[0] if raw["target_pref_name"].notna().any() else "?"
    df, audit = clean(raw)
    print(f"   {name}")
    print(f"   records {audit['raw_records']:,} -> exact= {audit.get('exact_relation',0):,}"
          f" -> pChEMBL {audit.get('has_pchembl',0):,}"
          f" -> valid {audit.get('validity_ok',0):,}"
          f" -> {audit.get('assay_type','?')} only {audit.get('single_assay_type',0):,}"
          f" -> compounds {audit.get('final_compounds',0):,}")
    if audit.get("inconsistent_dropped"):
        print(f"   dropped {audit['inconsistent_dropped']} compounds with "
              f"discordant replicates (SD > {REPLICATE_SD} log)")

    n = len(df)
    if n < 10:
        print("   too few compounds to analyse")
        return {"target_chembl_id": tid, "pref_name": name, "n_compounds": n,
                "verdict": "FAIL"}

    span = float(df["pchembl"].max() - df["pchembl"].min())
    cliffs = find_cliffs(df)
    strict = cliffs[cliffs["strict"]] if not cliffs.empty else pd.DataFrame()
    n_strict = len(strict)
    n_soft = len(cliffs)
    if n_strict:
        con = strict[strict["cliff_class"] == "constitutional"]
        ste = strict[strict["cliff_class"] == "stereochemical"]
    else:
        con = ste = pd.DataFrame()
    n_con, n_ste = len(con), len(ste)
    involved = (len(set(con["cid_a"]) | set(con["cid_b"])) if n_con else 0)
    involved_ste = (len(set(ste["cid_a"]) | set(ste["cid_b"])) if n_ste else 0)
    n_mmp = mmp_cliffs(df)

    # verdict rests on CONSTITUTIONAL cliffs - the classic QSAR failure mode
    ok = (n >= MIN_COMPOUNDS and span >= MIN_LOG_SPAN and n_con >= MIN_CLIFF_PAIRS)
    verdict = "PASS" if ok else "FAIL"

    print(f"   compounds {n:,} | pChEMBL {df['pchembl'].min():.1f}-{df['pchembl'].max():.1f}"
          f" (span {span:.1f})")
    print(f"   strict cliffs ({n_strict}): {n_con} constitutional "
          f"({involved} compounds), {n_ste} stereochemical ({involved_ste} compounds)")
    print(f"   soft cliffs (>={SOFT_DELTA} log): {n_soft} | MMP-based: {n_mmp}")
    print(f"   {verdict}")

    if n_strict:
        strict.to_csv(RESULTS / f"01_cliffs_{tid}.csv", index=False)
        df.to_csv(RESULTS / f"01_compounds_{tid}.csv", index=False)

    return {"target_chembl_id": tid, "pref_name": name, "n_compounds": n,
            "assay_type": audit.get("assay_type"), "pchembl_span": round(span, 2),
            "cliffs_constitutional": n_con, "cliffs_stereochemical": n_ste,
            "cmpds_constitutional": involved, "cmpds_stereochemical": involved_ste,
            "cliff_pairs_soft": n_soft, "mmp_cliffs": n_mmp, "verdict": verdict}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+")
    ap.add_argument("--preset", choices=PRESETS)
    ap.add_argument("--search")
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)

    if a.search:
        res = search_targets(a.search)
        print(res.to_string(index=False) if not res.empty else "no targets found")
        return

    targets = a.targets or (PRESETS[a.preset] if a.preset else None)
    if not targets:
        ap.error("give --targets, --preset or --search")

    rows = [screen(t, a.refresh) for t in targets]
    summary = pd.DataFrame(rows)
    out = RESULTS / "01_feasibility.csv"
    summary.to_csv(out, index=False)
    print("\n" + "=" * 78)
    print(summary.to_string(index=False))
    print(f"\ncriteria: >={MIN_COMPOUNDS} compounds, >={MIN_LOG_SPAN} log span, "
          f">={MIN_CLIFF_PAIRS} CONSTITUTIONAL cliffs")
    print(f"written: {out}")


if __name__ == "__main__":
    main()
