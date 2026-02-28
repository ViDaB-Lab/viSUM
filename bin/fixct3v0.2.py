#!/usr/bin/env python3
"""
fixct3v0.2.py

Purpose
-------
Convert Cenote-Taker3 {run}_virus_summary.tsv into a viSUM-style evidence CSV
with ICTV ranks + hallmark counts:

seqid,
d__Domain,r__Realm,k__Kingdom,p__Phylum,c__Class,o__Order,f__Family,g__Genus,s__Species,
ct3_virion_hallmark_count,ct3_rep_hallmark_count,ct3_RDRP_hallmark_count

Key behaviors
-------------
1) Parse CT3 `taxonomy_hierarchy` which uses tokens like:
   d_Viruses;-_Varidnaviria;k_Bamfordvirae;p_Preplasmiviricota;c_Polintoviricetes;o_Orthopolintovirales

   Notes:
   - domain token is usually `d_...`
   - realm token is often `-_...`  (yes, dash!)
   - then k_, p_, c_, o_, f_, g_, s_ may appear

2) Missing ranks are filled as 'unclassified'.

3) ICTV fallback (optional, default: on):
   If CT3 provides only weak taxonomy (e.g., realm is "unclassified virus"),
   try to map the deepest available label to a full ICTV lineage using the
   ICTV master species list (same concept as fixgenomadv0.2) :contentReference[oaicite:4]{index=4}

4) Always writes a CSV with header, even if CT3 summary is empty (pipeline-safe).

Usage
-----
python3 fixct3v0.2.py --ct3 run_virus_summary.tsv --ictv ICTV.csv --out prefix.ct3.vsum.csv

"""

import argparse
from typing import Dict, Optional, Tuple, List

import pandas as pd


# ----------------------------
# Helpers
# ----------------------------

ICTV_RANKS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

OUT_COLS = [
    "seqid",
    "d__Domain", "r__Realm", "k__Kingdom", "p__Phylum", "c__Class", "o__Order", "f__Family", "g__Genus", "s__Species",
    "ct3_virion_hallmark_count", "ct3_rep_hallmark_count", "ct3_RDRP_hallmark_count",
]


def clean_label(x: Optional[str]) -> str:
    """Normalize empty/NA-ish tokens to literal 'unclassified'."""
    if x is None:
        return "unclassified"
    s = str(x).strip()
    if s == "":
        return "unclassified"
    if s.lower() in {"na", "nan", "none"}:
        return "unclassified"
    return s


def lc(x: Optional[str]) -> str:
    """Lowercased clean label for matching."""
    return clean_label(x).lower()


def to_int_or_zero(x: Optional[str]) -> int:
    """CT3 counts might be NaN; convert safely to int."""
    try:
        if x is None:
            return 0
        if pd.isna(x):
            return 0
        return int(float(x))
    except Exception:
        return 0


# ----------------------------
# ICTV indexing (same idea as fixgenomadv0.2)
# ----------------------------

def build_ictv_index(ictv_csv: str) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    index[rank][lower(label)] -> full lineage dict {realm..species}
    """
    df = pd.read_csv(ictv_csv)
    df.columns = [c.lower() for c in df.columns]

    for r in ICTV_RANKS:
        if r not in df.columns:
            df[r] = ""

    index: Dict[str, Dict[str, Dict[str, str]]] = {r: {} for r in ICTV_RANKS}

    for _, row in df.fillna("").iterrows():
        lineage = {r: clean_label(row.get(r, "")) for r in ICTV_RANKS}
        for r in ICTV_RANKS:
            key = lc(lineage[r])
            if key != "unclassified":
                index[r][key] = lineage

    return index


def ictv_lookup(index: Dict[str, Dict[str, Dict[str, str]]], label: str) -> Optional[Dict[str, str]]:
    """Search most specific -> broadest."""
    key = lc(label)
    if key == "unclassified":
        return None

    for r in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]:
        if key in index[r]:
            return index[r][key]
    return None


# ----------------------------
# CT3 taxonomy parsing
# ----------------------------

CT3_PREFIX_TO_RANK = {
    "d_": "domain",
    "-_": "realm",   # CT3 uses -_ for realm (e.g. -_Varidnaviria)
    "r_": "realm",   # just in case a future CT3 version uses r_
    "k_": "kingdom",
    "p_": "phylum",
    "c_": "class",
    "o_": "order",
    "f_": "family",
    "g_": "genus",
    "s_": "species",
}

def parse_ct3_taxonomy_hierarchy(taxstr: Optional[str]) -> Dict[str, str]:
    """
    Parse CT3 taxonomy_hierarchy and return clean rank strings WITHOUT prefixes,
    e.g. {"domain":"Viruses","realm":"Varidnaviria",...}
    """
    out = {r: "unclassified" for r in ["domain"] + ICTV_RANKS}
    out["domain"] = "Viruses"

    if taxstr is None or (isinstance(taxstr, float) and pd.isna(taxstr)):
        return out

    s = str(taxstr).strip()
    if not s:
        return out

    parts = [p.strip() for p in s.split(";") if p.strip()]
    for p in parts:
        if len(p) < 3:
            continue
        pref = p[:2]
        if pref not in CT3_PREFIX_TO_RANK:
            continue
        rank = CT3_PREFIX_TO_RANK[pref]
        val = clean_label(p[2:])
        # CT3 sometimes uses "unclassified virus"
        if val.lower() == "unclassified virus":
            val = "unclassified"
        out[rank] = val

    return out


def deepest_non_unclassified(rankdict: Dict[str, str]) -> str:
    """
    Choose deepest informative label from CT3 ranks:
    species -> genus -> family -> ... -> realm
    """
    for r in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]:
        v = clean_label(rankdict.get(r, "unclassified"))
        if lc(v) != "unclassified":
            return v
    return "unclassified"


def fill_missing_with_ictv(rankdict: Dict[str, str], ictv_lineage: Dict[str, str]) -> Dict[str, str]:
    """
    Only fill ranks that are currently unclassified; never overwrite CT3 detail.
    """
    out = dict(rankdict)
    for r in ICTV_RANKS:
        if lc(out.get(r)) == "unclassified":
            out[r] = clean_label(ictv_lineage.get(r, "unclassified"))
    return out


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ct3", required=True, help="CT3 {run}_virus_summary.tsv")
    ap.add_argument("--ictv", required=True, help="ICTV master species list CSV")
    ap.add_argument("--out", required=True, help="Output evidence CSV")
    ap.add_argument("--fallback", choices=["none", "ictv"], default="ictv",
                    help="Fallback behavior if CT3 taxonomy is sparse/unclassified (default: ictv)")
    args = ap.parse_args()

    ictv_index = build_ictv_index(args.ictv)

    ct3 = pd.read_csv(args.ct3, sep="\t")
    ct3.columns = [c.strip() for c in ct3.columns]

    # Prefer column names over positions
    required = ["input_name", "taxonomy_hierarchy", "virion_hallmark_count", "rep_hallmark_count", "RDRP_hallmark_count"]
    missing = [c for c in required if c not in ct3.columns]
    if missing:
        raise ValueError(f"CT3 summary missing required columns: {missing}. Found: {list(ct3.columns)}")

    rows: List[Dict[str, object]] = []

    stats = {
        "n_total": 0,
        "n_used_ct3_tax": 0,
        "n_filled_by_ictv": 0,
        "n_unclassified_final": 0,
    }

    for _, r in ct3.iterrows():
        stats["n_total"] += 1

        seqid = str(r["input_name"]).strip()

        # Parse taxonomy hierarchy
        rankdict = parse_ct3_taxonomy_hierarchy(r.get("taxonomy_hierarchy"))

        # Decide whether taxonomy is "sparse" (i.e., basically unclassified)
        # If realm is unclassified and no family/order/class etc, we try fallback.
        has_any_detail = any(lc(rankdict.get(x)) != "unclassified" for x in ICTV_RANKS)
        if has_any_detail:
            stats["n_used_ct3_tax"] += 1

        if args.fallback == "ictv":
            label = deepest_non_unclassified(rankdict)
            if lc(label) != "unclassified":
                L = ictv_lookup(ictv_index, label)
                if L:
                    before = dict(rankdict)
                    rankdict = fill_missing_with_ictv(rankdict, L)
                    # count as filled only if something changed
                    if any(before[k] != rankdict[k] for k in ICTV_RANKS):
                        stats["n_filled_by_ictv"] += 1

        # Final check for totally unclassified (beyond domain)
        if all(lc(rankdict.get(x)) == "unclassified" for x in ICTV_RANKS):
            stats["n_unclassified_final"] += 1

        out = {
            "seqid": seqid,
            "d__Domain": f"d__{clean_label(rankdict.get('domain', 'Viruses'))}",
            "r__Realm": f"r__{clean_label(rankdict.get('realm'))}",
            "k__Kingdom": f"k__{clean_label(rankdict.get('kingdom'))}",
            "p__Phylum": f"p__{clean_label(rankdict.get('phylum'))}",
            "c__Class": f"c__{clean_label(rankdict.get('class'))}",
            "o__Order": f"o__{clean_label(rankdict.get('order'))}",
            "f__Family": f"f__{clean_label(rankdict.get('family'))}",
            "g__Genus": f"g__{clean_label(rankdict.get('genus'))}",
            "s__Species": f"s__{clean_label(rankdict.get('species'))}",
            "ct3_virion_hallmark_count": to_int_or_zero(r.get("virion_hallmark_count")),
            "ct3_rep_hallmark_count": to_int_or_zero(r.get("rep_hallmark_count")),
            "ct3_RDRP_hallmark_count": to_int_or_zero(r.get("RDRP_hallmark_count")),
        }
        rows.append(out)

    out_df = pd.DataFrame(rows, columns=OUT_COLS)
    out_df.to_csv(args.out, index=False)

    print("✅ Wrote:", args.out)
    print("📊 Stats:", stats)


if __name__ == "__main__":
    main()