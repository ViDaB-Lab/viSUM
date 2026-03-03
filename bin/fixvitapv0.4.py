#!/usr/bin/env python3
"""
fixvitapv0.4.py

Purpose
-------
Convert VITAP v1.10 'best_determined_lineages.tsv' into a viSUM-style evidence table.

What's new vs v0.3
------------------
1) Supports VITAP v1.10 column name:
     'lineage_score/participation_index'
   (older builds sometimes used 'lineage_score').

2) Uses --header_map (preferred) to filter ONLY query contigs.
   Still supports --seqids (one-per-line allowlist) for legacy usage.

3) Outputs a stable viSUM schema with additional VITAP metadata columns:
     vitap_pi, vitap_confidence_level

4) Robust lineage parsing:
   VITAP lineage string is typically:
     Species;Genus;Family;Order;Class;Phylum;Kingdom;Realm
   with '-' placeholders. We assign ranks from the RIGHT.

5) Optional ICTV fallback:
   If lineage is too uninformative (mostly '-' / unclassified),
   attempt to fill ranks from ICTV master species list.

Output
------
CSV columns:
  seqid,d__Domain,r__Realm,k__Kingdom,p__Phylum,c__Class,o__Order,f__Family,g__Genus,s__Species,vitap_pi,vitap_confidence_level

Notes
-----
- Domain is always d__Viruses for rows we keep (VITAP is used as viral taxonomy evidence).
- We intentionally do NOT try to force genus/species if VITAP doesn't provide them.
- "vitap_pi" is whatever VITAP reports (lineage_score/participation_index); typically 0..1.
"""

import argparse
from typing import Dict, Optional, List, Set
import pandas as pd


VIRUS_DOMAIN = "Viruses"
RANKS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

COLS_OUT = [
    "seqid",
    "d__Domain",
    "r__Realm",
    "k__Kingdom",
    "p__Phylum",
    "c__Class",
    "o__Order",
    "f__Family",
    "g__Genus",
    "s__Species",
    "vitap_pi",
    "vitap_confidence_level",
]


# ----------------------------
# Normalization helpers
# ----------------------------
def clean_token(x: Optional[str]) -> str:
    """Normalize empty/NA/dash tokens to literal 'unclassified'."""
    if x is None:
        return "unclassified"
    x = str(x).strip()
    if x == "" or x.lower() in {"-", "na", "nan", "none"}:
        return "unclassified"
    # Defensive: remove rare bracketed prefixes if any tool emits them
    if "]_" in x:
        x = x.split("]_")[-1].strip()
    return x


def lc(x: Optional[str]) -> str:
    return clean_token(x).lower()


# ----------------------------
# Allowlist readers
# ----------------------------
def read_seqids_file(seqids_path: str) -> Set[str]:
    """Read one seqid per line."""
    keep: Set[str] = set()
    with open(seqids_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                keep.add(s)
    return keep


def read_header_map(header_map_tsv: str) -> Set[str]:
    """
    Read viSUM header_map.tsv and build an allowlist.

    Expected columns from your NORMALIZE_FASTA:
      prefix, new_id, original_id, original_header, length

    We include:
      - new_id (definitely used downstream)
      - original_id (common FASTA first token)
      - optional: first token of original_header (extra safety)
    """
    hm = pd.read_csv(header_map_tsv, sep="\t", dtype=str).fillna("")
    hm.columns = [c.strip().lower() for c in hm.columns]

    # Robustly detect columns
    new_id_col = "new_id" if "new_id" in hm.columns else None
    orig_id_col = "original_id" if "original_id" in hm.columns else None
    orig_header_col = "original_header" if "original_header" in hm.columns else None

    if not new_id_col:
        raise ValueError(f"header_map missing required column 'new_id'. Found: {list(hm.columns)}")

    keep: Set[str] = set(hm[new_id_col].astype(str).str.strip().tolist())

    if orig_id_col:
        keep |= set(hm[orig_id_col].astype(str).str.strip().tolist())

    if orig_header_col:
        # Add the first token of the original header (FASTA ID)
        toks = hm[orig_header_col].astype(str).str.strip().str.split().str[0].fillna("")
        keep |= set(toks.tolist())

    # Remove empties
    keep.discard("")
    return keep


# ----------------------------
# ICTV indexing (optional fallback)
# ----------------------------
def build_ictv_index(ictv_csv: str) -> Dict[str, Dict[str, str]]:
    """
    Build a label->lineage dict for quick lookup by any rank label (lowercased).

    The ICTV CSV is expected to have columns named like:
      realm, kingdom, phylum, class, order, family, genus, species
    If some are missing, we create them as empty.
    """
    ictv = pd.read_csv(ictv_csv).fillna("")
    ictv.columns = [c.strip().lower() for c in ictv.columns]
    for r in RANKS:
        if r not in ictv.columns:
            ictv[r] = ""

    index: Dict[str, Dict[str, str]] = {}
    for _, row in ictv.iterrows():
        lineage = {r: clean_token(row.get(r, "")) for r in RANKS}
        # Store mapping by each rank label (most specific entries will exist too)
        for r in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]:
            name = lc(lineage.get(r))
            if name != "unclassified":
                index.setdefault(name, lineage)
    return index


def ictv_lookup(index: Dict[str, Dict[str, str]], label: str) -> Optional[Dict[str, str]]:
    label = lc(label)
    if label == "unclassified":
        return None
    return index.get(label)


# ----------------------------
# VITAP lineage parsing
# ----------------------------
def parse_vitap_lineage(raw: Optional[str]) -> Dict[str, str]:
    """
    Parse VITAP lineage into ranks.

    VITAP v1.10 typical example:
      -;-;-;Ortervirales;Revtraviricetes;Artverviricota;Pararnavirae;Riboviria

    Which corresponds to:
      species,genus,family,order,class,phylum,kingdom,realm

    We split by ';', clean tokens, then assign from RIGHT:
      realm <- last token
      kingdom <- second last
      ...
      species <- first (if present)
    """
    raw = "" if raw is None else str(raw)
    toks = [clean_token(t) for t in raw.split(";")]

    out = {r: "unclassified" for r in RANKS}
    if not toks:
        return out

    right = toks[::-1]  # realm-ish first
    for i, r in enumerate(RANKS):  # realm..species
        if i < len(right):
            out[r] = clean_token(right[i])
    return out


def best_label_for_fallback(parsed: Dict[str, str]) -> str:
    """Pick the most specific non-unclassified label."""
    for r in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]:
        if lc(parsed.get(r)) != "unclassified":
            return parsed[r]
    return "unclassified"


def should_fallback(parsed: Dict[str, str]) -> bool:
    """
    Decide whether parsed lineage is too uninformative.
    Heuristic:
      - if <= 1 informative rank, fallback
    """
    informative = [r for r in RANKS if lc(parsed.get(r)) != "unclassified"]
    return len(informative) <= 1


# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vitap", required=True, help="VITAP best_determined_lineages.tsv")
    ap.add_argument("--out", required=True, help="Output evidence CSV")

    # Filtering inputs (choose one)
    ap.add_argument("--header_map", help="header_map.tsv from NORMALIZE_FASTA (preferred)")
    ap.add_argument("--seqids", help="Allowlist file (one seqid per line)")

    # Optional fallback
    ap.add_argument("--ictv", help="ICTV master species list CSV (for fallback)")
    ap.add_argument("--fallback", choices=["none", "ictv"], default="none",
                    help="Fallback method when parsed lineage is weak (default: none)")

    args = ap.parse_args()

    if not args.header_map and not args.seqids:
        raise ValueError("Provide either --header_map or --seqids so we can filter to query contigs.")

    if args.fallback == "ictv" and not args.ictv:
        raise ValueError("--fallback ictv requires --ictv <ICTV master species list CSV>")

    # Build allowlist
    if args.header_map:
        keep = read_header_map(args.header_map)
    else:
        keep = read_seqids_file(args.seqids)

    # ICTV index if needed
    ictv_index = build_ictv_index(args.ictv) if args.fallback == "ictv" else {}

    df = pd.read_csv(args.vitap, sep="\t", dtype=str).fillna("")

    # Normalize column names
    cols = {c.strip().lower(): c for c in df.columns}
    id_col  = cols.get("genome_id")
    lin_col = cols.get("lineage")

    # Score column has changed names across versions
    score_col = (
        cols.get("lineage_score/participation_index")
        or cols.get("lineage_score")
        or cols.get("lineage_score_participation_index")
    )
    conf_col = cols.get("confidence_level")

    if not id_col or not lin_col:
        raise ValueError(f"Expected columns Genome_ID and lineage. Found: {list(df.columns)}")

    out_rows: List[Dict[str, str]] = []
    stats = {
        "n_total_rows": 0,
        "n_kept_by_allowlist": 0,
        "n_fallback_attempted": 0,
        "n_fallback_filled": 0,
    }

    for _, row in df.iterrows():
        stats["n_total_rows"] += 1
        seqid = str(row.get(id_col, "")).strip()

        # Filter out reference genomes / non-query rows
        if seqid not in keep:
            continue
        stats["n_kept_by_allowlist"] += 1

        raw_lineage = row.get(lin_col, "")
        parsed = parse_vitap_lineage(raw_lineage)

        # Optional ICTV fallback
        if args.fallback == "ictv" and should_fallback(parsed):
            stats["n_fallback_attempted"] += 1
            label = best_label_for_fallback(parsed)
            L = ictv_lookup(ictv_index, label)
            if L:
                parsed = {r: clean_token(L.get(r, "unclassified")) for r in RANKS}
                stats["n_fallback_filled"] += 1

        vitap_pi = row.get(score_col, "") if score_col else ""
        vitap_conf = row.get(conf_col, "") if conf_col else ""

        out_rows.append({
            "seqid": seqid,
            "d__Domain": f"d__{VIRUS_DOMAIN}",
            "r__Realm": f"r__{clean_token(parsed['realm'])}",
            "k__Kingdom": f"k__{clean_token(parsed['kingdom'])}",
            "p__Phylum": f"p__{clean_token(parsed['phylum'])}",
            "c__Class": f"c__{clean_token(parsed['class'])}",
            "o__Order": f"o__{clean_token(parsed['order'])}",
            "f__Family": f"f__{clean_token(parsed['family'])}",
            "g__Genus": f"g__{clean_token(parsed['genus'])}",
            "s__Species": f"s__{clean_token(parsed['species'])}",
            "vitap_pi": vitap_pi,
            "vitap_confidence_level": vitap_conf,
        })

    out_df = pd.DataFrame(out_rows, columns=COLS_OUT)

    # Always write a valid file (header-only is allowed)
    out_df.to_csv(args.out, index=False)

    print("✅ Wrote:", args.out)
    print("📊 Stats:", stats)


if __name__ == "__main__":
    main()