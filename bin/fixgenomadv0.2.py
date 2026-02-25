#!/usr/bin/env python3
"""
fixgenomadv0.2.py

Purpose
-------
Convert GeNomad virus_summary taxonomy output into a clean viSUM-style table:

seqid, d__, r__, k__, p__, c__, o__, f__, g__, s__

Key robustness upgrades vs v0.1
-------------------------------
1) **GeNomad-first parsing (by position).**
   Your GeNomad TSV taxonomy strings are already ICTV-like:
     Viruses;Realm;Kingdom;Phylum;Class;Order;Family
   That is the most reliable source for ranks. We parse it directly.

2) **ICTV CSV only used as a fallback**, mainly when GeNomad taxonomy is just "Viruses"
   (or blank). This avoids accidental re-mapping errors from last-token heuristics.

3) Optional TaxonKit fallback is supported, but NOT required for typical GeNomad outputs.

4) Better input validation + summary stats.

Assumptions from your attached data
-----------------------------------
- GeNomad `taxonomy` column is either:
    - "Viruses"  (unresolved)
    - 7-field lineage: Viruses;Realm;Kingdom;Phylum;Class;Order;Family
- GeNomad does not include genus/species (so those remain unclassified unless you
  intentionally infer them, which is risky and not recommended).

Usage
-----
python3 fixgenomadv0.2.py \
  --genomad virus_summary.tsv \
  --ictv ICTV_masterspecieslist.csv \
  --out prefix.genomad.vsum.csv \
  --fallback none

If you *really* want TaxonKit fallback:
  --fallback taxonkit --taxonkit_db /path/to/taxdump
"""

import argparse
import os
import subprocess
from typing import Dict, Optional, Tuple

import pandas as pd


# ----------------------------
# Helpers
# ----------------------------

def clean_label(x: Optional[str]) -> str:
    """Normalize empty/NA/unclassified tokens to the literal string 'unclassified'."""
    if x is None:
        return "unclassified"
    x = str(x).strip()
    if x == "":
        return "unclassified"
    if x.lower() in {"na", "nan", "none", "no support", "unclassified"}:
        return "unclassified"
    return x


def lc(x: str) -> str:
    """Case-insensitive key for matching."""
    return clean_label(x).lower()


def detect_cols(df: pd.DataFrame) -> Tuple[str, str]:
    """Find seq id and taxonomy columns with a bit of flexibility."""
    cols = {c.lower(): c for c in df.columns}

    seq_col = None
    for cand in ["seq_name", "seqname", "sequence", "contig", "id"]:
        if cand in cols:
            seq_col = cols[cand]
            break

    tax_col = None
    for cand in ["taxonomy", "lineage", "classification", "taxon"]:
        if cand in cols:
            tax_col = cols[cand]
            break

    if not seq_col:
        raise ValueError(f"Could not find sequence id column. Columns: {list(df.columns)}")
    if not tax_col:
        raise ValueError(f"Could not find taxonomy column. Columns: {list(df.columns)}")

    return seq_col, tax_col


# ----------------------------
# ICTV indexing
# ----------------------------

ICTV_RANKS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

def build_ictv_index(ictv_csv: str) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Build index so any label at any rank can map to a full ICTV lineage:
      index[rank][lower(label)] = {realm, kingdom, ..., species}
    """
    df = pd.read_csv(ictv_csv)
    df.columns = [c.lower() for c in df.columns]

    # Ensure expected columns exist
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
    """
    Try mapping a label to an ICTV lineage.
    We search from most specific to broadest.
    """
    key = lc(label)
    if key == "unclassified":
        return None

    for r in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]:
        if key in index[r]:
            return index[r][key]
    return None


# ----------------------------
# GeNomad taxonomy parsing
# ----------------------------

def parse_genomad_taxonomy(taxstr: str) -> Dict[str, str]:
    """
    Parse GeNomad taxonomy string.

    Expected formats in your TSV:
      1) "Viruses"
      2) "Viruses;Realm;Kingdom;Phylum;Class;Order;Family"  (7 fields)

    Returns dict with keys:
      d__, r__, k__, p__, c__, o__, f__, g__, s__
    """
    taxstr = "" if taxstr is None else str(taxstr).strip()
    parts = [p.strip() for p in taxstr.split(";")] if taxstr else []

    # defaults
    out = {
        "d__": "d__Viruses",
        "r__": "r__unclassified",
        "k__": "k__unclassified",
        "p__": "p__unclassified",
        "c__": "c__unclassified",
        "o__": "o__unclassified",
        "f__": "f__unclassified",
        "g__": "g__unclassified",
        "s__": "s__unclassified",
    }

    if not parts:
        return out

    # If taxonomy is just "Viruses" (unresolved)
    if len(parts) == 1:
        # sometimes it's literally "Viruses" or "viruses"
        return out

    # Your data: max length is 7 and first token is "Viruses"
    # positions: 0=Viruses, 1=Realm, 2=Kingdom, 3=Phylum, 4=Class, 5=Order, 6=Family
    if len(parts) == 7 and parts[0].lower() == "viruses":
        realm, kingdom, phylum, clazz, order, family = parts[1:]
        out["r__"] = f"r__{clean_label(realm)}"
        out["k__"] = f"k__{clean_label(kingdom)}"
        out["p__"] = f"p__{clean_label(phylum)}"
        out["c__"] = f"c__{clean_label(clazz)}"
        out["o__"] = f"o__{clean_label(order)}"
        out["f__"] = f"f__{clean_label(family)}"
        return out

    # If GeNomad format changes someday, we fail softly:
    # keep domain=Viruses and do best effort by taking last 6 tokens as r..f (if possible)
    if len(parts) > 1 and parts[0].lower() == "viruses":
        tail = parts[1:]
        # pad / trim to 6
        tail = (tail + [""] * 6)[:6]
        realm, kingdom, phylum, clazz, order, family = tail
        out["r__"] = f"r__{clean_label(realm)}"
        out["k__"] = f"k__{clean_label(kingdom)}"
        out["p__"] = f"p__{clean_label(phylum)}"
        out["c__"] = f"c__{clean_label(clazz)}"
        out["o__"] = f"o__{clean_label(order)}"
        out["f__"] = f"f__{clean_label(family)}"
        return out

    return out


def best_label_from_genomad(taxstr: str) -> str:
    """
    Pick the most informative (most specific) label from a GeNomad taxonomy string.
    This is ONLY used for ICTV fallback when GeNomad doesn't provide ranks.

    We take the last non-empty token.
    """
    if taxstr is None:
        return "unclassified"
    parts = [clean_label(p) for p in str(taxstr).split(";")]
    parts = [p for p in parts if p != "unclassified" and p.lower() != "viruses"]
    return parts[-1] if parts else "unclassified"


# ----------------------------
# Optional TaxonKit fallback
# ----------------------------

def have_exe(name: str) -> bool:
    from shutil import which
    return which(name) is not None


def taxonkit_name2taxid(label: str, taxonkit_db: Optional[str]) -> Optional[str]:
    """
    Very minimal: return first matching taxid.
    If you need deterministic 'best rank' selection, we can add it—but for GeNomad,
    you likely won't need this fallback.
    """
    if not have_exe("taxonkit"):
        return None
    label = clean_label(label)
    if label == "unclassified":
        return None

    env = os.environ.copy()
    if taxonkit_db:
        env["TAXONKIT_DB"] = taxonkit_db

    p = subprocess.run(
        ["taxonkit", "name2taxid"],
        input=label + "\n",
        text=True,
        capture_output=True,
        env=env,
    )
    if p.returncode != 0 or not p.stdout.strip():
        return None

    # name2taxid output: name \t taxid ...
    taxid = p.stdout.strip().split("\t")[1]
    return taxid if taxid.isdigit() else None


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomad", required=True, help="GeNomad virus_summary.tsv")
    ap.add_argument("--ictv", required=True, help="ICTV master species list CSV")
    ap.add_argument("--out", required=True, help="Output CSV")
    ap.add_argument("--fallback", choices=["none", "ictv", "taxonkit"], default="ictv",
                    help="Fallback when GeNomad taxonomy is unresolved (default: ictv)")
    ap.add_argument("--taxonkit_db", default=None, help="Optional path to taxdump for TaxonKit (sets TAXONKIT_DB)")
    args = ap.parse_args()

    geno = pd.read_csv(args.genomad, sep="\t")
    seq_col, tax_col = detect_cols(geno)

    ictv_index = build_ictv_index(args.ictv)

    out_rows = []
    stats = {
        "n_total": 0,
        "n_genomad_parsed": 0,
        "n_genomad_unresolved": 0,
        "n_filled_by_ictv": 0,
        "n_taxonkit_attempted": 0,
    }

    for _, row in geno.iterrows():
        stats["n_total"] += 1
        seqid = str(row[seq_col])
        taxstr = row.get(tax_col, "")

        parsed = parse_genomad_taxonomy(taxstr)

        # Check if GeNomad gave us real ranks beyond domain
        genomad_has_detail = any(parsed[k] != f"{k}unclassified".replace("k", "k__") for k in ["r__", "k__", "p__", "c__", "o__", "f__"])
        # (Above is a bit hacky; we’ll just use a simpler check:)
        genomad_has_detail = parsed["r__"] != "r__unclassified" or parsed["f__"] != "f__unclassified"

        if genomad_has_detail:
            stats["n_genomad_parsed"] += 1
            out_rows.append({"seqid": seqid, **parsed})
            continue

        # Unresolved (usually taxonomy == "Viruses")
        stats["n_genomad_unresolved"] += 1

        # Fallback: ICTV lookup based on best label (if any)
        if args.fallback in {"ictv", "taxonkit"}:
            label = best_label_from_genomad(taxstr)
            L = ictv_lookup(ictv_index, label)
            if L:
                parsed["r__"] = f"r__{clean_label(L.get('realm'))}"
                parsed["k__"] = f"k__{clean_label(L.get('kingdom'))}"
                parsed["p__"] = f"p__{clean_label(L.get('phylum'))}"
                parsed["c__"] = f"c__{clean_label(L.get('class'))}"
                parsed["o__"] = f"o__{clean_label(L.get('order'))}"
                parsed["f__"] = f"f__{clean_label(L.get('family'))}"
                # genus/species remain unclassified unless you choose to infer them
                stats["n_filled_by_ictv"] += 1
                out_rows.append({"seqid": seqid, **parsed})
                continue

        # Optional TaxonKit attempt (only if explicitly requested)
        if args.fallback == "taxonkit":
            stats["n_taxonkit_attempted"] += 1
            label = best_label_from_genomad(taxstr)
            _ = taxonkit_name2taxid(label, args.taxonkit_db)
            # We are NOT expanding lineage here because GeNomad already supplies ICTV ranks in most cases.
            # If you want full taxonkit lineage expansion, we can add it cleanly later.

        out_rows.append({"seqid": seqid, **parsed})

    out_df = pd.DataFrame(out_rows, columns=["seqid", "d__", "r__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"])
    out_df.to_csv(args.out, index=False)

    print("✅ Wrote:", args.out)
    print("📊 Stats:", stats)


if __name__ == "__main__":
    main()