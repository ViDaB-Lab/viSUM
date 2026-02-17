import pandas as pd
import subprocess
import argparse
import os
import time
from tqdm.contrib.concurrent import process_map

try:
    import taxopy
    TAXOPY_AVAILABLE = True
except ImportError:
    TAXOPY_AVAILABLE = False

parser = argparse.ArgumentParser(description="Fix Cenote-Taker3 classifications using ICTV + fallback (taxonkit or taxopy).")
parser.add_argument("--ct3", required=True, help="Path to Cenote-Taker3 summary TSV file")
parser.add_argument("--ictv", required=True, help="Path to ICTV master species CSV")
parser.add_argument("--out", default="ct3_fixed.csv", help="Output path for fixed CT3 classifications")
parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
parser.add_argument("--fallback", choices=["taxonkit", "taxopy"], default="taxonkit", help="Fallback method for unmatched taxa")
parser.add_argument("--taxdump_dir", help="Path to taxdump directory (only for taxopy)")
args = parser.parse_args()

ictv_df = pd.read_csv(args.ictv)
ictv_df.columns = [col.lower() for col in ictv_df.columns]
ictv_df = ictv_df.fillna("")

ictv_lookup = []
for _, row in ictv_df.iterrows():
    lineage = {
        "domain": "Viruses",
        "realm": row.get("realm", "").strip() or "unclassified",
        "kingdom": row.get("kingdom", "").strip() or "unclassified",
        "phylum": row.get("phylum", "").strip() or "unclassified",
        "class": row.get("class", "").strip() or "unclassified",
        "order": row.get("order", "").strip() or "unclassified",
        "family": row.get("family", "").strip() or "unclassified",
        "genus": row.get("genus", "").strip() or "unclassified",
        "species": row.get("species", "").strip() or "unclassified"
    }
    ictv_lookup.append(lineage)

if args.fallback == "taxopy":
    if not TAXOPY_AVAILABLE:
        raise ImportError("taxopy is not installed. Please install with pip install taxopy")
    if not args.taxdump_dir:
        raise ValueError("--taxdump_dir is required for taxopy fallback")
    taxdb = taxopy.TaxDb(
        nodes_dmp=os.path.join(args.taxdump_dir, "nodes.dmp"),
        names_dmp=os.path.join(args.taxdump_dir, "names.dmp")
    )
else:
    taxdb = None

prefix_map = {
    "domain": "d__", "realm": "r__", "kingdom": "k__", "phylum": "p__",
    "class": "c__", "order": "o__", "family": "f__", "genus": "g__", "species": "s__"
}
rank_keys = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

def clean(name):
    return "unclassified" if str(name).strip() in ["", "-"] else str(name).strip()

def get_fallback_lineage(x):
    seqid, name = x
    name = clean(name)
    try:
        if args.fallback == "taxonkit":
            proc = subprocess.run(["taxonkit", "name2taxid"], input=name, text=True, capture_output=True)
            taxid = proc.stdout.strip().split("\t")[-1]
            if not taxid.isdigit(): return default_lineage(seqid)
            proc = subprocess.run(["taxonkit", "lineage"], input=taxid, text=True, capture_output=True)
            parts = proc.stdout.strip().split("\t")[-1].split(";")
        elif args.fallback == "taxopy":
            tax = taxopy.Taxon(name, taxdb)
            lineage = tax.lineage
            parts = [taxopy.utils.get_name_by_taxid(t, taxdb) for t in lineage if taxopy.utils.get_rank_by_taxid(t, taxdb) in rank_keys]
        else:
            return default_lineage(seqid)
        while len(parts) < 8:
            parts.append("unclassified")
        lineage_dict = dict(zip(rank_keys, parts))
        lineage_dict["domain"] = "Viruses"
        return [seqid] + [f"{prefix_map[r]}{clean(lineage_dict.get(r, 'unclassified'))}" for r in prefix_map]
    except:
        return default_lineage(seqid)

def default_lineage(seqid):
    return [seqid] + [f"{prefix_map[r]}{'Viruses' if r=='domain' else 'unclassified'}" for r in prefix_map]

def best_ictv_match(tax_dict):
    for lineage in ictv_lookup:
        match = True
        for rank in reversed(rank_keys):
            prefix = rank[0] + "_"
            if prefix in tax_dict:
                if clean(lineage.get(rank, "")) != tax_dict[prefix]:
                    match = False
                    break
                else:
                    return lineage
    return None

ct3_rows = []
ct3_fallback = []
ct3_df = pd.read_csv(args.ct3, sep="\t")

for _, row in ct3_df.iterrows():
    try:
        seqid = str(row.iloc[1]).strip()
        class_str = row.iloc[12] if len(row) > 12 else ""
        if pd.isna(class_str) or class_str.strip() == "" or class_str.strip().lower() == "unclassified virus":
            ct3_fallback.append((seqid, "unclassified"))
            continue
        parts = class_str.replace("-_", "").split(";")
        tax = {p[:2]: clean(p[2:]) for p in parts if len(p) >= 3 and p[1] == '_'}
        lineage = best_ictv_match(tax)
        if lineage:
            row_data = [seqid] + [f"{prefix_map[r]}{clean(lineage[r])}" for r in prefix_map]
            ct3_rows.append(row_data)
        else:
            # Fallback to the deepest known taxon in CT3 string
            for rank_prefix in reversed(["s_", "g_", "f_", "o_", "c_", "p_", "k_", "r_"]):
                if rank_prefix in tax:
                    ct3_fallback.append((seqid, tax[rank_prefix]))
                    break
            else:
                ct3_fallback.append((seqid, "unclassified"))
    except Exception:
        ct3_fallback.append((seqid, "unclassified"))

if ct3_fallback:
    print(f"\n🌀 Running fallback on {len(ct3_fallback)} CT3 entries...")
    fallback_rows = process_map(get_fallback_lineage, ct3_fallback, max_workers=args.threads, chunksize=10, desc="CT3 fallback")
    ct3_rows.extend(fallback_rows)

columns = ["seqid"] + list(prefix_map.values())
pd.DataFrame(ct3_rows, columns=columns).to_csv(args.out, index=False)
print(f"\n✔ CT3 classification written to: {args.out}")

