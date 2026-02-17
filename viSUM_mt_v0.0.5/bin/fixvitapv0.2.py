import pandas as pd
import subprocess
import argparse
import os
from tqdm.contrib.concurrent import process_map

try:
    import taxopy
    TAXOPY_AVAILABLE = True
except ImportError:
    TAXOPY_AVAILABLE = False

parser = argparse.ArgumentParser(description="Fix VITAP classifications for vSUM format using ICTV + fallback.")
parser.add_argument("--vitap", required=True)
parser.add_argument("--ictv", required=True)
parser.add_argument("--out", default="vitap_fixed.csv")
parser.add_argument("--threads", type=int, default=4)
parser.add_argument("--fallback", choices=["taxonkit", "taxopy"], default="taxopy")
parser.add_argument("--taxdump_dir", help="Path to taxdump (required for taxopy)")
args = parser.parse_args()

ictv_df = pd.read_csv(args.ictv)
ictv_df.columns = [col.lower() for col in ictv_df.columns]
ictv_df = ictv_df.fillna("")

ictv_lookup = {}
for _, row in ictv_df.iterrows():
    for rank in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]:
        name = row.get(rank, "").strip().lower()
        if name:
            ictv_lookup.setdefault(name, {}).update({k: v for k, v in row.items()})

valid_names_by_rank = {rank: set(ictv_df[rank].dropna().unique().tolist()) for rank in ictv_df.columns}

taxdb = None
if args.fallback == "taxopy":
    if not TAXOPY_AVAILABLE:
        raise ImportError("taxopy not installed.")
    if not args.taxdump_dir:
        raise ValueError("--taxdump_dir is required for taxopy fallback")
    taxdb = taxopy.TaxDb(
        nodes_dmp=os.path.join(args.taxdump_dir, "nodes.dmp"),
        names_dmp=os.path.join(args.taxdump_dir, "names.dmp")
    )

prefix_map = {
    "domain": "d__", "realm": "r__", "kingdom": "k__", "phylum": "p__",
    "class": "c__", "order": "o__", "family": "f__", "genus": "g__", "species": "s__"
}
rank_order = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def clean(x):
    x = str(x).strip()
    if x.lower() in ["", "-", "na", "nan"]:
        return "unclassified"
    if "]_" in x:
        x = x.split("]_")[-1]  # Remove bracketed prefix like [Realm]_...
    return x


def fill_missing_ranks(lineage_dict):
    last_valid = None
    for rank in rank_order:
        val = lineage_dict.get(rank, "unclassified")
        if val.lower() == "unclassified":
            if last_valid and last_valid not in valid_names_by_rank.get(rank, set()):
                lineage_dict[rank] = f"{last_valid}_U"
        elif val.endswith("_U") and (val[:-2] not in valid_names_by_rank.get(rank, set())):
            if last_valid:
                lineage_dict[rank] = f"{last_valid}_U"
        else:
            last_valid = val
    return lineage_dict


def fill_from_name(name):
    name = clean(name).lower()
    if name in ictv_lookup:
        lineage = ictv_lookup[name]
        lineage = {rank: clean(lineage.get(rank, "unclassified")) for rank in rank_order}
        return fill_missing_ranks(lineage)
    if args.fallback == "taxopy":
        try:
            tax = taxopy.Taxon(name, taxdb)
            lineage = {}
            for tid in tax.lineage:
                r = taxopy.utils.get_rank_by_taxid(tid, taxdb)
                n = taxopy.utils.get_name_by_taxid(tid, taxdb)
                if r in rank_order:
                    lineage[r] = n
            return fill_missing_ranks({rank: clean(lineage.get(rank, "unclassified")) for rank in rank_order})
        except Exception:
            pass
    elif args.fallback == "taxonkit":
        try:
            tid = subprocess.run(["taxonkit", "name2taxid"], input=name, capture_output=True, text=True).stdout.strip().split("\t")[-1]
            if tid.isdigit():
                lineage = subprocess.run(["taxonkit", "lineage"], input=tid, capture_output=True, text=True).stdout.strip().split("\t")[-1]
                parts = [clean(p) for p in lineage.split(";")]
                while len(parts) < 8:
                    parts.append("unclassified")
                return fill_missing_ranks(dict(zip(rank_order, parts)))
        except Exception:
            pass
    name = name.replace(" ", "_")
    return fill_missing_ranks({rank: name + "_U" for rank in rank_order})


def process_row(x):
    seqid = str(x["Genome_ID"]).strip()
    raw = str(x["lineage"])
    parts = [clean(p) for p in raw.split(";")[::-1]]
    lineage = dict(zip(rank_order, parts))
    unique_vals = set(parts)

    if len(unique_vals) <= 3 or all(p == parts[0] for p in parts):
        fallback_lineage = fill_from_name(parts[0])
        fallback_lineage["species"] = "unclassified"
        lineage = fallback_lineage
    else:
        for rank in ["realm", "kingdom", "phylum", "class"]:
            if lineage.get("species", "").lower() == lineage.get(rank, "").lower():
                lineage["species"] = "unclassified"

    lineage["domain"] = "Viruses"
    lineage = fill_missing_ranks(lineage)
    return [seqid] + [f"{prefix_map[r]}{lineage[r]}" for r in ["domain"] + rank_order]


# Main execution

df = pd.read_csv(args.vitap, sep="\t")
results = process_map(process_row, [row for _, row in df.iterrows()], max_workers=args.threads)
columns = ["seqid"] + [prefix_map[r] for r in ["domain"] + rank_order]
pd.DataFrame(results, columns=columns).to_csv(args.out, index=False)
print(f"\u2714 Written: {args.out}")

