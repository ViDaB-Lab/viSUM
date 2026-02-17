import pandas as pd
import subprocess
import argparse
import os
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm

try:
    import taxopy
    TAXOPY_AVAILABLE = True
except ImportError:
    TAXOPY_AVAILABLE = False

parser = argparse.ArgumentParser(description="Fix GenoMad classifications using ICTV + fallback")
parser.add_argument("--genomad", required=True, help="Path to GenoMad virus summary TSV file")
parser.add_argument("--ictv", required=True, help="Path to ICTV master species CSV")
parser.add_argument("--out", default="genomad_fixed.csv", help="Output file for standardized classifications")
parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads")
parser.add_argument("--fallback", choices=["taxonkit", "taxopy"], default="taxonkit", help="Fallback method")
parser.add_argument("--taxdump_dir", help="Path to taxdump dir for taxopy")
args = parser.parse_args()

print("📥 Loading ICTV taxonomy...")
ictv_df = pd.read_csv(args.ictv)
ictv_df.columns = [col.lower() for col in ictv_df.columns]

ictv_dict = {rank: {} for rank in ["species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"]}

for _, row in ictv_df.fillna("").iterrows():
    lineage = {rank: row.get(rank, "").strip() for rank in ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]}
    for rank in ictv_dict:
        key = lineage[rank].strip().lower()
        if key:
            ictv_dict[rank][key] = lineage

if args.fallback == "taxopy":
    if not TAXOPY_AVAILABLE:
        raise ImportError("Please install taxopy for fallback")
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
rank_keys = ["realm", "kingdom", "phylum", "class", "order", "family"]

print("📄 Reading GenoMad summary file...")
geno_df = pd.read_csv(args.genomad, sep="\t")

def clean(x):
    x = str(x).strip()
    return "unclassified" if x.lower() in ["", "na", "nan", "no support", "unclassified"] else x

def fill_lineage_from_dict(label):
    label = clean(label).lower()
    for rank in reversed(rank_keys):
        if label in ictv_dict[rank]:
            lineage = ictv_dict[rank][label].copy()
            return ["d__Viruses"] + [prefix_map[k] + clean(lineage.get(k, "unclassified")) for k in rank_keys] + ["g__unclassified", "s__unclassified"]
    return None

def use_taxonomy_as_is(seqid, taxstr):
    parts = [clean(x) for x in taxstr.split(";")]
    lineage = ["d__Viruses"]
    for i, rank in enumerate(rank_keys):
        value = parts[i + 1] if len(parts) > i + 1 else "unclassified"
        lineage.append(prefix_map[rank] + clean(value))
    lineage += ["g__unclassified", "s__unclassified"]
    return [seqid] + lineage

def get_fallback_lineage(x):
    seqid, label, taxstr = x
    try:
        if args.fallback == "taxonkit":
            p = subprocess.run(["taxonkit", "name2taxid"], input=label, capture_output=True, text=True)
            taxid = p.stdout.strip().split("\t")[-1]
            if not taxid.isdigit():
                return use_taxonomy_as_is(seqid, taxstr)
            p = subprocess.run(["taxonkit", "lineage"], input=taxid, capture_output=True, text=True)
            parts = p.stdout.strip().split("\t")[-1].split(";")
        elif args.fallback == "taxopy":
            tax = taxopy.Taxon(label, taxdb)
            lineage = tax.lineage
            parts = [taxopy.utils.get_name_by_taxid(t, taxdb) for t in lineage if taxopy.utils.get_rank_by_taxid(t, taxdb) in rank_keys]
        else:
            return use_taxonomy_as_is(seqid, taxstr)
        while len(parts) < 6:
            parts.append("unclassified")
        d = dict(zip(rank_keys, parts))
        return [seqid] + ["d__Viruses"] + [prefix_map[k] + clean(d.get(k, "unclassified")) for k in rank_keys] + ["g__unclassified", "s__unclassified"]
    except:
        return use_taxonomy_as_is(seqid, taxstr)

print("🔍 Processing taxonomy...")
out_rows = []
fallback_needed = []

for _, row in tqdm(geno_df.iterrows(), total=len(geno_df)):
    seqid = row["seq_name"]
    taxstr = str(row.get("taxonomy", ""))
    parts = [p for p in taxstr.split(";") if clean(p) != "unclassified"]
    label = parts[-1] if parts else "unclassified"

    lineage = fill_lineage_from_dict(label)
    if lineage:
        out_rows.append([seqid] + lineage)
    else:
        fallback_needed.append((seqid, label, taxstr))

if fallback_needed:
    print(f"🧪 Running fallback or 'as-is' formatting on {len(fallback_needed)} entries...")
    fallback_results = process_map(get_fallback_lineage, fallback_needed, max_workers=args.threads, chunksize=10, desc="Fallback")
    out_rows.extend(fallback_results)

columns = ["seqid", "d__"] + [prefix_map[k] for k in rank_keys] + ["g__", "s__"]
pd.DataFrame(out_rows, columns=columns).to_csv(args.out, index=False)
print(f"\n✅ GenoMad classifications written to: {args.out}")

