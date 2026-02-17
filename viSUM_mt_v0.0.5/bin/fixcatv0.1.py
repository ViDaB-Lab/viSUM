#!/usr/bin/env python3

import csv
import sys
import pandas as pd

def load_ictv_lookup(ictv_file):
    ictv_df = pd.read_csv(ictv_file)
    ictv_df.columns = [col.strip().lower() for col in ictv_df.columns]

    # Build valid sets for each rank to prevent cross-rank assignment
    valid_values = {rank: set(ictv_df[rank].dropna().unique()) for rank in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']}

    # Realm lookup
    realm_lookup = {}
    for rank in valid_values:
        realm_lookup[rank] = dict(ictv_df[[rank, 'realm']].dropna().drop_duplicates().values.tolist())

    # Backfilling lookup: lower -> higher
    full_lookup = {rank: {} for rank in valid_values}
    for higher_rank in valid_values:
        for lower_rank in valid_values:
            if higher_rank == lower_rank:
                continue
            sub_df = ictv_df[[lower_rank, higher_rank]].dropna().drop_duplicates()
            full_lookup[lower_rank].update(dict(sub_df.values.tolist()))

    return full_lookup, realm_lookup, valid_values

def infer_tax_from_lower(tax_dict, desired_rank, full_lookup, valid_values):
    for lower_rank in ['species', 'genus', 'family', 'order', 'class', 'phylum']:
        val = tax_dict.get(lower_rank)
        if val and val in full_lookup[lower_rank]:
            candidate = full_lookup[lower_rank][val]
            if candidate in valid_values[desired_rank]:
                return candidate
    return "unclassified"

def infer_realm(tax_dict, realm_lookup):
    for rank in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']:
        val = tax_dict.get(rank)
        if val and val in realm_lookup[rank]:
            return realm_lookup[rank][val]
    return "unclassified"

def backfill_taxonomy(tax_dict, full_lookup, valid_values):
    tax_dict = tax_dict.copy()
    ranks = list(valid_values.keys())
    for i in range(len(ranks) - 1, 0, -1):
        lower, higher = ranks[i], ranks[i - 1]
        if higher not in tax_dict or tax_dict[higher].lower() in ['unclassified', 'na', '']:
            candidate = tax_dict.get(lower)
            if candidate and candidate in full_lookup[lower]:
                proposed = full_lookup[lower][candidate]
                if proposed in valid_values[higher]:
                    tax_dict[higher] = proposed
    return tax_dict

def clean_name(name):
    name = name.strip()
    if name.lower() in ['na', 'no support', '']:
        return None
    return name

def main(cat_file, ictv_file, output_file):
    full_lookup, realm_lookup, valid_values = load_ictv_lookup(ictv_file)
    rank_names = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    output_rows = []

    with open(cat_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 6:
                continue

            contig_id = row[0].strip()
            tax_dict = {}

            for i, val in enumerate(row[5:5 + len(rank_names)]):
                rank = rank_names[i]
                name = val.split(":")[0].strip()
                cleaned = clean_name(name)
                if cleaned:
                    tax_dict[rank] = cleaned

            # Backfill and restrict to proper rank values
            tax_dict = backfill_taxonomy(tax_dict, full_lookup, valid_values)

            d = f"d__{tax_dict.get('domain', 'unclassified')}"
            r = f"r__{infer_realm(tax_dict, realm_lookup)}"
            k = f"k__{tax_dict.get('kingdom', infer_tax_from_lower(tax_dict, 'kingdom', full_lookup, valid_values))}"
            p = f"p__{tax_dict.get('phylum', 'unclassified')}"
            c = f"c__{tax_dict.get('class', 'unclassified')}"
            o = f"o__{tax_dict.get('order', 'unclassified')}"
            f_ = f"f__{tax_dict.get('family', 'unclassified')}"
            g = f"g__{tax_dict.get('genus', 'unclassified')}"
            s = f"s__{tax_dict.get('species', 'unclassified')}"

            output_rows.append([contig_id, d, r, k, p, c, o, f_, g, s])

    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['contig_id', 'd__', 'r__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'])
        writer.writerows(output_rows)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python fixcat_vsum_ictvfill.py <cat_virus_subset.csv> <ICTV_masterlist.csv> <output.csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])

