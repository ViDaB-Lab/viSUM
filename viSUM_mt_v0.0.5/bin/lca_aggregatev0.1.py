#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import sys

def load_ictv(path):
    legal_realm, legal_kingdom, legal_phylum = set(), set(), set()
    legal_class, legal_order, legal_family = set(), set(), set()
    legal_genus = set()
    orphan_kingdom, orphan_phylum, orphan_class = set(), set(), set()
    orphan_orders, orphan_families, orphan_genus = set(), set(), set()

    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            realm = row.get('Realm','').strip()
            if not realm or realm.lower() == 'realm':
                continue
            kingdom = row.get('Kingdom','').strip()
            phylum  = row.get('Phylum','').strip()
            cls     = row.get('Class','').strip()
            order   = row.get('Order','').strip()
            family  = row.get('Family','').strip()
            genus   = row.get('Genus','').strip()

            legal_realm.add(realm)
            if kingdom:
                legal_kingdom.add(f"{realm}|{kingdom}")
            else:
                orphan_kingdom.add(realm)

            if phylum:
                legal_phylum.add(f"{realm}|{kingdom}|{phylum}")
            else:
                orphan_phylum.add(f"{realm}|{kingdom}")

            if cls:
                legal_class.add(f"{realm}|{kingdom}|{phylum}|{cls}")
            else:
                orphan_class.add(f"{realm}|{kingdom}|{phylum}")

            if order:
                legal_order.add(f"{realm}|{kingdom}|{phylum}|{cls}|{order}")
            else:
                orphan_orders.add(f"{realm}|{kingdom}|{phylum}|{cls}")

            if family:
                legal_family.add(f"{realm}|{kingdom}|{phylum}|{cls}|{order}|{family}")
            else:
                orphan_families.add(f"{realm}|{kingdom}|{phylum}|{cls}|{order}")

            if genus:
                legal_genus.add(f"{realm}|{kingdom}|{phylum}|{cls}|{order}|{family}|{genus}")
            else:
                orphan_genus.add(f"{realm}|{kingdom}|{phylum}|{cls}|{order}|{family}")

    return (
        legal_realm, legal_kingdom, legal_phylum,
        legal_class, legal_order, legal_family, legal_genus,
        orphan_kingdom, orphan_phylum, orphan_class,
        orphan_orders, orphan_families, orphan_genus
    )

def choose_best(weight_map):
    best, bestw = None, 0.0
    for full, w in weight_map.items():
        val = full.split('__',1)[1] if '__' in full else full
        if val == "unclassified":
            continue
        if w > bestw:
            best, bestw = full, w
    if best is None:
        return "unclassified", 0.0
    return best, bestw

def main():
    p = argparse.ArgumentParser(description="Aggregate taxonomic calls (family+genus filtering only).")
    p.add_argument('--input',     required=True, help="CSV: prog,d__,r__,k__,p__,c__,o__,f__,g__,s__")
    p.add_argument('--ictv',      required=True, help="ICTV master species list (CSV)")
    p.add_argument('--threshold', type=float, default=0.5, help="min family support fraction")
    p.add_argument('--weights',   nargs='+', required=True, help="prog=weight tokens")
    args = p.parse_args()

    weight = { prog: float(val) for tok in args.weights for prog, val in [tok.split('=',1)] }

    (
        legal_realm, legal_kingdom, legal_phylum,
        legal_class, legal_order, legal_family, legal_genus,
        orphan_kingdom, orphan_phylum, orphan_class,
        orphan_orders, orphan_families, orphan_genus
    ) = load_ictv(args.ictv)

    sum_d = {}; sum_r = {}; sum_k = {}
    sum_p = {}; sum_c = {}; sum_o = {}
    sum_f = {}; sum_g = {}
    total_w = 0.0

    with open(args.input) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            prog, *vals = row
            if len(vals) < 9:
                vals += ['unclassified'] * (9 - len(vals))
            d, r, k, p_, c, o, fam, g, _ = vals[:9]

            if not d or d.strip() in ["d__", ""]:
                d = "d__Viruses"

            w = weight.get(prog, 1.0)
            total_w += w
            for m, key in (
                (sum_d, d), (sum_r, r), (sum_k, k),
                (sum_p, p_), (sum_c, c), (sum_o, o),
                (sum_f, fam), (sum_g, g)
            ):
                m[key] = m.get(key, 0.0) + w

    d_best, dw = choose_best(sum_d)
    r_best, rw = choose_best(sum_r)
    k_best, kw = choose_best(sum_k)
    p_best, pw = choose_best(sum_p)
    c_best, cw = choose_best(sum_c)
    o_best, ow = choose_best(sum_o)

    strip = lambda x: x.split('__',1)[1] if '__' in x else x
    d_val = strip(d_best) or "Viruses"
    r_val = strip(r_best)
    k_val, p_val = strip(k_best), strip(p_best)
    c_val, o_val = strip(c_best), strip(o_best)

    valid_fams = {}
    invalid_w = 0.0
    for fam_name, w in sum_f.items():
        fam_val = strip(fam_name)
        lineage = f"{r_val}|{k_val}|{p_val}|{c_val}|{o_val}|{fam_val}"
        if lineage in legal_family or f"{r_val}|{k_val}|{p_val}|{c_val}|{o_val}" in orphan_families:
            valid_fams[fam_name] = w
        else:
            invalid_w += w
    valid_fams['unclassified'] = valid_fams.get('unclassified',0.0) + invalid_w

    f_best, fw = choose_best(valid_fams)
    g_best, gw = choose_best(sum_g)
    f_val, g_val = strip(f_best), strip(g_best)
    s_val = "unclassified"

    if f_val != "unclassified":
        fullg = f"{r_val}|{k_val}|{p_val}|{c_val}|{o_val}|{f_val}|{g_val}"
        if fullg not in legal_genus and f"{r_val}|{k_val}|{p_val}|{c_val}|{o_val}|{f_val}" not in orphan_genus:
            g_val = "unclassified"

    if total_w and (fw/total_w) < args.threshold:
        f_val = g_val = "unclassified"
    if gw != total_w:
        g_val = "unclassified"

    ranks = ["d__", "r__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    values = [d_val, r_val, k_val, p_val, c_val, o_val, f_val, g_val, s_val]
    final_taxonomy = [f"{r}{v}" for r, v in zip(ranks, values)]

    # NEW: Avoid _U fill if only d__Viruses is known
    informative = [v for v in values if v != "unclassified"]
    if len(informative) > 1:  # i.e., more than just d__Viruses
        last_valid = None
        for i in range(len(final_taxonomy)):
            val = final_taxonomy[i].split("__", 1)[1]
            if val == "unclassified":
                if last_valid:
                    base = last_valid.split("_U")[0]
                    final_taxonomy[i] = f"{ranks[i]}{base}_U"
            else:
                last_valid = val

    frac = lambda x: x/total_w if total_w else 0.0
    ftax = ",".join(final_taxonomy)
    frank = ";;".join([
        f"d__{frac(dw):.2f}", f"r__{frac(rw):.2f}",
        f"k__{frac(kw):.2f}", f"p__{frac(pw):.2f}",
        f"c__{frac(cw):.2f}", f"o__{frac(ow):.2f}",
        f"f__{frac(fw):.2f}", f"g__{frac(gw):.2f}",
    ])

    print(ftax)
    print(frank)

if __name__ == "__main__":
    main()

