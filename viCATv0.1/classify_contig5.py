#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
import argparse
from collections import defaultdict

# rank ? prefix for unclassified
RANK_PREFIX = {
    1: 'r__', 2: 'k__', 3: 'p__', 4: 'c__',
    5: 'o__', 6: 'f__', 7: 'g__', 8: 's__',
}

def classify_contig(hits_tsv, f, normalize):
    # 1) read hits, group by ORF
    hits_per_orf = defaultdict(list)
    with open(hits_tsv, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for row in reader:
            if len(row) < 6:
                continue
            orf      = row[0]
            bitscore = float(row[5])
            lineage  = row[-1].split(';')
            hits_per_orf[orf].append((bitscore, lineage))

    # 2) aggregate into bsums
    bsums    = defaultdict(float)
    total_bs = 0.0
    max_rank = 0

    if normalize:
        # each ORF casts 1 vote split evenly
        for orf, orf_hits in hits_per_orf.items():
            weight = 1.0 / len(orf_hits)
            total_bs += 1.0
            for bitscore, lineage in orf_hits:
                max_rank = max(max_rank, len(lineage))
                for i, tax in enumerate(lineage, start=1):
                    if tax and tax.lower() != 'unclassified':
                        bsums[(i, tax)] += weight
    else:
        # classic: each ORF contributes its best hit bitscore
        for orf, orf_hits in hits_per_orf.items():
            bitscore, lineage = orf_hits[0]
            total_bs += bitscore
            max_rank = max(max_rank, len(lineage))
            for i, tax in enumerate(lineage, start=1):
                if tax and tax.lower() != 'unclassified':
                    bsums[(i, tax)] += bitscore

    # 3) decide LCA: skip empty ranks, stop when support < f
    classification    = {}
    classification_bs = {}
    best_depth = 0

    for rank in range(1, max_rank + 1):
        # collect sums for this rank
        sums = {tax: bs for (r, tax), bs in bsums.items() if r == rank}
        if not sums:
            # no contributions here ? skip
            continue
        best_tax, best_bs = max(sums.items(), key=lambda x: x[1])
        if best_bs / total_bs >= f:
            best_depth = rank
            classification[rank]    = best_tax
            classification_bs[rank] = best_bs
        else:
            # had contributions but not enough support ? stop
            break

    # 4) build full LCA path, filling unclassified where no call
    lca = []
    for r in range(1, max_rank + 1):
        lca.append(classification.get(r, RANK_PREFIX.get(r, '') + 'unclassified'))

    # 5) get support bitscore/votes for deepest accepted rank
    support_bs = classification_bs.get(best_depth, 0.0)

    print(f"{';'.join(lca)}\t{support_bs}\t{total_bs}")


def main():
    p = argparse.ArgumentParser(
        description="Contig LCA classifier (CAT style), handles leading empty ranks"
    )
    p.add_argument("hits_tsv",
                   help="TSV of top hit(s) for all ORFs of one contig")
    p.add_argument("f", type=float,
                   help="fraction cutoff for LCA (e.g. 0.5)")
    p.add_argument("-n", "--normalize", action="store_true",
                   help="if set, split each ORFs vote equally across its hits")
    args = p.parse_args()

    classify_contig(args.hits_tsv, args.f, args.normalize)


if __name__ == "__main__":
    main()
