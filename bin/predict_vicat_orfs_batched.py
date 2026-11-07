#!/usr/bin/env python3
"""Run the existing viCAT GV/RV meta predictors in bounded, ordered batches."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def batches(records, max_bases: int, max_records: int):
    """Never split a record; an oversized record is processed alone."""
    batch, bases = [], 0
    for record in records:
        size = len(record.seq)
        if batch and (bases + size > max_bases or len(batch) >= max_records):
            yield batch
            batch, bases = [], 0
        batch.append(record)
        bases += size
        if bases >= max_bases or len(batch) >= max_records:
            yield batch
            batch, bases = [], 0
    if batch:
        yield batch


def predict(args) -> None:
    # Same parser, finders, parameters and writers as the pinned meta-mode CLI.
    # Only scheduling differs: each map receives one bounded batch.
    from pyrodigal.cli import parse, zopen
    if args.caller == "gv":
        from pyrodigal_gv import ViralGeneFinder
    else:
        from pyrodigal_rv import ViralGeneFinder

    if len({args.input.resolve(), args.proteins.resolve(), args.gff.resolve()}) != 3:
        raise ValueError("Input FASTA and the two output files must be distinct")
    finder = ViralGeneFinder(
        meta=True, closed=False, mask=False, training_info=None,
        min_gene=90, min_edge_gene=60, max_overlap=60,
    )

    def call(record):
        return record.id, finder.find_genes(record.seq)

    records_done = genes_done = bases_done = batch_count = 0
    largest_batch = largest_record = largest_batch_records = 0
    with zopen(args.input) as source, args.proteins.open("w") as proteins, \
            args.gff.open("w") as gff, ThreadPoolExecutor(max_workers=args.threads) as pool:
        for batch in batches(parse(source), args.batch_bases, args.batch_records):
            batch_count += 1
            batch_bases = sum(len(record.seq) for record in batch)
            largest_batch = max(largest_batch, batch_bases)
            largest_record = max(largest_record, max(len(record.seq) for record in batch))
            largest_batch_records = max(largest_batch_records, len(batch))
            for sequence_id, predictions in pool.map(call, batch):
                predictions.write_gff(gff, sequence_id)
                predictions.write_translations(proteins, sequence_id, include_stop=True)
                genes_done += len(predictions)
            records_done += len(batch)
            bases_done += batch_bases
            if batch_count == 1 or batch_count % 100 == 0:
                proteins.flush()
                gff.flush()
                print(f"viCAT {args.caller}: batches={batch_count} sequences={records_done} "
                      f"bases={bases_done} genes={genes_done}", flush=True)
    print(f"viCAT {args.caller} complete: threads={args.threads} batches={batch_count} "
          f"sequences={records_done} bases={bases_done} genes={genes_done} "
          f"max_batch_bases={largest_batch} max_batch_records={largest_batch_records} "
          f"max_sequence_bases={largest_record}", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caller", choices=("gv", "rv"), required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--threads", type=positive_int, required=True)
    parser.add_argument("--batch-bases", type=positive_int, default=1_000_000)
    parser.add_argument("--batch-records", type=positive_int, default=256)
    predict(parser.parse_args(argv))


if __name__ == "__main__":
    main()
