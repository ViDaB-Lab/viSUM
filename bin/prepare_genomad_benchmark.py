#!/usr/bin/env python3
"""Audit author metadata and prepare fail-closed, reproducible viSUM inputs.

Standard library only. No network, classifier inference, or pipeline execution.
"""
from __future__ import annotations

import argparse
import collections
import contextlib
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sqlite3
import sys
import tarfile

TRUTH = 'prepared/truth_metadata.tsv'
CLUSTERS = 'source/reference_sequences/reference_clusters.tsv'
SUMMARY = 'prepared/setup_summary.json'
EXCLUDED = 'source/benchmark_data/excluded_from_benchmark.txt'
WEIGHTS = [f'author_weight_split_{i}' for i in range(6)]
REQUIRED = ['sequence_id', 'reference_cluster', 'source_class', *WEIGHTS, 'length_nt']


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


@contextlib.contextmanager
def source_files(source):
    """Read named members without extracting untrusted archive paths."""
    source = Path(source)
    if source.is_dir():
        yield lambda name: (source / name).open(encoding='utf-8', newline='')
    else:
        with tarfile.open(source, 'r:*') as tar:
            members = {}
            for member in tar:
                if member.name in members:
                    raise ValueError(f'Duplicate archive member: {member.name}')
                members[member.name] = member

            def opener(name):
                member = members.get(name)
                if member is None or not member.isfile():
                    raise ValueError(f'Missing regular metadata member: {name}')
                return io.TextIOWrapper(tar.extractfile(member), encoding='utf-8', newline='')
            yield opener


def table_writer(handle, fields, delimiter='\t'):
    writer = csv.writer(handle, delimiter=delimiter, lineterminator='\n')
    writer.writerow(fields)
    return writer


def parent_id(identifier):
    match = re.fullmatch(r'(.+)\|fragment_([1-9][0-9]*)', identifier)
    if not match or any(c.isspace() for c in identifier):
        raise ValueError(f'Invalid author fragment identifier: {identifier}')
    return match[1]


def load_routes(path):
    routes = {}
    if path is None:
        return routes
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = {'parent_sequence_id', 'input_type', 'evidence_source', 'evidence_record'}
        if (not fields.issubset(reader.fieldnames or []) or
                len(reader.fieldnames) != len(set(reader.fieldnames))):
            raise ValueError(f'Routing map requires {sorted(fields)}')
        for row in reader:
            key = row['parent_sequence_id']
            if not key or key in routes:
                raise ValueError(f'Blank/duplicate routing key: {key}')
            if row['input_type'] not in {'dna', 'rna'}:
                raise ValueError(f'Unresolved/invalid route: {key}')
            if not row['evidence_source'].strip() or not row['evidence_record'].strip():
                raise ValueError(f'Route lacks evidence provenance: {key}')
            routes[key] = (row['input_type'], row['evidence_source'], row['evidence_record'])
    return routes


def make_index(opener, db, routes):
    with opener(SUMMARY) as f:
        expected = json.load(f)
    with opener(EXCLUDED) as f:
        excluded = {line.strip() for line in f if line.strip()}
    db.executescript('''
        CREATE TABLE clusters (id TEXT PRIMARY KEY, split INTEGER);
        CREATE TABLE members (id TEXT PRIMARY KEY, cluster TEXT);
        CREATE TABLE truth (id TEXT PRIMARY KEY, parent TEXT, class TEXT,
            length INTEGER, cluster TEXT, split INTEGER, route TEXT,
            raw TEXT, seen INTEGER DEFAULT 0, batch TEXT);
    ''')
    with opener(CLUSTERS) as f:
        for row in csv.reader(f, delimiter='\t'):
            if len(row) != 4 or int(row[2]) not in range(1, 6):
                raise ValueError('Malformed author cluster record')
            ids = row[3].split(',')
            if int(row[1]) != len(ids) or any(not key for key in ids):
                raise ValueError(f'Cluster member-count mismatch: {row[0]}')
            db.execute('INSERT INTO clusters VALUES (?,?)', (row[0], int(row[2])))
            db.executemany('INSERT INTO members VALUES (?,?)', ((key, row[0]) for key in ids))
    db.commit()
    counts = collections.Counter()
    parents = {}
    bases = 0
    used_routes = set()
    with opener(TRUTH) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        if fields is None or len(fields) != len(set(fields)) or not set(REQUIRED).issubset(fields):
            raise ValueError('Missing/duplicate truth columns')
        for n, row in enumerate(reader, 1):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f'Malformed truth row {n}')
            key, cls = row['sequence_id'], row['source_class']
            parent = parent_id(key)
            if cls not in {'chromosome', 'plasmid', 'virus'} or key in excluded:
                raise ValueError(f'Invalid class or excluded fragment: {key}')
            length = int(row['length_nt'])
            if length <= 0:
                raise ValueError(f'Invalid length: {key}')
            for col in WEIGHTS:
                weight = float(row[col])
                if not math.isfinite(weight) or weight < 0:
                    raise ValueError(f'Invalid weight: {key}/{col}')
            if parent not in parents:
                record = db.execute('SELECT m.cluster,c.split FROM members m JOIN clusters c '
                                    'ON c.id=m.cluster WHERE m.id=?', (parent,)).fetchone()
                if record is None:
                    raise ValueError(f'Parent absent from reference clusters: {parent}')
                parents[parent] = (record[0], record[1], cls)
            cluster, split, parent_class = parents[parent]
            if cluster != row['reference_cluster'] or parent_class != cls:
                raise ValueError(f'Conflicting parent/cluster/class: {key}')
            route = 'dna' if cls != 'virus' else routes.get(parent, ('unresolved',))[0]
            if parent in routes:
                if cls != 'virus':
                    raise ValueError(f'Routing map must contain viral parents only: {parent}')
                used_routes.add(parent)
            db.execute('INSERT INTO truth(id,parent,class,length,cluster,split,route,raw) '
                       'VALUES (?,?,?,?,?,?,?,?)',
                       (key, parent, cls, length, cluster, split, route, json.dumps(row)))
            counts[cls] += 1
            bases += length
            if n % 10000 == 0:
                db.commit()
            if n % 100000 == 0:
                print(f'Audited {n:,} truth rows', file=sys.stderr, flush=True)
    db.commit()
    if not counts:
        raise ValueError('Empty truth table')
    if set(routes) - used_routes:
        raise ValueError('Routing map contains parents absent from viral truth')
    if (sum(counts.values()) != expected['retained_sequences'] or
            dict(counts) != expected['retained_by_source_class'] or
            bases != expected['total_nucleotides'] or
            len(excluded) != expected['author_exclusions_removed']):
        raise ValueError('Truth totals/exclusion counts disagree with setup summary')
    db.execute('CREATE INDEX truth_parent ON truth(parent)')
    db.commit()
    return fields, dict(counts), bases


def fasta_records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        name, seq = None, []
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                header = line[1:].strip()
                if not header:
                    raise ValueError('Empty FASTA header')
                name, seq = header.split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError('Sequence before FASTA header')
                seq.append(line.strip())
        if name is not None:
            yield name, ''.join(seq)


def build_fasta(db, fields, args, out):
    """Manifest appears only after complete one-to-one FASTA reconciliation."""
    with contextlib.ExitStack() as stack:
        handles = {}
        counts = collections.Counter()
        composition = collections.Counter()
        composition_bases = collections.Counter()
        emitted = table_writer(stack.enter_context((out / 'execution_map.tsv').open('w', newline='', encoding='utf-8')),
                               ['benchmark_batch', 'sequence_id', 'parent_sequence_id', 'input_type',
                                'evaluation_stratum', 'source_class', 'length_nt', 'sequence_sha256'])
        truths = table_writer(stack.enter_context((out / 'truth.routed.tsv').open('w', newline='', encoding='utf-8')),
                              fields + ['parent_sequence_id', 'author_cluster_split', 'input_type',
                                        'benchmark_batch', 'sequence_sha256'])
        for key, seq in fasta_records(args.fasta):
            row = db.execute('SELECT parent,length,split,route,raw,seen FROM truth WHERE id=?', (key,)).fetchone()
            if row is None or row[5]:
                raise ValueError(f'Unexpected/duplicate FASTA ID: {key}')
            if len(seq) != row[1] or not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+', seq):
                raise ValueError(f'FASTA length/alphabet mismatch: {key}')
            bucket = int(hashlib.sha256((args.seed + '\0' + row[0]).encode()).hexdigest(), 16) % args.shards
            raw = json.loads(row[4])
            sequence_hash = hashlib.sha256(seq.encode('ascii')).hexdigest()
            planned = execution_routes(row[3], getattr(args, 'unresolved_policy', 'error'))
            batches = []
            for route, stratum in planned:
                cohort = 'UNRESOLVED' if row[3] == 'unresolved' else 'KNOWN'
                batch = f'GNP_{route.upper()}_{cohort}_{bucket + 1:04d}'
                batches.append(batch)
                if batch not in handles:
                    handles[batch] = stack.enter_context((out / f'{batch}.fasta').open('w', encoding='ascii', newline='\n'))
                handles[batch].write(f'>{key}\n{seq}\n')
                emitted.writerow([batch,key,row[0],route,stratum,raw['source_class'],len(seq),sequence_hash])
                counts[batch] += 1
                composition[batch,raw['source_class']] += 1
                composition_bases[batch,raw['source_class']] += len(seq)
            # Truth stays one row per original input. Paired executions live only
            # in execution_map.tsv, never as extra biological observations.
            truths.writerow([raw[col] for col in fields] + [row[0], row[2], row[3], ';'.join(batches),sequence_hash])
            db.execute('UPDATE truth SET seen=1,batch=? WHERE id=?', (';'.join(batches), key))
            if sum(counts.values()) % 10000 < len(planned):
                db.commit()
        db.commit()
        missing = db.execute('SELECT COUNT(*) FROM truth WHERE seen=0').fetchone()[0]
        if missing:
            raise ValueError(f'{missing} truth IDs missing from FASTA')
    with (out / 'manifest.csv').open('w', encoding='utf-8', newline='') as f:
        w = table_writer(f, ['prefix', 'type', 'fasta'], ',')
        for batch in sorted(counts):
            w.writerow([batch, batch.split('_')[1].lower(), f'{batch}.fasta'])
    with (out / 'batch_composition.tsv').open('w', encoding='utf-8', newline='') as f:
        w = table_writer(f, ['batch', 'source_class', 'sequences', 'bases'])
        w.writerows((batch,cls,composition[batch,cls],composition_bases[batch,cls])
                    for batch,cls in sorted(composition))
    return dict(sorted(counts.items()))


def execution_routes(route, policy):
    if route in ('dna', 'rna'):
        return [(route, 'source_routed')]
    if route == 'unresolved' and policy == 'both':
        return [('dna', 'unresolved_dna'), ('rna', 'unresolved_rna')]
    raise ValueError('Unresolved input requires explicit --unresolved-policy both')


def prepare(args):
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    report = {'status': 'INCOMPLETE', 'pipeline_launched': False}
    try:
        routes = load_routes(args.routing)
        policy = getattr(args, 'unresolved_policy', 'error')
        if policy not in ('error', 'both'):
            raise ValueError('Invalid unresolved policy')
        with contextlib.closing(sqlite3.connect(out / 'audit.sqlite')) as db:
            with source_files(args.metadata) as opener:
                fields, counts, bases = make_index(opener, db, routes)
            report.update(counts=counts, total_nt=bases,
                          preparation_script_sha256=sha256(__file__),
                          viral_parents=db.execute("SELECT COUNT(DISTINCT parent) FROM truth WHERE class='virus'").fetchone()[0],
                          author_cluster_split_counts=dict(db.execute('SELECT split,COUNT(*) FROM truth GROUP BY split')),
                          routing_sha256=sha256(args.routing) if args.routing else None,
                          seed=args.seed, shards_per_type=args.shards,
                          metadata_source=str(Path(args.metadata).resolve()),
                          metadata_archive_sha256=sha256(args.metadata) if Path(args.metadata).is_file() else None,
                          weights='Original training weights preserved verbatim; not designated evaluation weights',
                          calibration='Enabled; benchmark uses fixed metagenome composition and default filters',
                          independence='Author training/evaluation collection; not an independent holdout')
            if Path(args.metadata).is_dir():
                report['metadata_sha256'] = {name: sha256(Path(args.metadata) / name)
                                             for name in (TRUTH, CLUSTERS, SUMMARY, EXCLUDED)}
            unresolved = db.execute("SELECT COUNT(*) FROM truth WHERE route='unresolved'").fetchone()[0]
            report['unresolved_fragments'] = unresolved
            report['unresolved_policy'] = policy
            report['unique_input_sequences'] = sum(counts.values())
            report['planned_execution_sequences'] = sum(counts.values()) + (unresolved if policy == 'both' else 0)
            report['route_fragment_counts'] = dict(db.execute('SELECT route,COUNT(*) FROM truth GROUP BY route'))
            report['evaluation_rule'] = 'One truth row per original fragment; paired unresolved DNA/RNA outcomes are separate strata, never an OR/union or best-route selection'
            with (out / 'routing_required.tsv').open('w', newline='', encoding='utf-8') as f:
                w = table_writer(f, ['parent_sequence_id', 'reference_cluster', 'fragment_count',
                                     'input_type', 'evidence_source', 'evidence_record'])
                for row in db.execute("SELECT parent,cluster,COUNT(*) FROM truth WHERE route='unresolved' GROUP BY parent,cluster ORDER BY parent"):
                    w.writerow([*row, '', '', ''])
            with (out / 'routing_used.tsv').open('w', newline='', encoding='utf-8') as f:
                w = table_writer(f, ['parent_sequence_id', 'input_type', 'evidence_source', 'evidence_record'])
                for key in sorted(routes):
                    w.writerow([key, *routes[key]])
            if args.command == 'build':
                if unresolved and policy != 'both':
                    raise ValueError(f'{unresolved} viral fragments have no source-supported DNA/RNA route; see routing_required.tsv')
                report['batches'] = build_fasta(db, fields, args, out)
                if sum(report['batches'].values()) != report['planned_execution_sequences']:
                    raise ValueError('Execution count does not reconcile with unique inputs and paired routes')
                report['input_fasta_sha256'] = sha256(args.fasta)
                report['status'] = 'INPUTS_PREPARED_NOT_LAUNCHED'
            else:
                report['status'] = ('METADATA_VALID_DUAL_ROUTING_PLANNED' if policy == 'both' else 'METADATA_VALID_ROUTING_INCOMPLETE') if unresolved else 'METADATA_AND_ROUTING_VALID'
    except Exception as exc:
        report['status'] = 'FAILED'
        report['error'] = str(exc)
        raise
    finally:
        (out / 'preparation_summary.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    with (out / 'checksums.sha256').open('w', encoding='utf-8', newline='\n') as f:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.name != 'checksums.sha256':
                f.write(f'{sha256(path)}  {path.name}\n')
    print(json.dumps(report, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['inspect', 'build'])
    parser.add_argument('--metadata', type=Path, required=True, help='Metadata tar.gz or benchmark root directory')
    parser.add_argument('--out', type=Path, required=True, help='New directory; existing directories are refused')
    parser.add_argument('--routing', type=Path, help='Source-cited viral parent-to-input-type TSV')
    parser.add_argument('--unresolved-policy', choices=['error', 'both'], default='error',
                        help='Explicitly authorize paired DNA/RNA executions of unresolved viral parents')
    parser.add_argument('--fasta', type=Path, help='Prepared FASTA, optionally gzip-compressed; required for build')
    parser.add_argument('--shards', type=int, default=32, help='Hash shards per input type, not a size limit')
    parser.add_argument('--seed', default='visum-alpha-gnp-v1')
    args = parser.parse_args()
    if args.shards < 1 or args.shards > 256:
        parser.error('--shards must be between 1 and 256')
    if args.command == 'build' and args.fasta is None:
        parser.error('build requires --fasta')
    try:
        prepare(args)
    except (ValueError, OSError, KeyError, sqlite3.IntegrityError, tarfile.TarError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
