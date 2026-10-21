#!/usr/bin/env python3
"""Evaluate frozen viSUM endpoints on prepared geNomad-derived execution strata.

No inference, relabelling, best-route selection or pooled dual-route sensitivity.
SQLite keeps the complete execution map off the Python heap; bundles are read
one batch at a time. Native-tool comparative scoring is a separate analysis.
"""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sqlite3
import tempfile

from benchmark_bundle import Bundle, ENDPOINTS, sha
from prepare_genomad_benchmark import execution_routes


def rows(path, delimiter='\t'):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def evaluate(prepared, results, output):
    prepared, results, output = map(lambda p: Path(p).resolve(), (prepared, results, output))
    if output.exists(): raise FileExistsError('Use a new evaluation directory')
    checks = {}
    for line in (prepared / 'checksums.sha256').read_text().splitlines():
        digest, name = line.split('  ', 1)
        path = (prepared / name).resolve()
        if path.parent != prepared or name in checks: raise ValueError('Unsafe/duplicate preparation checksum')
        checks[name] = digest
    for name in ('preparation_summary.json', 'truth.routed.tsv', 'execution_map.tsv', 'manifest.csv'):
        if name not in checks or sha(prepared / name) != checks[name]:
            raise ValueError('Prepared metadata checksum mismatch: ' + name)
    preparation = json.loads((prepared / 'preparation_summary.json').read_text())
    if preparation['status'] != 'INPUTS_PREPARED_NOT_LAUNCHED': raise ValueError('Incomplete preparation')
    batches = list(rows(prepared / 'manifest.csv', ','))
    if len({b['prefix'] for b in batches}) != len(batches): raise ValueError('Duplicate batch')
    batch_types = {b['prefix']: b['type'] for b in batches}
    if set(batch_types) != set(preparation['batches']): raise ValueError('Preparation batch coverage mismatch')
    if {p.name.removesuffix('_results') for p in results.glob('*_results')} != set(batch_types):
        raise ValueError('Missing or unexpected result sample directories')
    counts, identities = Counter(), {}
    with tempfile.TemporaryDirectory(prefix='visum-eval-') as temporary:
        db = sqlite3.connect(Path(temporary) / 'evaluation.sqlite')
        try:
            db.executescript('CREATE TABLE truth(id TEXT PRIMARY KEY,parent TEXT,cls TEXT,route TEXT,hash TEXT,length INTEGER);'
                'CREATE TABLE expected(id TEXT,route TEXT,stratum TEXT,PRIMARY KEY(id,route));'
                'CREATE TABLE execution(batch TEXT,id TEXT,route TEXT,stratum TEXT,cls TEXT,PRIMARY KEY(batch,id),UNIQUE(id,route));')
            n = 0
            for row in rows(prepared / 'truth.routed.tsv'):
                if row['source_class'] not in ('virus', 'chromosome', 'plasmid'):
                    raise ValueError('Invalid source class')
                if row['source_class'] != 'virus' and row['input_type'] != 'dna':
                    raise ValueError('Genomic negatives must use DNA')
                db.execute('INSERT INTO truth VALUES (?,?,?,?,?,?)', (row['sequence_id'], row['parent_sequence_id'],
                    row['source_class'], row['input_type'], row['sequence_sha256'], int(row['length_nt'])))
                for route, stratum in execution_routes(row['input_type'], preparation['unresolved_policy']):
                    db.execute('INSERT INTO expected VALUES (?,?,?)', (row['sequence_id'], route, stratum))
                n += 1
                if n % 10000 == 0: db.commit()
            emitted = 0
            for row in rows(prepared / 'execution_map.tsv'):
                truth = db.execute('SELECT parent,cls,hash,length FROM truth WHERE id=?', (row['sequence_id'],)).fetchone()
                expected = db.execute('SELECT stratum FROM expected WHERE id=? AND route=?',
                    (row['sequence_id'], row['input_type'])).fetchone()
                if (truth != (row['parent_sequence_id'], row['source_class'], row['sequence_sha256'], int(row['length_nt']))
                        or expected != (row['evaluation_stratum'],)
                        or batch_types.get(row['benchmark_batch']) != row['input_type']):
                    raise ValueError('Execution map disagrees with frozen truth/routing')
                db.execute('INSERT INTO execution VALUES (?,?,?,?,?)', (row['benchmark_batch'], row['sequence_id'],
                    row['input_type'], row['evaluation_stratum'], row['source_class']))
                emitted += 1
                if emitted % 10000 == 0: db.commit()
            if (n <= 0 or n != preparation['unique_input_sequences'] or emitted != preparation['planned_execution_sequences']
                    or db.execute('SELECT COUNT(*) FROM expected').fetchone()[0] != emitted):
                raise ValueError('Unique input / paired execution coverage mismatch')
            output.mkdir(parents=True)
            columns = ['benchmark_batch', 'sequence_id', 'input_type', 'evaluation_stratum', 'source_class', *ENDPOINTS]
            with (output / 'input_calls.tsv').open('w', encoding='utf-8', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
                writer.writeheader()
                for batch in batches:
                    sample = batch['prefix']
                    bundle = Bundle(results / (sample + '_results') / 'viharmony', sample)
                    mapping, selected = bundle.endpoints()
                    expected_rows = {r[0]: r[1:] for r in db.execute(
                        'SELECT id,route,stratum,cls FROM execution WHERE batch=?', (sample,))}
                    if set(mapping.values()) != set(expected_rows) or len(mapping) != preparation['batches'][sample]:
                        raise ValueError('Result dispositions do not cover prepared input batch exactly')
                    identities[sample] = sha(bundle.manifest_path)
                    for normalized, original in mapping.items():
                        route, stratum, cls = expected_rows[original]
                        calls = {name: int(normalized in selected[name]) for name in ENDPOINTS}
                        writer.writerow(dict(benchmark_batch=sample, sequence_id=original, input_type=route,
                            evaluation_stratum=stratum, source_class=cls, **calls))
                        for endpoint, called in calls.items():
                            counts[(stratum, route, cls, endpoint, 'denominator')] += 1
                            counts[(stratum, route, cls, endpoint, 'called')] += called
            metric_rows = []
            for stratum, route, cls, endpoint, kind in sorted(counts):
                if kind != 'denominator': continue
                denominator = counts[(stratum, route, cls, endpoint, kind)]
                called = counts[(stratum, route, cls, endpoint, 'called')]
                metric_rows.append(dict(evaluation_stratum=stratum, input_type=route, source_class=cls,
                    endpoint=endpoint, denominator=denominator, called=called, rate=called / denominator,
                    interpretation=('original_label_sensitivity' if cls == 'virus' else 'original_label_FPR')
                        if endpoint in ('primary', 'viral_candidates') else 'TE_export_rate_not_viral_detection'))
            with (output / 'metrics.tsv').open('w', encoding='utf-8', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(metric_rows[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader(); writer.writerows(metric_rows)
            report = dict(status='PASS', endpoint_version='visum-alpha-endpoints-v1', unique_inputs=n,
                execution_sequences=emitted, prepared_metadata_sha256={name: checks[name] for name in
                    ('preparation_summary.json', 'truth.routed.tsv', 'execution_map.tsv', 'manifest.csv')},
                result_manifest_sha256=identities, truth_labels_changed=False, paired_arms_pooled=False,
                input_fasta_reverification='Use verify_genomad_preparation.py separately before launching; this evaluator verifies metadata/output checksums.',
                comparative_native_tool_metrics='Not generated here; retain raw native outputs for separate comparative scoring.')
            (output / 'evaluation_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
            return report
        finally:
            db.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    evaluate(args.prepared, args.results, args.outdir)
