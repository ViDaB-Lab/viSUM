#!/usr/bin/env python3
"""Export/checksum a frozen-vOTU compatibility patch; apply only to a new copy.

This is an opt-in standalone prototype, NOT a portal data correction or a new
taxonomy rule. It restores the frozen reference's normalized vOTU assignments.
No network, protein, database-building, or deployment support is present.
"""
import csv
import hashlib
import json
from pathlib import Path
import tempfile

from build_vicat_taxonomy_lookup import RANKS, build_metadata_sql, clean_text

PATCH_COLUMNS = ('uvig', 'portal_votu_key', 'frozen_votu_key')


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def sql(value):
    return "'" + str(value).replace("'", "''") + "'"


def raw_reader(path):
    return f"read_csv({sql(path)}, delim='\t', header=true, all_varchar=true, sample_size=10000, nullstr='\\N')"


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n', encoding='utf-8')


def write_tsv(path, columns, rows):
    with Path(path).open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(columns)
        writer.writerows(rows)


def export_patch(original, converted, output, work, expected_rows, expected_changes,
                 source_identity, connect):
    original, converted, output, work = map(Path, (original, converted, output, work))
    output.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)
    patch_file = output / 'votu_compatibility_patch.tsv'
    differences = output / 'votu_differences.tsv'
    manifest_file = output / 'patch_manifest.json'
    if any(path.exists() for path in (patch_file, differences, manifest_file)):
        raise FileExistsError('Refusing to overwrite an existing patch/export')
    ranks = [rank for rank, _, _ in RANKS]
    with tempfile.TemporaryDirectory(prefix='patch-export-', dir=work.resolve()) as temporary:
        with connect(work) as connection:
            connection.execute('ATTACH ' + sql(Path(temporary) / 'audit.duckdb') + ' AS audit')
            metrics = {}
            for label, path in (('original', original), ('portal', converted)):
                connection.execute(f'CREATE TABLE audit.{label} AS ' + build_metadata_sql(path))
                count, unique = connection.execute(f'SELECT count(*), count(DISTINCT uvig) FROM audit.{label}').fetchone()
                if count != expected_rows or count != unique:
                    raise ValueError('Unexpected row count or duplicate UViGs in ' + label)
                metrics[label + '_rows'] = count
                metrics[label + '_votu_groups'] = connection.execute(f'SELECT count(DISTINCT votu_key) FROM audit.{label}').fetchone()[0]
            changed_ranks = ' OR '.join(f'a.{rank} IS DISTINCT FROM b.{rank}' for rank in ranks)
            connection.execute(f'''CREATE TABLE audit.joined AS SELECT coalesce(a.uvig,b.uvig) AS uvig,
                a.uvig IS NULL AS missing_original, b.uvig IS NULL AS missing_portal,
                a.votu_key AS frozen_votu_key, b.votu_key AS portal_votu_key,
                ({changed_ranks}) AS rank_changed,
                a.taxonomy_method IS DISTINCT FROM b.taxonomy_method AS method_changed
                FROM audit.original a FULL OUTER JOIN audit.portal b USING (uvig)''')
            counts = connection.execute('''SELECT
                count(*) FILTER (WHERE missing_original), count(*) FILTER (WHERE missing_portal),
                count(*) FILTER (WHERE rank_changed), count(*) FILTER (WHERE method_changed),
                count(*) FILTER (WHERE frozen_votu_key IS DISTINCT FROM portal_votu_key)
                FROM audit.joined''').fetchone()
            metrics.update(zip(('missing_original', 'missing_portal', 'rank_changed_rows',
                                'method_changed_rows', 'votu_changed_rows'), counts))
            if any(counts[:4]):
                raise ValueError('Patch refused: missing UViGs or taxonomy/method changes: ' + json.dumps(metrics))
            if counts[4] != expected_changes:
                raise ValueError(f'Expected exactly {expected_changes} vOTU differences, observed {counts[4]}')
            # Bound the in-memory export: production expectation is 1,409 records.
            if not 0 <= expected_changes <= 10000:
                raise ValueError('Patch export is restricted to at most 10,000 reviewed differences')
            connection.execute('''CREATE TABLE audit.pairs AS SELECT frozen_votu_key,portal_votu_key,count(*) AS uvigs
                FROM audit.joined GROUP BY frozen_votu_key,portal_votu_key''')
            metrics['split_original_groups'] = connection.execute('''SELECT count(*) FROM (
                SELECT frozen_votu_key FROM audit.pairs GROUP BY frozen_votu_key HAVING count(*)>1)''').fetchone()[0]
            metrics['merged_portal_groups'] = connection.execute('''SELECT count(*) FROM (
                SELECT portal_votu_key FROM audit.pairs GROUP BY portal_votu_key HAVING count(*)>1)''').fetchone()[0]
            rows = connection.execute('''SELECT uvig,portal_votu_key,frozen_votu_key FROM audit.joined
                WHERE frozen_votu_key IS DISTINCT FROM portal_votu_key ORDER BY uvig''').fetchall()
            review_rows = connection.execute(f'''SELECT j.uvig,j.portal_votu_key,j.frozen_votu_key,
                a.taxonomy_method, {', '.join('a.' + rank for rank in ranks)}
                FROM audit.joined j JOIN audit.original a USING (uvig)
                WHERE j.frozen_votu_key IS DISTINCT FROM j.portal_votu_key ORDER BY j.uvig''').fetchall()
    write_tsv(patch_file, PATCH_COLUMNS, rows)
    write_tsv(differences, PATCH_COLUMNS + ('taxonomy_method',) + tuple(ranks),
              (tuple('NA' if value is None else value for value in row) for row in review_rows))
    manifest = dict(schema='visum-metavr-votu-compatibility-v1',
        purpose='Restore frozen reference grouping, NOT declare portal assignments incorrect',
        original_metadata_sha256=digest(original), converted_portal_sha256=digest(converted),
        source_identity=source_identity, expected_rows=expected_rows, expected_changes=expected_changes,
        patch_sha256=digest(patch_file), differences_sha256=digest(differences),
        columns=list(PATCH_COLUMNS), before=metrics, full_source_equivalence=False,
        protein_equivalence='NOT TESTED', deployment_performed=False)
    json_write(manifest_file, manifest)
    return manifest


def load_patch(patch_file, manifest_file, converted):
    manifest = json.loads(Path(manifest_file).read_text(encoding='utf-8'))
    if manifest.get('schema') != 'visum-metavr-votu-compatibility-v1':
        raise ValueError('Unsupported patch manifest schema')
    if digest(patch_file) != manifest['patch_sha256']:
        raise ValueError('Compatibility patch checksum mismatch')
    if digest(converted) != manifest['converted_portal_sha256']:
        raise ValueError('Converted metadata checksum differs from the patch source')
    with Path(patch_file).open(encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if reader.fieldnames != list(PATCH_COLUMNS):
            raise ValueError('Unexpected patch columns')
        rows = list(reader)
    if len(rows) != manifest['expected_changes'] or len(rows) > 10000:
        raise ValueError('Unexpected number of compatibility patch records')
    if len({row['uvig'] for row in rows}) != len(rows):
        raise ValueError('Duplicate UViGs in compatibility patch')
    for row in rows:
        if set(row) != set(PATCH_COLUMNS) or any(not value or not value.strip() for value in row.values()):
            raise ValueError('Missing or malformed compatibility patch value')
        if row['portal_votu_key'] == row['frozen_votu_key']:
            raise ValueError('Compatibility patch contains an unchanged assignment')
    return manifest, rows


def apply_patch_copy(converted, patch_file, manifest_file, destination, work, connect):
    converted, patch_file, manifest_file, destination, work = map(Path,
        (converted, patch_file, manifest_file, destination, work))
    if destination.exists() or destination.with_name(destination.name + '.part').exists():
        raise FileExistsError('Refusing to overwrite patched metadata or a partial output')
    if destination.resolve() in (converted.resolve(), patch_file.resolve(), manifest_file.resolve()):
        raise ValueError('Patched metadata must be a new private copy')
    manifest, rows = load_patch(patch_file, manifest_file, converted)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.part')
    with connect(work) as connection:
        connection.execute('CREATE TEMP TABLE patch (uvig VARCHAR,portal_votu_key VARCHAR,frozen_votu_key VARCHAR)')
        if rows:
            connection.executemany('INSERT INTO patch VALUES (?,?,?)',
                [tuple(row[name] for name in PATCH_COLUMNS) for row in rows])
        matches, wrong = connection.execute('''SELECT count(m.uvig),
            count(*) FILTER (WHERE m.votu_key IS DISTINCT FROM p.portal_votu_key)
            FROM patch p LEFT JOIN (''' + build_metadata_sql(converted) + ''') m USING (uvig)''').fetchone()
        if matches != len(rows) or wrong:
            raise ValueError('Patch precondition failed: UViG absent/duplicated or portal assignment changed')
        query = f'''SELECT r.* REPLACE (
            CASE WHEN p.uvig IS NOT NULL THEN p.frozen_votu_key ELSE r.votu END AS votu)
            FROM {raw_reader(converted)} r LEFT JOIN patch p ON {clean_text('r.uvig')}=p.uvig'''
        connection.execute(f'COPY ({query}) TO {sql(partial)} '
                           "(FORMAT CSV, DELIMITER '\t', HEADER true, COMPRESSION gzip, NULL '\\N')")
    partial.replace(destination)
    return dict(applied_rows=len(rows), patched_sha256=digest(destination),
                patch_sha256=manifest['patch_sha256'], output=str(destination))


def verify_non_votu_columns(converted, patched, work, expected_rows, connect):
    """Full raw-field set comparison: not merely normalized taxonomy equality."""
    with connect(Path(work)) as connection:
        before = [row[0] for row in connection.execute(f'DESCRIBE SELECT * FROM {raw_reader(converted)}').fetchall()]
        after = [row[0] for row in connection.execute(f'DESCRIBE SELECT * FROM {raw_reader(patched)}').fetchall()]
        if before != after or 'uvig' not in before or 'votu' not in before:
            raise ValueError('Patched metadata schema differs')
        # A bounded-memory hash join compares every retained raw column by UViG.
        compared = [name for name in before if name != 'votu']
        def quote(name):
            return '"' + name.replace('"', '""') + '"'
        changed = ' OR '.join('a.' + quote(name) + ' IS DISTINCT FROM b.' + quote(name) for name in compared)
        count, different = connection.execute(f'''SELECT count(*),count(*) FILTER (WHERE {changed})
            FROM {raw_reader(converted)} a FULL OUTER JOIN {raw_reader(patched)} b USING (uvig)''').fetchone()
        if count != expected_rows or different:
            raise ValueError('Non-vOTU fields changed or rows missing/duplicated')
        return dict(compared_columns=compared, joined_rows=count, changed_rows=different)
