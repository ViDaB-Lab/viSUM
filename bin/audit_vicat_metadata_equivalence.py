#!/usr/bin/env python3
"""Prove metadata equivalence modulo a BIJECTIVE renaming of vOTU labels.

Labels are not ignored: every UViG and consumed taxonomy/method value must
match, and the complete old/new vOTU partitions must map one-to-one. Splits,
merges, missing/duplicate IDs, and changed taxonomy/methods always reject.
This audits metadata only; it never updates references or biological policy.
"""
from pathlib import Path
import tempfile

from build_vicat_taxonomy_lookup import RANKS, build_metadata_sql


def sql(value):
    return "'" + str(value).replace("'", "''") + "'"


def examples(connection, query):
    result = connection.execute(query)
    names = [column[0] for column in result.description]
    return [dict(zip(names, row)) for row in result.fetchall()]


def audit_metadata_equivalence(original, converted, work, expected_rows, connect):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    ranks = [rank for rank, _, _ in RANKS]
    result = dict(status='INCOMPATIBLE', builder_semantics_equivalent=False,
        exact_builder_metadata_match=False, protein_equivalence='NOT TESTED',
        full_source_equivalence=False, database_build_performed=False,
        proof_scope='all normalized UViGs/ranks/methods and full bidirectional vOTU partition mapping')
    # Disk-backed audit tables permit bounded memory. These are private temporary
    # metadata tables, NOT a viCAT database or representative-protein build.
    with tempfile.TemporaryDirectory(prefix='metavr-partition-audit-', dir=work.resolve()) as temporary:
        with connect(work) as connection:
            connection.execute('ATTACH ' + sql(Path(temporary) / 'audit.duckdb') + ' AS audit')
            for label, path in (('original', original), ('portal', converted)):
                connection.execute(f'CREATE TABLE audit.{label}_metadata AS ' + build_metadata_sql(Path(path)))
                result[label + '_rows'] = connection.execute(f'SELECT count(*) FROM audit.{label}_metadata').fetchone()[0]
                result[label + '_duplicate_uvig_groups'] = connection.execute(f'''SELECT count(*) FROM (
                    SELECT uvig FROM audit.{label}_metadata GROUP BY uvig HAVING count(*) > 1)''').fetchone()[0]
                result[label + '_votu_groups'] = connection.execute(f'''SELECT count(*) FROM (
                    SELECT votu_key FROM audit.{label}_metadata GROUP BY votu_key)''').fetchone()[0]
            if (result['original_rows'] != expected_rows or result['portal_rows'] != expected_rows or
                    result['original_duplicate_uvig_groups'] or result['portal_duplicate_uvig_groups']):
                result['reason'] = 'Unexpected normalized row count or duplicate UViG identifiers; no renaming can be approved'
                return result
            changed_ranks = ' OR '.join(f'a.{rank} IS DISTINCT FROM b.{rank}' for rank in ranks)
            connection.execute(f'''CREATE TABLE audit.joined AS
                SELECT coalesce(a.uvig, b.uvig) AS uvig,
                    a.uvig IS NULL AS missing_from_original,
                    b.uvig IS NULL AS missing_from_portal,
                    a.votu_key AS original_votu, b.votu_key AS portal_votu,
                    ({changed_ranks}) AS rank_changed,
                    a.taxonomy_method IS DISTINCT FROM b.taxonomy_method AS method_changed
                FROM audit.original_metadata a FULL OUTER JOIN audit.portal_metadata b USING (uvig)''')
            counts = connection.execute('''SELECT
                count(*) FILTER (WHERE missing_from_original),
                count(*) FILTER (WHERE missing_from_portal),
                count(*) FILTER (WHERE NOT missing_from_original AND NOT missing_from_portal AND rank_changed),
                count(*) FILTER (WHERE NOT missing_from_original AND NOT missing_from_portal AND method_changed),
                count(*) FILTER (WHERE NOT missing_from_original AND NOT missing_from_portal AND original_votu IS DISTINCT FROM portal_votu)
                FROM audit.joined''').fetchone()
            names = ('missing_from_original', 'missing_from_portal', 'rank_changed_rows',
                     'method_changed_rows', 'votu_label_changed_rows')
            result.update(zip(names, counts))
            connection.execute('''CREATE TABLE audit.votu_pairs AS
                SELECT original_votu, portal_votu, count(*) AS uvig_count
                FROM audit.joined WHERE NOT missing_from_original AND NOT missing_from_portal
                GROUP BY original_votu, portal_votu''')
            result['mapping_pairs'] = connection.execute('SELECT count(*) FROM audit.votu_pairs').fetchone()[0]
            result['split_original_groups'] = connection.execute('''SELECT count(*) FROM (
                SELECT original_votu FROM audit.votu_pairs GROUP BY original_votu HAVING count(*) > 1)''').fetchone()[0]
            result['merged_portal_groups'] = connection.execute('''SELECT count(*) FROM (
                SELECT portal_votu FROM audit.votu_pairs GROUP BY portal_votu HAVING count(*) > 1)''').fetchone()[0]
            result['renaming_examples'] = examples(connection, '''SELECT original_votu, portal_votu, uvig_count
                FROM audit.votu_pairs WHERE original_votu IS DISTINCT FROM portal_votu LIMIT 20''')
            result['split_examples'] = examples(connection, '''SELECT original_votu, count(*) AS portal_groups
                FROM audit.votu_pairs GROUP BY original_votu HAVING count(*) > 1 LIMIT 20''')
            result['merge_examples'] = examples(connection, '''SELECT portal_votu, count(*) AS original_groups
                FROM audit.votu_pairs GROUP BY portal_votu HAVING count(*) > 1 LIMIT 20''')
            if any(result[name] for name in ('missing_from_original', 'missing_from_portal', 'rank_changed_rows', 'method_changed_rows')):
                fields = ['votu_key', 'taxonomy_method'] + ranks
                result['taxonomic_difference_examples'] = examples(connection, f'''SELECT
                    j.uvig, j.missing_from_original, j.missing_from_portal,
                    {', '.join('a.' + field + ' AS original_' + field + ', b.' + field + ' AS portal_' + field for field in fields)}
                    FROM audit.joined j LEFT JOIN audit.original_metadata a USING (uvig)
                    LEFT JOIN audit.portal_metadata b USING (uvig)
                    WHERE j.missing_from_original OR j.missing_from_portal OR j.rank_changed OR j.method_changed LIMIT 20''')
            invalid = any(result[name] for name in ('missing_from_original', 'missing_from_portal',
                'rank_changed_rows', 'method_changed_rows', 'split_original_groups', 'merged_portal_groups'))
            bijection = result['mapping_pairs'] == result['original_votu_groups'] == result['portal_votu_groups']
            if not invalid and bijection:
                result['builder_semantics_equivalent'] = True
                result['exact_builder_metadata_match'] = result['votu_label_changed_rows'] == 0
                result['status'] = 'EXACT_MATCH' if result['exact_builder_metadata_match'] else 'EQUIVALENT_VOTU_RELABELING'
                result['reason'] = ('All consumed values and group labels match' if result['exact_builder_metadata_match'] else
                    'Labels differ, but every consumed taxonomy/method value matches and every complete vOTU group maps bijectively')
            else:
                result['reason'] = 'UViGs, consumed taxonomy/method values, or vOTU group membership differ; reference substitution remains blocked'
    return result
