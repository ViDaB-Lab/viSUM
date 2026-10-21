"""Integrity checks for documentation assets, not classifier changes."""
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from recover_vicat_nonviral_snapshot import GROUPS, recover, recover_catalogue, summary_metrics

SNAPSHOT = Path(__file__).resolve().parents[1] / 'assets/vicat/nonviral_classaware_v1_20260904'


def test_recovered_manifest_is_complete_versioned_and_build_compatible():
    import re
    sys.path.insert(0, str(SNAPSHOT.parents[2] / 'bin'))
    from classify_vicat_nonviral_references import read_manifest
    path = SNAPSHOT / 'assembly_manifest.tsv'
    rows = read_manifest(path)
    assert len(rows) == len({r['assembly_accession'] for r in rows}) == 480
    assert all(re.fullmatch(r'GCF_\d+\.\d+', r['assembly_accession']) for r in rows)
    assert all(r['organism_name'] and r['taxid'].isdigit() for r in rows)
    assert Counter(r['cellular_group'] for r in rows) == GROUPS
    provenance = json.loads((SNAPSHOT / 'provenance.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == provenance['manifest']['sha256']


def test_build_evidence_matches_recorded_hashes_and_totals():
    provenance = json.loads((SNAPSHOT / 'provenance.json').read_text())
    for name in ('vicat_nonviral_build_config.tsv', 'vicat_nonviral_database_metadata.tsv',
                 'vicat_nonviral_cluster_summary.tsv', 'vicat_nonviral_classification_summary.tsv'):
        data = (SNAPSHOT / 'build_evidence' / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == provenance['database_files_sha256'][name]
    stats = provenance['inventory']
    summary = summary_metrics((SNAPSHOT / 'build_evidence/vicat_nonviral_classification_summary.tsv').read_text())
    assert stats['included_protein_count'] == summary['included_protein_count', 'all'] == 3731712
    assert summary['excluded_protein_count', 'all'] == 0
    assert sum(stats['proteins_by_group'].values()) == sum(stats['proteins_by_class'].values()) == 3731712
    with (SNAPSHOT / 'build_evidence/vicat_nonviral_cluster_summary.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert sum(int(r['input_proteins']) for r in rows) == 3731712
    assert sum(int(r['representatives']) for r in rows) == provenance['build']['representatives'] == 2708478


def test_unmatched_summary_cannot_be_confused_with_build_summary():
    provenance = json.loads((SNAPSHOT / 'provenance.json').read_text())
    unmatched = provenance['unmatched_summary']
    assert unmatched['status'] == 'unmatched_not_a_build_input_or_authoritative_summary'
    data = (SNAPSHOT / unmatched['filename']).read_bytes()
    assert hashlib.sha256(data).hexdigest() == unmatched['sha256']
    other = summary_metrics(data.decode())
    matching = summary_metrics((SNAPSHOT / 'build_evidence/vicat_nonviral_classification_summary.tsv').read_text())
    differences = {category: matching[key, category] - value
                   for (key, category), value in other.items() if key == 'reference_class_count'}
    assert differences == {'CELLULAR_CHROMOSOME': 0, 'CELLULAR_UNPLACED': 64,
                           'PLASMID': 150, 'PLASTID': 186, 'MITOCHONDRIAL': 3, 'SHARED_NONVIRAL': 1}
    assert sum(differences.values()) == other['excluded_protein_count', 'all'] == 404


def test_recovery_rejects_wrong_source_before_reading_other_files(tmp_path):
    source = tmp_path / 'wrong.gz'
    source.write_bytes(b'wrong snapshot')
    before = source.read_bytes()
    with pytest.raises(ValueError, match='checksum'):
        recover(source, tmp_path / 'missing.tar.gz', tmp_path / 'missing.tsv')
    assert source.read_bytes() == before


@pytest.mark.parametrize('accession,second_organism,message', [
    ('GCF_000000001', 'First organism', 'versioned'),
    ('GCF_000000001.1', 'Different organism', 'Conflicting'),
])
def test_recovery_rejects_invalid_or_conflicting_annotations(tmp_path, accession, second_organism, message):
    source = tmp_path / 'source.tsv.gz'
    header = 'source_accession\tcellular_group\torganism_name\ttaxid\treference_class\tclassification_note\n'
    rows = [f'{accession}\tbacteria\t{organism}\t123\tCELLULAR_CHROMOSOME\t\n'
            for organism in ('First organism', second_organism)]
    source.write_bytes(gzip.compress((header + ''.join(rows)).encode()))
    with pytest.raises(ValueError, match=message):
        recover_catalogue(source)


def test_summary_rejects_duplicate_keys():
    with pytest.raises(ValueError, match='Duplicate'):
        summary_metrics('metric\tcategory\tvalue\nx\tall\t1\nx\tall\t2\n')
