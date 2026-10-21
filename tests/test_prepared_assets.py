"""Readiness checks never install, modify assets, or classify sequences."""

import base64
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
sys.path.insert(0, str(ROOT / 'tests'))
import validate_prepared_assets as validation
import test_reference_database_validation as reference_fixtures


class PreparedAssetsTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.fasta = self.asset('input.fasta', b'>c1\nACGT\n')

    def asset(self, name, contents=b'fixture\n'):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(contents)
        return path

    def spec(self, assets=None, kind='dna'):
        return {'samples': [{'prefix': 'sample', 'type': kind, 'fasta': str(self.fasta)}], 'assets': assets or {}}

    def fixture(self, tool):
        root = self.root / tool
        root.mkdir(exist_ok=True)
        spec = {'path': str(root)}
        if tool in ('genomad', 'vcontact3'):
            builder = reference_fixtures.ReferenceValidationTests()
            builder.root = root
            if tool == 'genomad':
                builder.genomad_fixture()
            else:
                builder.vc_fixture()
                spec.update(release='232', domain='both')
            return spec
        names = {
            'virsorter2': ['hmm/viral/combined.hmm', 'conda_envs/env/conda-meta/history'] +
                [f'group/{g}/model' for g in ('dsDNAphage', 'ssDNA', 'NCLDV', 'RNA', 'lavidaviridae')],
            'deepmicroclass2': ['predict.py', 'models.py', 'utils.py', 'selftest.py', 'model/300/DeepMicroClass-best.ckpt'] +
                [f'model/8-class/DeepMicroClass-best-{n}-8class.ckpt' for n in (500, 1000, 3000)],
            'virbot': ['virbot/VirBot.py', 'virbot/__init__.py'] + ['virbot/data/ref/' + n for n in
                ('VirBot.hmm', 'VirBot.dmnd', 'VirBot_hmm_threshold.txt', 'VirBot_hmm_taxa_full.txt', 'VirBot_RNAvirus_acc.txt')],
            'gianthunter': ['database.dmnd', 'names.csv', 'nodes.csv', 'ProkaryoticGroup.pkl', 'proteins.csv',
                'RefVirus.dmnd', 'RefVirus.faa', 'RefVirus_anno.pkl', 'taxid.csv', 'transformer.pth'],
            'checkv': ['genome_db/checkv_reps.' + ext for ext in ('faa', 'fna', 'tsv', 'dmnd')] +
                ['hmm_db/checkv_hmms.tsv', 'hmm_db/genome_lengths.tsv', 'hmm_db/checkv_hmms/test.hmm'],
            'vicat': ['diamond/IMGVR5_UViG_representatives.dmnd', 'taxonomy/IMGVR5_UViG.vicat_taxonomy_lookup.parquet'],
            'vicat_nonviral': ['vicat_nonviral.dmnd', 'vicat_nonviral_representatives.faa.gz',
                'vicat_nonviral_representative_metadata.parquet', 'vicat_nonviral_database_metadata.tsv',
                'vicat_nonviral_cluster_summary.tsv', 'SHA256SUMS', '.visum_db_complete'],
            'vitap': ['VMR_genome_fixture.' + ext for ext in ('gff', 'fasta', 'faa', 'dmnd')] +
                ['uniref90.dmnd', 'uniref90.accession2taxid', 'taxdmp/names.dmp', 'taxdmp/nodes.dmp', 'taxonomy.csv'] +
                [rank + '_genome.threshold' for rank in ('Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')],
            'deep6': ['Master/deep6.py', 'Master/deep6_functions.py'],
            'cenotetaker3': ['viral_cdds_and_pfams_191028.txt'] +
                [f'hmmscan_DBs/v3.1.1/{h}.h3m' for h in ('Virion_HMMs', 'DNA_rep_HMMs', 'RDRP_HMMs', 'Useful_Annotation_HMMs', 'phrogs_for_ct')] +
                [f'mmseqs_DBs/{db}{suffix}' for db in ('ct3_hallmark.taxDB', 'refseq_virus_prot_taxDB', 'CDD') for suffix in ('', '.dbtype')] +
                [f'mmseqs_DBs/{db}{suffix}' for db in ('ct3_hallmark.taxDB', 'refseq_virus_prot_taxDB') for suffix in ('_mapping', '_taxonomy')],
        }[tool]
        for name in names:
            self.asset(f'{tool}/{name}')
        if tool == 'deep6':
            spec['models'] = str(root / 'Models')
            for length in (250, 500, 1000, 1500):
                self.asset(f'deep6/Models/model_{length}_fixture.h5', b'\x89HDF\r\n\x1a\nfixture')
        if tool == 'cenotetaker3':
            spec['hmm_version'] = 'v3.1.1'
        return spec

    def test_complete_assets_for_all_twelve_tools(self):
        for tool in ('genomad', 'virsorter2', 'cenotetaker3', 'deep6', 'deepmicroclass2', 'virbot',
                     'gianthunter', 'checkv', 'vicat', 'vicat_nonviral', 'vitap', 'vcontact3'):
            with self.subTest(tool=tool):
                asset = self.fixture(tool)
                self.assertEqual(validation.validate_asset(tool, asset), str(Path(asset['path']).resolve()))

    def test_required_assets_cannot_be_missing_or_empty(self):
        spec = self.fixture('genomad')
        (Path(spec['path']) / 'genomad_db_h.index').write_bytes(b'')
        with self.assertRaisesRegex(ValueError, 'Run Nextflow with --setup'):
            validation.validate(self.spec({'genomad': spec}))

    def test_all_missing_assets_are_reported_without_creating_directories(self):
        assets = {name: {'path': str(self.root / name)} for name in ('genomad', 'checkv')}
        with self.assertRaises(ValueError) as error:
            validation.validate(self.spec(assets))
        self.assertIn('genomad:', str(error.exception))
        self.assertIn('checkv:', str(error.exception))
        self.assertFalse((self.root / 'genomad').exists())
        self.assertFalse((self.root / 'checkv').exists())

    def test_read_only_validation_keeps_contents_and_mtimes(self):
        asset = self.fixture('genomad')
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()}
        validation.validate(self.spec({'genomad': asset}))
        after = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(before, after)

    def test_molecule_specific_assets_only_required_for_present_types(self):
        missing = {'path': str(self.root / 'absent')}
        self.assertEqual(validation.validate(self.spec({'deep6': missing, 'virbot': missing}))['assets'], {})
        self.assertEqual(validation.validate(self.spec({'deepmicroclass2': missing, 'gianthunter': missing}, 'rna'))['assets'], {})
        with self.assertRaises(ValueError):
            validation.validate(self.spec({'deep6': missing}, 'rna'))

    def test_empty_batch_duplicate_prefix_bad_type_and_missing_input_fail(self):
        cases = [dict(samples=[], assets={}), self.spec(), self.spec(kind='protein'), self.spec()]
        cases[1]['samples'] *= 2
        cases[3]['samples'][0]['fasta'] = str(self.root / 'absent.fasta')
        for case in cases:
            with self.assertRaises((OSError, ValueError)):
                validation.validate(case)

    def test_hdf5_duplicate_models_and_lfs_pointers_are_rejected(self):
        deep6 = self.fixture('deep6')
        self.asset('deep6/Models/model_250_other.h5')
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            validation.validate_asset('deep6', deep6)
        virbot = self.fixture('virbot')
        self.asset('virbot/virbot/data/ref/VirBot.dmnd', b'version https://git-lfs.github.com/spec/v1\n')
        with self.assertRaisesRegex(ValueError, 'Git-LFS'):
            validation.validate_asset('virbot', virbot)

    def test_vicat_accepts_runtime_and_development_layouts_not_combined_only(self):
        asset = self.fixture('vicat')
        root = Path(asset['path'])
        (root / 'diamond/IMGVR5_UViG_representatives.dmnd').unlink()
        self.asset('vicat/vicat_viral_cellular.dmnd')
        with self.assertRaises(ValueError):
            validation.validate_asset('vicat', asset)
        self.asset('vicat/IMGVR5_UViG_representatives.dmnd')
        with self.assertRaisesRegex(ValueError, 'complete viCAT'):
            validation.validate_asset('vicat', asset)
        self.asset('vicat/IMGVR5_UViG.vicat_taxonomy_lookup.parquet')
        validation.validate_asset('vicat', asset)

    def test_cli_reports_success_and_returns_failure_without_partial_report(self):
        report = self.root / 'report.json'
        spec = self.spec({'genomad': self.fixture('genomad')})
        encoded = base64.b64encode(json.dumps(spec).encode()).decode()
        command = [sys.executable, '-B', str(ROOT / 'bin/validate_prepared_assets.py'), '--spec-base64', encoded, '--report', str(report)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(report.read_text())['validation'], 'passed')
        report.unlink()
        (Path(spec['assets']['genomad']['path']) / 'version.txt').unlink()
        result = subprocess.run(command, capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(report.exists())
        self.assertIn('--setup', result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
