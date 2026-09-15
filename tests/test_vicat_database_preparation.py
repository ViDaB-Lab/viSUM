import argparse
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import prepare_vicat_nonviral_build as preparation


class NonviralSourcePreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.manifest = self.root / 'manifest.tsv'
        self.manifest.write_text('assembly_accession\tcellular_group\nGCF_001.1\tbacteria\n')
        self.protein = self.root / 'package/ncbi_dataset/data/GCF_001.1/protein.faa'
        self.protein.parent.mkdir(parents=True)
        self.protein.write_text('>WP_001\nMKKLL\n')
        feature = self.root / 'metadata/bacteria/GCF_001.1/GCF_001.1.feature_table.txt.gz'
        feature.parent.mkdir(parents=True)
        with gzip.open(feature, 'wt') as handle:
            handle.write('# feature\tseq_type\tgenomic_accession\tproduct_accession\nCDS\tchromosome\tNC_001\tWP_001\n')
        self.args = argparse.Namespace(classified_dir=None, manifest=str(self.manifest),
            package_root=str(self.root / 'package'), metadata_root=str(self.root / 'metadata'),
            work_dir=str(self.root / 'work'))

    def test_real_classification_and_stable_reuse(self):
        result = preparation.classified_inputs(self.args)
        with gzip.open(result / 'cellular_chromosome.faa.gz', 'rt') as handle:
            self.assertIn('NONVIRAL|CELLULAR_CHROMOSOME|GCF_001.1|WP_001', handle.read())
        with patch.object(preparation.subprocess, 'run', side_effect=AssertionError('unexpected rebuild')):
            self.assertEqual(result, preparation.classified_inputs(self.args))

    def test_changed_inputs_rejected(self):
        preparation.classified_inputs(self.args)
        self.protein.write_text('>WP_001\nMKKLLL\n')
        with self.assertRaisesRegex(ValueError, 'inputs changed'):
            preparation.classified_inputs(self.args)

    def test_corrupted_classification_rejected(self):
        result = preparation.classified_inputs(self.args)
        (result / 'cellular_chromosome.faa.gz').write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'output changed'):
            preparation.classified_inputs(self.args)

    def test_missing_source_rejected_before_classification(self):
        self.protein.unlink()
        with self.assertRaises(FileNotFoundError):
            preparation.classified_inputs(self.args)
        self.assertFalse((self.root / 'work').exists())

    def test_manifest_path_escape_rejected(self):
        self.manifest.write_text('assembly_accession\tcellular_group\n../elsewhere\tbacteria\n')
        with self.assertRaisesRegex(ValueError, 'path components'):
            preparation.classified_inputs(self.args)

    def test_preclassified_inputs_bypass_classification(self):
        self.args.classified_dir = str(self.root / 'classified')
        with patch.object(preparation.subprocess, 'run', side_effect=AssertionError('unexpected classification')):
            self.assertEqual(preparation.classified_inputs(self.args), self.root / 'classified')

    def test_builder_invoked_with_existing_defaults(self):
        classified = self.root / 'classified'
        arguments = ['prepare', '--classified-dir', str(classified), '--destination', str(self.root / 'db'),
                     '--work-dir', str(self.root / 'work'), '--threads', '2']
        with patch.object(sys, 'argv', arguments), patch.object(preparation.subprocess, 'run') as run:
            preparation.main()
        command = run.call_args.args[0]
        self.assertIn('build_vicat_nonviral_database.sh', command[1])
        self.assertEqual(command[command.index('--classified-dir') + 1], str(classified))
        self.assertNotIn('--min-seq-id', command)


class PreparationWorkflowTests(unittest.TestCase):
    def test_both_databases_connected_and_builds_serialized(self):
        workflow = (ROOT / 'visum_nextflow.nf').read_text()
        self.assertIn('ch_vicat_nonviral_database_request = PREPARE_VICAT_DATABASE.out.database', workflow)
        self.assertIn('PREPARE_VICAT_NONVIRAL_DATABASE(ch_vicat_nonviral_database_request)', workflow)
        self.assertIn('Choose --vicat_db OR the two MetaVR source files', workflow)
        module = (ROOT / 'modules/local/vicat_nonviral_database.nf').read_text()
        for text in ('flock 9', 'prepare_vicat_nonviral_build.py', 'cache false', 'path(build_helpers)',
                     'IS DISTINCT FROM', 'reference_class IS NULL', 'sha256sum --check'):
            self.assertIn(text, module)

    def test_shell_syntax(self):
        bash = 'C:/Program Files/Git/bin/bash.exe' if os.name == 'nt' else shutil.which('bash')
        if not bash or not Path(bash).exists():
            self.skipTest('Bash unavailable')
        for name in ('vicat_database', 'vicat_nonviral_database'):
            source = (ROOT / f'modules/local/{name}.nf').read_text()
            shell = source.split('    """', 1)[1].rsplit('    """', 1)[0]
            shell = re.sub(r'(?<!\\)\$\{[^}]+\}', 'placeholder', shell).replace('\\$', '$')
            result = subprocess.run([bash, '-n'], input=shell, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        for name in ('build_vicat_database.sh', 'build_vicat_nonviral_database.sh'):
            result = subprocess.run([bash, '-n'], input=(ROOT / 'bin' / name).read_text(),
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_viral_resume_is_fingerprinted(self):
        builder = (ROOT / 'bin/build_vicat_database.sh').read_text()
        self.assertLess(builder.index('sha256sum --check -'), builder.index('checkpoint_done()'))
        self.assertIn('cmp -s "$WORK_DIR/build_config.current.tsv"', builder)
        self.assertIn('vicat_build_config.tsv', builder)


if __name__ == '__main__':
    unittest.main()
