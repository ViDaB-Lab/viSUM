"""Portable checks for the Linux harness; do not pretend to execute Nextflow."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import run_vicat_nonviral_linux_validation as probe


class NativeProbeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_fixture_has_checksums_for_two_sources_and_all_six_classes(self):
        manifest = probe.make_remote_fixture(self.root / 'remote')
        self.assertIn('GCF_000000001.1', manifest)
        checks = next((self.root / 'remote').rglob('md5checksums.txt'))
        for line in checks.read_text().splitlines():
            expected, name = line.split('  ./')
            self.assertEqual(probe.acquisition.digest(checks.parent / name, 'md5'), expected)
        import gzip
        feature = gzip.decompress(next((self.root / 'remote').rglob('*_feature_table.txt.gz')).read_bytes()).decode()
        for kind in ('chromosome', 'unplaced scaffold', 'plasmid', 'plastid', 'mitochondrion'):
            self.assertIn(kind, feature)
        self.assertEqual(feature.count('\tP5\n'), 2)

    def test_real_nonviral_module_biology_and_environment_are_not_substituted(self):
        manifest = probe.make_remote_fixture(self.root / 'remote')
        project = self.root / 'project'
        project.mkdir()
        changed = probe.prepare_fixture_project(project, 'http://127.0.0.1:12345/', manifest)
        self.assertEqual(len(changed), 3)
        self.assertNotIn('modules/local/vicat_nonviral_database.nf', changed)
        for name in ('modules/local/vicat_nonviral_database.nf', 'workflows/database_setup.nf',
                     'bin/classify_vicat_nonviral_references.py', 'bin/build_vicat_nonviral_database.sh',
                     'bin/prepare_vicat_nonviral_cluster_metadata.py', 'envs/vicat_db_build.yml'):
            self.assertEqual((project / name).read_bytes(), (probe.ROOT / name).read_bytes())

    def test_config_and_command_are_private_setup_only_with_explicit_resume(self):
        manifest = probe.make_remote_fixture(self.root / 'remote')
        project = self.root / 'project'
        project.mkdir()
        probe.prepare_fixture_project(project, 'http://127.0.0.1:12345/', manifest)
        config = probe.fixture_config(project, self.root)
        for value in ('params.run_genomad = false', 'params.run_vitap = false', 'conda.useMamba = true',
                      'params.max_cpus = 2', 'params.max_memory = \'8 GB\'', 'RUN_GENOMAD',
                      'beforeScript', 'conda.cacheDir'):
            self.assertIn(value, config)
        command = probe.native_command(Path('/nextflow'), project, self.root, self.root / 'c.config', 'reuse',
                                       '12345678-1234-1234-1234-123456789abc')
        self.assertIn('--setup', command)
        self.assertEqual(command[-2:], ['-resume', '12345678-1234-1234-1234-123456789abc'])
        self.assertNotIn('--input', command)
        self.assertNotIn('-profile', command)
        with self.assertRaisesRegex(ValueError, 'resume session'):
            probe.native_command(Path('/nextflow'), project, self.root, self.root / 'c.config', 'reuse', 'latest')

    def test_trace_rejects_cached_setup_and_unrelated_tasks(self):
        work = self.root / 'work'
        row = dict(task_id='1', hash='aa/bb', name='DATABASE_SETUP:PREPARE_VICAT_DATABASE (fixture)',
                   status='COMPLETED', exit='0', workdir=str(work / 'aa/bb'))
        second = dict(row, name='DATABASE_SETUP:PREPARE_VICAT_NONVIRAL_DATABASE (fixture)')
        path = self.root / 'trace.tsv'

        def write(rows):
            with path.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)

        write([row, second])
        self.assertEqual(len(probe.check_trace(path, work)), 2)
        for invalid in ([row, dict(second, status='CACHED')], [row, dict(second, name='RUN_GENOMAD')],
                        [row, second, row], [row, dict(second, workdir=str(self.root / 'outside'))]):
            write(invalid)
            with self.assertRaises(AssertionError):
                probe.check_trace(path, work)

    def test_database_snapshot_detects_byte_or_mtime_changes(self):
        file = self.root / 'reference'
        file.write_text('first')
        before = probe.snapshot_database(self.root)
        file.write_text('different')
        self.assertNotEqual(before, probe.snapshot_database(self.root))

    def test_nonviral_setup_metadata_publication_is_explicitly_refreshed(self):
        module = (probe.ROOT / 'modules/local/vicat_nonviral_database.nf').read_text()
        publication = module.split('    publishDir', 1)[1].split('    input:', 1)[0]
        self.assertIn('overwrite: true', publication)
        self.assertIn('failOnError: true', publication)
        self.assertIn("pattern: 'vicat_nonviral_database_setup_metadata.tsv'", publication)

    def test_private_environment_does_not_point_to_installed_caches(self):
        old = self.root / 'installed'
        framework = old / 'framework' / probe.FRAMEWORK_VERSION
        framework.mkdir(parents=True)
        (framework / f'nextflow-{probe.FRAMEWORK_VERSION}-one.jar').write_text('fixture')
        output = self.root / 'new'
        environment = probe.private_environment(output, old, Path('/usr/bin/nextflow'))
        for key in ('NXF_HOME', 'NXF_PLUGINS_DIR', 'CONDA_PKGS_DIRS', 'MAMBA_ROOT_PREFIX', 'XDG_CACHE_HOME'):
            self.assertTrue(str(environment[key]).startswith(str(output)), key)
        self.assertNotIn('NXF_CACHE_DIR', environment)
        self.assertNotIn('NXF_PARAMS_FILE', environment)

    def test_package_integrity_detects_changes_and_path_escape(self):
        path = self.root / 'source.txt'
        path.write_text('content')
        manifest = self.root / 'NONVIRAL_PACKAGE.json'
        manifest.write_text(json.dumps(dict(files={'source.txt': probe.sha(path)})))
        probe.verify_package(self.root)
        path.write_text('changed')
        with self.assertRaisesRegex(ValueError, 'missing/changed'):
            probe.verify_package(self.root)
        manifest.write_text(json.dumps(dict(files={'../outside': 'not-a-digest'})))
        with self.assertRaisesRegex(ValueError, 'missing/changed'):
            probe.verify_package(self.root)

    def test_review_includes_failure_logs_but_not_environments_or_large_sources(self):
        output = self.root / 'results'
        output.mkdir()
        (output / 'report.json').write_text('{"status":"FAIL"}')
        for name in ('native/first_setup.console.log', 'native/work/aa/bb/.command.err',
                     'native/conda/env-abc/large', 'live/sources/downloads/a/protein.faa.gz'):
            path = output / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('fixture')
        archive = probe.make_review(output)
        import tarfile
        with tarfile.open(archive) as handle:
            names = handle.getnames()
        self.assertIn('validation_results/native/work/aa/bb/.command.err', names)
        self.assertNotIn('validation_results/native/conda/env-abc/large', names)
        self.assertFalse(any(name.endswith('protein.faa.gz') for name in names))


if __name__ == '__main__':
    unittest.main()
