"""Harness contracts only: these tests do not simulate proof of Nextflow caching."""
import copy
import gzip
import json
import shutil
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_rank_cap_linux_validation as runner


class RankCapHarnessTests(unittest.TestCase):
    def command(self, rank=None, session=None):
        return runner.native_command(Path('/nextflow'), Path('/private'), Path('/private/test.config'),
            Path('/private/params.json'), 'first_family' if session is None else 'resume', rank, session,
            project=Path('/package'))

    def test_default_family_is_not_a_test_override(self):
        command = self.command()
        self.assertNotIn('--harmonizer_max_rank', command)
        self.assertNotIn('-resume', command)
        self.assertLess(command.index('-C'), command.index('run'))
        self.assertLess(command.index('-log'), command.index('run'))
        self.assertEqual(command[command.index('-work-dir') + 1], str(Path('/private/work')))

    def test_resume_always_has_the_original_explicit_session(self):
        session = '529b89df-9cd1-44d3-98ba-efeda26a9b3d'
        command = self.command('genus', session)
        self.assertEqual(command[command.index('-resume') + 1], session)
        self.assertEqual(command[command.index('--harmonizer_max_rank') + 1], 'genus')
        with self.assertRaises(ValueError):
            self.command('genus', 'latest')

    def rows(self, stage):
        rows = []
        for name in ('STAGE_SAVED_INPUTS', 'VIHARMONY'):
            for sample in ('DNA', 'RNA'):
                state = 'COMPLETED' if stage == 'first_family' or (name == 'VIHARMONY' and stage in {'genus','species'}) else 'CACHED'
                rows.append(dict(name=f'{name} ({sample})', status=state))
        return rows

    def test_all_five_expected_cache_states_are_checked(self):
        for stage, _, _ in runner.STAGES:
            with self.subTest(stage=stage):
                runner.verify_trace(self.rows(stage), stage, ('DNA','RNA'))

    def test_cache_miss_failure_duplicate_missing_or_search_task_is_rejected(self):
        cases = []
        rows = self.rows('unchanged_resume'); rows[0]['status'] = 'COMPLETED'; cases.append(('unchanged_resume',rows))
        rows = self.rows('genus'); rows[-1]['status'] = 'CACHED'; cases.append(('genus',rows))
        rows = self.rows('first_family'); rows[-1]['status'] = 'FAILED'; cases.append(('first_family',rows))
        cases.extend(('first_family', rows) for rows in (self.rows('first_family')[:-1],
                      self.rows('first_family') + [self.rows('first_family')[0]],
                      self.rows('first_family') + [dict(name='RUN_GENOMAD (DNA)',status='COMPLETED')]))
        for stage, rows in cases:
            with self.subTest(rows=rows), self.assertRaises(AssertionError):
                runner.verify_trace(rows, stage, ('DNA','RNA'))

    def test_actual_task_command_rank_and_workers_are_not_assumed(self):
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / 'command.sh'
            script.write_text('python3 /package/bin/run_viharmony.py --max-rank family --workers 8\n')
            runner.verify_task_command(script, 'family', 8)
            with self.assertRaises(AssertionError): runner.verify_task_command(script, 'genus', 8)
            with self.assertRaises(AssertionError): runner.verify_task_command(script, 'family', 4)

    def test_audit_fingerprint_preserves_values_columns_and_multiplicity(self):
        with tempfile.TemporaryDirectory() as temporary:
            a, b = (Path(temporary) / n for n in ('a_audit.tsv.gz', 'b_audit.tsv.gz'))
            a.write_bytes(gzip.compress(b'a\tb\n1\t2\n3\t4\n',mtime=0))
            b.write_bytes(gzip.compress(b'b\ta\n4\t3\n2\t1\n',mtime=1))
            self.assertEqual(runner.fingerprint(a),runner.fingerprint(b))
            for value in (b'b\ta\n4\t3\n2\t1\n2\t1\n', b'b\ta\n4\t3\n9\t1\n', b'c\ta\n4\t3\n2\t1\n'):
                b.write_bytes(gzip.compress(value))
                self.assertNotEqual(runner.fingerprint(a),runner.fingerprint(b))

    def test_manifest_validation_does_not_accept_changed_sequence_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            fasta = folder / 'sample.fasta'; fasta.write_text('>a\nACGT\n')
            manifest = dict(outputs=dict(fasta=dict(file=fasta.name,sha256=runner.sha(fasta))))
            (folder / 'sample.harmonizer_manifest.json').write_text(json.dumps(manifest))
            expected = {'family': {'sample': runner.scientific_snapshot(folder,'sample')}}
            runner.verify_expected(folder,'sample','family',expected)
            fasta.write_text('>a\nACTT\n')
            with self.assertRaises(AssertionError): runner.verify_expected(folder,'sample','family',expected)

    def test_native_fixture_imports_only_real_harmonizer_and_saved_copy_process(self):
        project = Path(__file__).resolve().parents[1]
        path = project / 'native_validation.nf'
        if not path.is_file(): path = project / 'tests/integrated_reconciliation_native.nf'
        text = path.read_text()
        self.assertEqual(text.count('include {'),1)
        self.assertIn('include { VIHARMONY }',text)
        for forbidden in ('RUN_GENOMAD','RUN_VIRSORTER2','RUN_VICAT_DIAMOND','DATABASE_SETUP'):
            self.assertNotIn(forbidden,text)

    def bundle(self, directory, rank):
        directory.mkdir(parents=True,exist_ok=True)
        path=directory/'sample.taxonomy.tsv'
        path.write_text('sequence_id\ttaxonomy\na\t'+rank+'\n')
        manifest=dict(policy=dict(harmonizer_max_rank=rank),outputs=dict(taxonomy=dict(file=path.name,sha256=runner.sha(path))))
        (directory/'sample.harmonizer_manifest.json').write_text(json.dumps(manifest))

    def copy_bundle(self, work, published):
        # Unit-test fixture only, NOT the native publication mechanism.
        published.mkdir(exist_ok=True)
        for path in work.iterdir(): shutil.copyfile(path,published/path.name)

    def test_internally_valid_old_published_bundle_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            work,published=(Path(temporary)/n for n in ('work','published'))
            self.bundle(work,'family'); self.bundle(published,'species')
            runner.scientific_snapshot(published,'sample')  # Old files are internally valid.
            with self.assertRaisesRegex(AssertionError,'stale'):
                runner.verify_publication(work,published,'sample')

    def test_new_manifest_with_stale_table_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            work,published=(Path(temporary)/n for n in ('work','published'))
            self.bundle(work,'family'); self.bundle(published,'species')
            shutil.copyfile(work/'sample.harmonizer_manifest.json',published/'sample.harmonizer_manifest.json')
            with self.assertRaisesRegex(AssertionError,'checksum'):
                runner.verify_publication(work,published,'sample')

    def test_cached_family_restoration_requires_refresh_of_published_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            family,species,published=(Path(temporary)/n for n in ('family','species','published'))
            self.bundle(family,'family'); self.bundle(species,'species')
            self.copy_bundle(family,published)
            self.assertEqual(runner.verify_publication(family,published,'sample')['published_files_verified'],2)
            self.copy_bundle(species,published)
            runner.verify_publication(species,published,'sample')
            with self.assertRaisesRegex(AssertionError,'stale'):
                runner.verify_publication(family,published,'sample')
            self.copy_bundle(family,published)
            runner.verify_publication(family,published,'sample')

    def test_verifier_never_repairs_or_changes_published_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            work,published=(Path(temporary)/n for n in ('work','published'))
            self.bundle(work,'family'); self.bundle(published,'species')
            before={p.name:runner.sha(p) for p in published.iterdir()}
            with self.assertRaises(AssertionError): runner.verify_publication(work,published,'sample')
            self.assertEqual(before,{p.name:runner.sha(p) for p in published.iterdir()})

    def test_real_module_explicitly_refreshes_and_fails_on_publish_error(self):
        project=Path(__file__).resolve().parents[1]
        source=(project/'modules/local/viharmony.nf').read_text()
        self.assertIn("mode: 'copy', overwrite: true, failOnError: true",source)


if __name__ == '__main__':
    unittest.main()
