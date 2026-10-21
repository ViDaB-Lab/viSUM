"""User-facing help/docs and release-integrity contracts, no network or tools."""
import hashlib
import json
from pathlib import Path
import re
import sys
import tarfile
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'tests'),str(ROOT/'bin')]
import package_beta_candidate as freeze
from verify_alpha_candidate import verify

USER_DOCS=('docs/README.md','docs/usage.md','docs/database-setup.md',
    'docs/vicat-database-preparation.md','docs/vitap-database-lifecycle.md',
    'docs/resources.md','docs/interpretation.md','docs/troubleshooting.md',
    'docs/checkv-evidence-context.md','docs/provirus-boundary-policy.md',
    'docs/reporting-bundle.md','docs/te-candidate-routing.md')

class UserDocumentationTests(unittest.TestCase):
    def test_version_and_identity_are_consistent(self):
        config=(ROOT/'visum.config').read_text()
        workflow=(ROOT/'visum_nextflow.nf').read_text()
        citation=(ROOT/'CITATION.cff').read_text()
        self.assertIn("version     = '0.9.0-beta.2'",config)
        self.assertIn('viSUM 0.9.0-beta.2',workflow)
        self.assertIn('version: 0.9.0-beta.2',citation)
        self.assertNotIn('<Your Name>',config)
        self.assertNotRegex(citation,r'(?m)^doi:')

    def test_readme_preserves_motivation_license_citation_and_changes(self):
        text=(ROOT/'README.md').read_text()
        for phrase in ('## Motivation',"## What's changed since beta.1",'help my students',
            '10.5281/zenodo.22774539','0.9.0-beta.1','CITATION.cff','LICENSE',
            'not a beta.2 DOI','local review candidate'):
            if phrase=='not a beta.2 DOI':
                self.assertIn('No beta.2 DOI has been assigned',text)
            else:self.assertIn(phrase,text)
        self.assertEqual(hashlib.sha256((ROOT/'LICENSE').read_bytes()).hexdigest(),
            'd2354af2e50346e7dc2e760056480c70140d9127655bceb1fba062631b5d099a')

    def test_user_guides_do_not_require_benchmark_or_private_test_steps(self):
        forbidden=('benchmark','CELL_NEG_','BMARK_','student replay','native Linux v2',
            'pre-alpha','isolated alpha','validation gate','NOT TESTED','pending manifest')
        for name in USER_DOCS:
            text=(ROOT/name).read_text()
            for phrase in forbidden:self.assertNotIn(phrase.lower(),text.lower(),(name,phrase))

    def test_defaults_are_usable_without_location_and_report_overrides(self):
        readme=(ROOT/'README.md').read_text()
        self.assertIn('nextflow run visum_nextflow.nf --setup\n',readme)
        self.assertIn('nextflow run visum_nextflow.nf --setup --run_vicat true -profile vicat_build',readme)
        for flag in ('--dbdir','--tooldir','-with-trace','-with-report'):
            self.assertNotIn(flag,readme)

    def test_user_docs_link_targets_exist_and_code_fences_balance(self):
        for name in ('README.md','CHANGELOG.md',*USER_DOCS,'docs/development/README.md',
                     'docs/development/beta2-release-review.md'):
            path=ROOT/name;text=path.read_text()
            self.assertEqual(sum(line.startswith('```') for line in text.splitlines())%2,0,name)
            for link in re.findall(r'\]\(([^)]+)\)',text):
                target=link.split('#',1)[0]
                if not target or '://' in target:continue
                self.assertTrue((path.parent/target).exists(),(name,link))

    def test_user_guides_describe_reuse_scope_and_safety(self):
        usage=(ROOT/'docs/usage.md').read_text()
        self.assertIn('nextflow log',usage)
        self.assertIn('.nextflow/',usage)
        self.assertIn('work/',usage)
        self.assertIn('cross-version reuse',usage)
        resource=(ROOT/'docs/resources.md').read_text()
        self.assertIn('not an established minimum',resource)
        self.assertIn('YOUR_PARTITION',resource)
        self.assertIn('YOUR_ACCOUNT',resource)

    def test_help_matches_parameter_names_and_current_source_approval(self):
        workflow=(ROOT/'visum_nextflow.nf').read_text()
        config=(ROOT/'visum.config').read_text()
        help_text=workflow.split('def visumHelp() {',1)[1].split('include { DATABASE_SETUP }',1)[0]
        for name in ('setup','harmonizer_max_rank','harmonizer_cpus','harmonizer_memory',
                     'vicat_metavr_auto_download','vicat_nonviral_auto_download'):
            self.assertIn('--'+name,help_text)
            self.assertRegex(config,rf'(?m)^\s*{name}\s*=')
        self.assertNotIn('pending until full comparison',help_text)
        self.assertIn('same --run_* flags',help_text)

    def test_specific_student_recovery_and_legacy_files_are_absent(self):
        for name in ('tests/larvae_recovery.py','tests/package_larvae_recovery.py',
                     'tests/test_larvae_recovery.py','bin/fixgenomadv0.2.py',
                     'bin/fixct3v0.2.py','bin/fixvitapv0.4.py'):
            self.assertFalse((ROOT/name).exists(),name)

    def test_documented_workflow_flags_exist_in_config(self):
        names=set(re.findall(r'(?m)^\s*([a-z][a-z0-9_]*)\s*=',(ROOT/'visum.config').read_text()))
        for name in ('README.md',*USER_DOCS):
            for flag in re.findall(r'--([a-z][a-z0-9_]*)',(ROOT/name).read_text()):
                # Git installation and SLURM site-config examples are not workflow flags.
                if flag in ('account','branch') or flag.endswith('_'):continue
                self.assertIn(flag,names,(name,flag))

class BetaPackageTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.base=Path(self.temp.name);self.project=self.base/'project';self.project.mkdir()
        for name in freeze.TOP:(self.project/name).write_bytes(b'fixture\r\n')
        for name in freeze.TREES:(self.project/name).mkdir()
        (self.project/'bin/science.py').write_text('VALUE=1\n')
        (self.project/'bin/tool.sh').write_bytes(b'#!/bin/sh\r\nexit 0\r\n')
        self.validation=self.base/'validation.json'
        self.validation.write_text(json.dumps(dict(status='PASS',scientific_sources_unchanged={
            'bin/science.py':hashlib.sha256((self.project/'bin/science.py').read_bytes()).hexdigest()})))

    def test_determinism_license_citation_and_extracted_integrity(self):
        first=freeze.package(self.project,self.base/'first.tar.gz',self.validation)
        second=freeze.package(self.project,self.base/'second.tar.gz',self.validation)
        self.assertEqual(first['sha256'],second['sha256'])
        with tarfile.open(first['archive']) as archive:
            for name in ('LICENSE','CITATION.cff','CHANGELOG.md'):
                self.assertIn(freeze.RELEASE+'/'+name,archive.getnames())
            self.assertEqual(archive.getmember(freeze.RELEASE+'/bin/tool.sh').mode,0o755)
            self.assertNotIn(b'\r',archive.extractfile(freeze.RELEASE+'/bin/tool.sh').read())
            archive.extractall(self.base/'unpacked',filter='data')
        target=self.base/'unpacked'/freeze.RELEASE
        self.assertEqual(verify(target)['status'],'PASS')
        (target/'LICENSE').write_text('changed')
        with self.assertRaises(ValueError):verify(target)

    def test_no_overwrite_or_unvalidated_source_changes(self):
        output=self.base/'first.tar.gz'
        freeze.package(self.project,output,self.validation)
        with self.assertRaises(FileExistsError):freeze.package(self.project,output,self.validation)
        (self.project/'bin/science.py').write_text('VALUE=2\n')
        with self.assertRaisesRegex(ValueError,'scientific source changed'):
            freeze.package(self.project,self.base/'changed.tar.gz',self.validation)
        self.validation.write_text('{"status":"FAIL"}')
        with self.assertRaises(ValueError):freeze.package(self.project,self.base/'failed.tar.gz',self.validation)

if __name__=='__main__':unittest.main(verbosity=2)
