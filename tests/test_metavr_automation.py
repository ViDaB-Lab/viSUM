"""Approval integrity and actual private source-preparation protocol, tiny/offline."""
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'bin'),str(ROOT/'tests')]
import prepare_vicat_metavr_sources as sources
import test_metavr_source_patch_integration as fixtures
from run_metavr_automation_validation import validate


class ApprovalTests(unittest.TestCase):
    def test_packaged_approval_matches_native_evidence(self):
        manifest=sources.load_manifest(ROOT/'assets/vicat/metavr_portal_v1.json')
        evidence=json.loads((ROOT/'assets/vicat/metavr_native_validation_v1.json').read_text())
        self.assertEqual(evidence['protein']['fingerprint'],manifest['comparison']['proteins'])
        self.assertEqual(evidence['protein']['portal_sha256'],manifest['files']['proteins']['sha256'])
        self.assertEqual(evidence['protein']['original_sha256'],manifest['frozen']['protein_sha256'])
        self.assertEqual(evidence['metadata']['prepared_metadata_sha256'],manifest['comparison']['metadata']['prepared_metadata_sha256'])
        self.assertEqual(evidence['metadata']['records'],manifest['frozen']['uvigs'])
        self.assertEqual(evidence['metadata']['applied_rows'],1409)
        self.assertEqual(evidence['metadata']['non_votu_changed_rows'],0)
        self.assertFalse(evidence['database_build_performed'])

    def test_changed_evidence_rejected_before_acquisition(self):
        with tempfile.TemporaryDirectory() as temporary:
            target=Path(temporary)/'assets'
            shutil.copytree(ROOT/'assets/vicat',target)
            (target/'metavr_native_validation_v1.json').write_text('{}')
            with patch.object(sources,'urlopen',side_effect=AssertionError('No network')):
                with self.assertRaisesRegex(ValueError,'validation evidence'):
                    sources.prepare_sources(target/'metavr_portal_v1.json',Path(temporary)/'cache',
                        json.loads((target/'metavr_portal_v1.json').read_text())['frozen'])


class ProtocolTests(unittest.TestCase):
    def setUp(self):
        self.fixture=fixtures.IntegratedPatchTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.fixture.temporary.cleanup)
        self.f=self.fixture.fixture
        self.manifest=self.f.compare()
        expected=self.f.root/'expected.tsv.gz'
        sources.prepare_metadata(self.f.template,self.f.manifest,self.f.new_metadata,expected,self.f.root/'expected-work')
        self.manifest['comparison']['metadata']=dict(records=2,prepared_metadata_sha256=sources.digest(expected))
        self.project=self.f.root/'private/package'
        assets=self.project/'assets/vicat'
        assets.mkdir(parents=True)
        shutil.copytree(self.f.root/'bundle',assets/'bundle')
        sources.atomic_json(assets/'metavr_portal_v1.json',self.manifest)

    def run_protocol(self):
        arguments=(self.f.new_proteins,self.f.new_metadata,self.project.parent/'results',self.project,0)
        if os.name=='nt':
            # Windows unit hosts may lack symlink privilege. Copy ONLY tiny
            # synthetic fixtures here; the Linux runner always links real files.
            def tiny_copy(link,target):
                if Path(target).stat().st_size>4096: raise AssertionError('Tiny fixtures only')
                shutil.copyfile(target,link)
            with patch.object(Path,'symlink_to',tiny_copy): return validate(*arguments)
        return validate(*arguments)

    def test_actual_prepare_and_reuse_runner_with_synthetic_sources(self):
        report=self.run_protocol()
        self.assertEqual(report['status'],'PASS',report.get('traceback'))
        self.assertTrue(report['cached_hashes_and_timestamps_unchanged'])
        self.assertTrue(all(report['inputs_unchanged'].values()))
        self.assertFalse(any(report['forbidden_boundary_calls'].values()))

    def test_wrong_approved_prepared_checksum_fails(self):
        self.manifest['comparison']['metadata']['prepared_metadata_sha256']='0'*64
        sources.atomic_json(self.project/'assets/vicat/metavr_portal_v1.json',self.manifest)
        report=self.run_protocol()
        self.assertEqual(report['status'],'FAIL')
        self.assertIn('differs from approved',report['error'])

    def test_changed_protein_checksum_fails_without_network(self):
        self.manifest['files']['proteins']['sha256']='0'*64
        sources.atomic_json(self.project/'assets/vicat/metavr_portal_v1.json',self.manifest)
        report=self.run_protocol()
        self.assertEqual(report['status'],'FAIL')
        self.assertIn('Completed MetaVR source changed',report['error'])

    def test_existing_output_refused(self):
        output=self.project.parent/'results'
        output.mkdir()
        with self.assertRaises(FileExistsError):
            validate(self.f.new_proteins,self.f.new_metadata,output,self.project,0)


if __name__=='__main__': unittest.main(verbosity=2)
