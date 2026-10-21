"""Small deterministic archive/integrity checks; no pipeline execution."""
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import package_alpha_candidate as freeze

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from verify_alpha_candidate import verify


class CandidatePackageTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.project = self.root / 'source'
        self.project.mkdir()
        for name in freeze.TOP:
            (self.project / name).write_bytes(b'fixture\r\n')
        for folder in freeze.TREES:
            (self.project / folder).mkdir()
        (self.project / 'bin/science.py').write_text('VALUE = 1\n')
        (self.project / 'bin/tool.sh').write_bytes(b'#!/bin/sh\r\nexit 0\r\n')
        (self.project / 'envs/tool.yml').write_text('name: fixture\n')
        self.validation = self.root / 'validation.json'
        self.validation.write_text(json.dumps(dict(status='PASS', scientific_sources_unchanged={
            'bin/science.py': hashlib.sha256((self.project / 'bin/science.py').read_bytes()).hexdigest()})))

    def package(self, name):
        with patch.object(freeze.subprocess, 'check_output', return_value='a' * 40 + '\n'):
            return freeze.package(self.project, self.root / name, self.validation)

    def test_reproducible_rc3_and_extracted_integrity(self):
        first, second = self.package('first.tar.gz'), self.package('second.tar.gz')
        self.assertEqual(first['sha256'], second['sha256'])
        self.assertEqual(first['release_id'], 'visum-alpha-rc3-20261020')
        with tarfile.open(first['archive']) as archive:
            names = archive.getnames()
            self.assertTrue(all(name.startswith(freeze.RELEASE + '/') for name in names))
            self.assertEqual(archive.getmember(freeze.RELEASE + '/bin/tool.sh').mode, 0o755)
            self.assertNotIn(b'\r', archive.extractfile(freeze.RELEASE + '/bin/tool.sh').read())
            archive.extractall(self.root / 'extracted', filter='data')
        extracted = self.root / 'extracted' / freeze.RELEASE
        self.assertEqual(verify(extracted)['status'], 'PASS')
        (extracted / 'bin/science.py').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'missing/changed'):
            verify(extracted)

    def test_no_overwrite_of_archive_or_checksum(self):
        self.package('existing.tar.gz')
        with self.assertRaises(FileExistsError):
            self.package('existing.tar.gz')
        (self.root / 'sidecar.tar.gz.sha256').write_text('existing')
        with self.assertRaises(FileExistsError):
            self.package('sidecar.tar.gz')

    def test_scientific_change_or_nonpassing_validation_blocks_freeze(self):
        (self.project / 'bin/science.py').write_text('VALUE = 2\n')
        with self.assertRaisesRegex(ValueError, 'scientific source changed'):
            self.package('changed.tar.gz')
        self.validation.write_text('{"status":"FAIL"}')
        with self.assertRaisesRegex(ValueError, 'passing validation'):
            self.package('failed.tar.gz')


if __name__ == '__main__':
    unittest.main()
