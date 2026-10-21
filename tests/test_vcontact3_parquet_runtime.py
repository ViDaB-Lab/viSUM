"""Parquet readiness tests; no vConTACT3 execution or reference downloads."""
import contextlib
import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
import check_vcontact3_parquet as check

AVAILABLE = all(importlib.util.find_spec(name) for name in ('pandas', 'pyarrow'))
REQUIRED = os.environ.get('VISUM_REQUIRE_PARQUET_RUNTIME_TEST') == '1'


class Vcontact3ParquetRuntimeTests(unittest.TestCase):
    def test_missing_dependency_fails_with_actionable_message(self):
        with patch.object(check.importlib, 'import_module', side_effect=ImportError('No module named pyarrow')):
            with self.assertRaisesRegex(RuntimeError, 'selected vConTACT3 environment'):
                check.check_parquet_runtime()

    def test_cli_returns_failure_instead_of_silently_skipping(self):
        error = io.StringIO()
        with patch.object(check, 'check_parquet_runtime', side_effect=RuntimeError('missing reader')):
            with contextlib.redirect_stderr(error):
                self.assertEqual(check.main(), 1)
        self.assertIn('Parquet runtime check failed', error.getvalue())

    def test_cli_rejects_broken_backend_not_just_missing_import(self):
        error = io.StringIO()
        with patch.object(check, 'check_parquet_runtime', side_effect=ValueError('binary incompatibility')):
            with contextlib.redirect_stderr(error):
                self.assertEqual(check.main(), 1)
        self.assertIn('binary incompatibility', error.getvalue())

    @unittest.skipUnless(AVAILABLE or REQUIRED, 'Requires pandas and pyarrow; set VISUM_REQUIRE_PARQUET_RUNTIME_TEST=1 to require them')
    def test_real_parquet_roundtrip_uses_default_pandas_reader(self):
        versions = check.check_parquet_runtime()
        self.assertTrue(versions['pandas'])
        self.assertTrue(versions['pyarrow'])

    @unittest.skipUnless(AVAILABLE or REQUIRED, 'Requires pandas and pyarrow')
    def test_real_cli_roundtrip(self):
        result = subprocess.run([sys.executable, '-B', str(ROOT / 'bin/check_vcontact3_parquet.py')],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('PASS: vConTACT3 Parquet round-trip', result.stdout)


if __name__ == '__main__':
    unittest.main(verbosity=2)
