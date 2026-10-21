"""Tiny environment/socket regressions; never execute TEsorter or a dataset."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / 'modules/local/run_tesorter.nf'
BASH = shutil.which('bash') if os.name == 'posix' else next(
    (str(p) for p in (Path('C:/Program Files/Git/bin/bash.exe'),
                     Path('C:/Program Files/Git/usr/bin/bash.exe')) if p.is_file()), None)


def bootstrap():
    text = MODULE.read_text(encoding='utf-8')
    match = re.search(r"beforeScript\s+'''(.*?)'''", text, re.S)
    if not match:
        raise AssertionError('TEsorter short-temp bootstrap is missing')
    return textwrap.dedent(match.group(1)).strip()


@unittest.skipUnless(BASH, 'Requires Bash')
class TEsorterTemporaryRootTests(unittest.TestCase):
    def test_bootstrap_is_valid_bash(self):
        result = subprocess.run([BASH, '-n'], input=bootstrap(), text=True,
                                capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_long_inherited_roots_are_overridden_without_changing_parent(self):
        environment = dict(os.environ)
        environment.update({key: '/very-long-inherited-root/' + 'x' * 160
                            for key in ('TMPDIR', 'TEMP', 'TMP')})
        before = environment.copy()
        command = bootstrap() + '\nprintf "%s\\n" "$TMPDIR" "$TEMP" "$TMP"'
        result = subprocess.run([BASH, '-c', command], env=environment,
                                text=True, capture_output=True, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.splitlines(), ['/tmp', '/tmp', '/tmp'])
        self.assertEqual(environment, before)

    @unittest.skipUnless(sys.platform == 'linux', 'Requires Linux AF_UNIX/forkserver')
    def test_forkserver_pool_works_with_long_inherited_temp_root(self):
        # Explicit forkserver also exercises this failure mode on Python 3.11+
        # instead of depending on the interpreter's platform default method.
        probe = '''
import multiprocessing as mp
import tempfile
assert tempfile.gettempdir() == '/tmp', tempfile.gettempdir()
with mp.get_context('forkserver').Pool(2) as pool:
    assert pool.map(abs, [-1, 2, -3]) == [1, 2, 3]
print('PASS short-temp forkserver pool')
'''
        with tempfile.TemporaryDirectory(prefix='visum-temp-regression-') as scratch:
            long_root = Path(scratch) / ('long-' + 'x' * 120)
            long_root.mkdir()
            environment = dict(os.environ)
            environment.update({key: str(long_root) for key in ('TMPDIR', 'TEMP', 'TMP')})
            # Prove the fixture itself triggers the OS socket-path restriction.
            negative = subprocess.run(
                [sys.executable, '-B', '-c',
                 "import os,socket; s=socket.socket(socket.AF_UNIX); s.bind(os.path.join(os.environ['TMPDIR'],'socket'))"],
                env=environment, text=True, capture_output=True, timeout=15)
            self.assertNotEqual(negative.returncode, 0)
            self.assertIn('AF_UNIX path too long', negative.stderr)
            positive = subprocess.run(
                [BASH, '-c', bootstrap() + '\nexec "$@"', 'tesorter-temp-test',
                 sys.executable, '-B', '-c', probe], env=environment,
                text=True, capture_output=True, timeout=30)
            self.assertEqual(positive.returncode, 0, positive.stdout + positive.stderr)
            self.assertIn('PASS short-temp forkserver pool', positive.stdout)
            self.assertEqual(list(long_root.iterdir()), [])


if __name__ == '__main__':
    unittest.main(verbosity=2)
