"""Execute the selection guards from the actual Nextflow shell template."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE = (ROOT / 'modules/local/vitap_database.nf').read_text(encoding='utf-8')
BASH = (str(Path('C:/Program Files/Git/bin/bash.exe')) if os.name == 'nt'
        else shutil.which('bash'))


def section(start, end):
    return MODULE.split(start, 1)[1].split(end, 1)[0].replace('\\$', '$')


@unittest.skipUnless(BASH and Path(BASH).exists(), 'Bash required')
class VitapDatabaseSelectionTests(unittest.TestCase):
    def run_shell(self, code):
        with tempfile.TemporaryDirectory() as work:
            env = os.environ.copy()
            if os.name == 'nt':
                env['MSYS'] = 'winsymlinks:nativestrict'
            result = subprocess.run([BASH, '-c', 'export PATH=/usr/bin:/bin:$PATH\nset -eu\n' + code], cwd=work,
                                    capture_output=True, text=True, env=env)
        return result

    def require_symlinks(self):
        result = self.run_shell('mkdir target\nln -s target link\n[[ -L link ]]')
        if result.returncode:
            self.skipTest('Native symlink creation unavailable: ' + result.stderr.strip())

    def test_unmodified_template_has_valid_shell_syntax(self):
        script = MODULE.split('    """', 1)[1].rsplit('    """', 1)[0]
        import re
        script = re.sub(r'(?<!\\)\$\{[^}]+\}', 'placeholder', script)
        script = script.replace('\\$', '$')
        result = subprocess.run([BASH, '-n'], input=script, text=True,
                                capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_only_implicit_reuse_takes_current_shortcut(self):
        self.require_symlinks()
        guard = section(
            '    flock 9\n', '    if [[ -n "\\$VMR_SOURCE"')
        for update, source, label, reused in [
            ('false', '', '', True), ('false', 'requested.xlsx', '', False),
            ('false', '', 'requested-release', False), ('true', '', '', False),
        ]:
            with self.subTest(update=update, source=source, label=label):
                result = self.run_shell(f'''
DB_ROOT="$PWD/db"
mkdir -p "$DB_ROOT/release"
ln -s release "$DB_ROOT/current"
UPDATE_DATABASE='{update}'
VMR_SOURCE='{source}'
REQUESTED_LABEL='{label}'
validate_database() {{ return 0; }}
emit_existing_managed_database() {{ echo REUSED; }}
{guard}
echo CONTINUE
''')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual('REUSED' in result.stdout, reused)
                self.assertEqual('CONTINUE' in result.stdout, not reused)

    def test_requested_release_requires_matching_recorded_digest(self):
        self.require_symlinks()
        guard = '    if validate_database "$TARGET_DATABASE";' + section(
            '    if validate_database "\\$TARGET_DATABASE";',
            '    if [[ -e "\\$TARGET_DATABASE"')
        for stored in ['requested-hash', 'different-hash', None]:
            with self.subTest(stored=stored):
                metadata = ('' if stored is None else
                            f"printf 'header\\nrelease\\tfile\\turl\\t{stored}\\n' > \"$TARGET_DATABASE/.visum_vmr_metadata.tsv\"")
                result = self.run_shell(f'''
DB_ROOT="$PWD/db"
RELEASE=new
TARGET_DATABASE="$DB_ROOT/DB_new"
mkdir -p "$TARGET_DATABASE" "$DB_ROOT/DB_old"
ln -s DB_old "$DB_ROOT/current"
trap 'echo CURRENT=$(readlink "$DB_ROOT/current")' EXIT
VMR_SHA=requested-hash
{metadata}
validate_database() {{ return 0; }}
emit_existing_managed_database() {{ echo REUSED; }}
{guard}
''')
                matches = stored == 'requested-hash'
                self.assertEqual(result.returncode == 0, matches, result.stderr)
                self.assertIn('CURRENT=DB_new' if matches else 'CURRENT=DB_old', result.stdout)

    def test_fingerprint_guard_without_symlink_privileges(self):
        guard = "        RECORDED_SHA=''" + section(
            "        RECORDED_SHA=''", '        ln -sfn "DB_\\$RELEASE"')
        for stored in ['requested-hash', 'different-hash', None]:
            with self.subTest(stored=stored):
                metadata = ('' if stored is None else
                            f"printf 'header\\nr\\tf\\tu\\t{stored}\\n' > \"$TARGET_DATABASE/.visum_vmr_metadata.tsv\"")
                result = self.run_shell(f'''
TARGET_DATABASE="$PWD/db"
mkdir -p "$TARGET_DATABASE"
VMR_SHA=requested-hash
{metadata}
{guard}
echo MATCH
''')
                self.assertEqual(result.returncode == 0, stored == 'requested-hash', result.stderr)

    def test_interrupted_work_requires_matching_digest(self):
        guard = '    BUILD_ROOT="$DB_ROOT/.build-$RELEASE"' + section(
            '    BUILD_ROOT="\\$DB_ROOT/.build-\\$RELEASE"',
            '    if [[ ! -e "\\$BUILD_ROOT/VMR_Genome"')
        for stored in ['requested-hash', 'different-hash', None, 'absent']:
            with self.subTest(stored=stored):
                setup = '' if stored == 'absent' else 'mkdir -p "$DB_ROOT/.build-new"'
                if stored not in (None, 'absent'):
                    setup += f'\nprintf "%s\\n" "{stored}" > "$DB_ROOT/.build-new/.visum_build_vmr.sha256"'
                result = self.run_shell(f'''
DB_ROOT="$PWD/db"
RELEASE=new
VMR_SHA=requested-hash
{setup}
{guard}
''')
                self.assertEqual(result.returncode == 0, stored in ('requested-hash', 'absent'), result.stderr)


if __name__ == '__main__':
    unittest.main()
