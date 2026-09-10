"""Linux routing test with fake prepare modules; never downloads databases.

Run with Nextflow on PATH. No analysis process is mocked: scheduling one is a
test failure. Real installers and database builds require a separate smoke test.
"""
import csv
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')


@unittest.skipUnless(os.name == 'posix' and NEXTFLOW, 'Requires Linux/Unix with Nextflow on PATH')
class SetupNextflowIntegrationTests(unittest.TestCase):
    def test_all_prepares_and_no_analysis_without_fasta(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            shutil.copy(ROOT / 'visum_nextflow.nf', project)
            shutil.copy(ROOT / 'visum.config', project)
            for directory in ('workflows', 'modules', 'envs', 'bin'):
                shutil.copytree(ROOT / directory, project / directory)
            setup = (project / 'workflows/database_setup.nf').read_text()
            expected = re.findall(r'include \{ (\w+) \}', setup)
            for process, module in re.findall(r"include \{ (\w+) \} from '../modules/local/(\w+)'", setup):
                module_path = project / f'modules/local/{module}.nf'
                original = module_path.read_text()
                inputs = original.split('    input:', 1)[1].split('    output:', 1)[0]
                outputs = original.split('    output:', 1)[1]
                outputs = re.split(r'\n    (?:script|shell):', outputs, maxsplit=1)[0]
                outputs = outputs.split('    // A shell block', 1)[0]
                names = re.findall(r"path\('([^']+)'\)", outputs)
                commands = '\n'.join(f"printf 'fixture\\n' > {name}" for name in names)
                module_path.write_text(f'process {process} {{\n    input:{inputs}\n    output:{outputs}\n    script:\n    """\n{commands}\n    """\n}}\n')
            (project / 'nonviral').mkdir()
            config = project / 'test.config'
            config.write_text('''
params.run_vicat = true
params.run_vitap = true
params.run_vcontact3 = true
params.run_tesorter = false
params.vicat_db = "${projectDir}/viral"
params.vicat_nonviral_db = "${projectDir}/nonviral"
conda.enabled = false
dag.enabled = false
report.enabled = false
timeline.enabled = false
trace.enabled = true
trace.file = "${projectDir}/setup-trace.tsv"
process.executor = 'local'
process.errorStrategy = 'terminate'
''')
            result = subprocess.run([NEXTFLOW, 'run', str(project / 'visum_nextflow.nf'),
                '-c', str(project / 'visum.config'), '-c', str(config), '--setup',
                '--outdir', str(project / 'results'), '-ansi-log', 'false'],
                cwd=project, capture_output=True, text=True, timeout=180)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (project / 'setup-trace.tsv').open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            observed = {row['name'].split(' (')[0].split(':')[-1] for row in rows}
            self.assertEqual(observed, set(expected))
            self.assertEqual(len(rows), len(expected))
            self.assertTrue(all(row['status'] == 'COMPLETED' for row in rows))


if __name__ == '__main__':
    unittest.main()
