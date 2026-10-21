"""Portable fail-closed checks for the tiny Linux resource verification harness."""
import copy
import hashlib
import json
from pathlib import Path
import re
import tempfile
import unittest

import verify_resource_defaults_linux as probe


class ResourceProbeTests(unittest.TestCase):
    def test_exact_uploaded_docker_warning_output_is_preserved_and_parsed(self):
        fixture = json.loads((probe.ROOT / 'tests/fixtures/resource_defaults/docker_warning_20261015.json').read_text())
        output = fixture['console_output']
        self.assertEqual(hashlib.sha256(output.encode('utf-8')).hexdigest(), fixture['sha256'])
        self.assertEqual(fixture['sha256'], 'b1619fe666428d9ed839613ae63b39a7a7e12201ef16077586f503696ece890d')
        with self.assertRaises(json.JSONDecodeError):
            json.loads(output)  # Reproduce the original harness failure.
        config, warnings = probe.parse_config_output(output)
        self.assertEqual(warnings, ["WARN: Unrecognized config option 'docker.pull'"])
        self.assertEqual(config, json.loads(output.split('\n', 1)[1]))
        self.assertTrue(all(config['params'][name] is None for name in probe.time_parameters()))
        self.assertEqual(config['params']['harmonizer_memory'], '32 GB')
        # Parsing must not suppress the real configuration problem.
        with self.assertRaisesRegex(AssertionError, 'unsupported docker.pull'):
            probe.validate_config(config, 'docker', probe.time_parameters())
        del config['docker']['pull']
        probe.validate_config(config, 'docker', probe.time_parameters())

    def test_parser_accepts_plain_json_and_records_only_explicit_warning_lines(self):
        config = self.config()
        config['literal'] = 'WARN: {this stays inside the JSON string}'
        plain = json.dumps(config, indent=2)
        self.assertEqual(probe.parse_config_output(plain), (config, []))
        output = ' \nWARN: first\r\n\nWARN: second\n' + plain + '\nWARN: third\n'
        self.assertEqual(probe.parse_config_output(output),
                         (config, ['WARN: first', 'WARN: second', 'WARN: third']))

    def test_parser_rejects_errors_malformed_json_extra_documents_and_unknown_noise(self):
        plain = json.dumps(self.config())
        bad = [prefix + '\n' + plain for prefix in ('ERROR: failed', 'INFO: unknown', 'Traceback:', 'noise')]
        bad += ['WARN: no JSON\n', 'WARN: warning\n{"params":', '[]', 'null',
                plain + '\nERROR: failed', plain + '\nnoise', plain + '\n{}']
        for output in bad:
            with self.subTest(output=output[:100]), self.assertRaises(ValueError):
                probe.parse_config_output(output)

    def test_rerun_selects_only_unfinished_checks_and_removes_unsupported_option(self):
        self.assertEqual(probe.UNFINISHED_PROFILES, ('docker', 'singularity'))
        self.assertEqual(probe.PROBE_CASES, ('defaults', 'parameters', 'site'))
        source = Path(probe.__file__).read_text()
        self.assertIn('for profile in UNFINISHED_PROFILES:', source)
        self.assertIn('for case in PROBE_CASES:', source)
        shell = (probe.ROOT / 'tests/run_resource_defaults_linux.sh').read_text()
        self.assertNotIn('tests/test_resource_workflow_contracts.py', shell)
        config = (probe.ROOT / 'visum.config').read_text()
        self.assertNotRegex(config, r'(?m)^\s*docker\.pull\s*=')
        self.assertIn('docker.enabled    = true', config)
        self.assertIn("docker.runOptions = '-u $(id -u):$(id -g)'", config)

    def test_generated_graph_copies_all_time_directives_but_no_analysis_bodies(self):
        graph, names = probe.probe_workflow()
        self.assertEqual(len(names), 27)
        self.assertEqual(set(names.values()), set(probe.time_parameters()))
        self.assertEqual(len(re.findall(r'(?m)^process RESOURCE_', graph)), 27)
        self.assertEqual(graph.count("printf '%s\\n'"), 27)
        for forbidden in ('include {', 'conda ', 'genomad end-to-end', 'diamond blastp',
                          'run_viharmony.py', 'wget ', 'curl ', '--setup', '-resume', 'sbatch'):
            self.assertNotIn(forbidden, graph)
        self.assertIn('memory params.harmonizer_memory', graph)
        self.assertIn('Math.min(params.harmonizer_cpus as int, params.max_cpus as int)', graph)

    def config(self):
        return dict(params=dict.fromkeys(probe.time_parameters()) | dict(
            harmonizer_memory='32 GB', harmonizer_cpus=1, max_memory='64 GB'),
            process=dict(time='{ params.utility_time }'))

    def test_default_config_and_all_profile_budgets(self):
        for profile in probe.PROFILES:
            value = self.config()
            value['params']['max_memory'] = {'local_full': '128 GB', 'vicat_build': '400 GB'}.get(profile, '64 GB')
            probe.validate_config(value, profile, probe.time_parameters())

    def test_missing_null_finite_profile_time_wrong_memory_and_workers_are_rejected(self):
        variants = []
        missing = self.config()
        del missing['params']['genomad_time']
        variants.append(missing)
        for key, value in (('genomad_time', '12h'), ('harmonizer_memory', '4 GB'),
                           ('harmonizer_cpus', 8), ('max_memory', '48 GB')):
            config = self.config()
            config['params'][key] = value
            variants.append(config)
        config = self.config()
        config['process']['time'] = '12h'
        variants.append(config)
        for value in variants:
            with self.subTest(config=value), self.assertRaises(AssertionError):
                probe.validate_config(value, 'default', probe.time_parameters())

    def observations(self, directory, case):
        names = {'RESOURCE_UTILITY': 'utility_time', 'RESOURCE_VIHARMONY': 'harmonizer_time'}
        values = {
            'RESOURCE_UTILITY': dict(probe='RESOURCE_UTILITY', cpus=1, memory_bytes=64 * 1024 ** 2,
                                    time_ms={'defaults': None, 'parameters': 10800000, 'site': 18000000}[case]),
            'RESOURCE_VIHARMONY': dict(probe='RESOURCE_VIHARMONY', cpus=1 if case == 'defaults' else 2,
                                      memory_bytes=(32 if case == 'defaults' else 48) * probe.GIB,
                                      time_ms={'defaults': None, 'parameters': 10800000, 'site': 25200000}[case]),
        }
        for name, value in values.items():
            probe.write(directory / (name + '.json'), json.dumps(value))
        return names, values

    def test_all_three_observation_cases_and_corrupt_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            for case in ('defaults', 'parameters', 'site'):
                directory = Path(temporary) / case
                names, values = self.observations(directory, case)
                self.assertEqual(probe.validate_observations(directory, names, case), values)
                for key, wrong in (('cpus', 99), ('memory_bytes', 4 * probe.GIB), ('time_ms', -1), ('probe', 'wrong')):
                    bad = copy.deepcopy(values['RESOURCE_VIHARMONY'])
                    bad[key] = wrong
                    probe.write(directory / 'RESOURCE_VIHARMONY.json', json.dumps(bad))
                    with self.subTest(case=case, field=key), self.assertRaises(AssertionError):
                        probe.validate_observations(directory, names, case)
                probe.write(directory / 'RESOURCE_VIHARMONY.json', json.dumps(values['RESOURCE_VIHARMONY']))
                (directory / 'RESOURCE_UTILITY.json').unlink()
                with self.assertRaises(AssertionError):
                    probe.validate_observations(directory, names, case)

    def test_runner_has_private_runtime_and_never_launches_production_workflow(self):
        shell = (probe.ROOT / 'tests/run_resource_defaults_linux.sh').read_text()
        source = Path(probe.__file__).read_text()
        self.assertIn('mktemp -d "$SNAPSHOT_DIR/../resource-defaults-run.', shell)
        self.assertIn('export NXF_HOME="$RUN_DIR/nextflow-home"', shell)
        self.assertIn('NXF_OFFLINE=true', shell)
        self.assertIn('NXF_SYNTAX_PARSER=v2', shell)
        self.assertIn("-XX:ActiveProcessorCount=2", shell)
        self.assertIn("executor = 'local'", source)
        self.assertNotIn('visum_nextflow.nf', shell)
        self.assertNotIn('visum_nextflow.nf', source)
        for forbidden in ('pip install', 'conda install', 'mamba install', 'git pull', 'sbatch ', 'pkill '):
            self.assertNotIn(forbidden, shell)


if __name__ == '__main__':
    unittest.main(verbosity=2)
