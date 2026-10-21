import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULES = ROOT / "modules" / "local"


class ResourceWorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = (ROOT / "visum.config").read_text(encoding="utf-8")
        cls.workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        cls.launcher = (ROOT / "visum").read_text(encoding="utf-8")

    def read_module(self, name: str) -> str:
        return (MODULES / name).read_text(encoding="utf-8")

    def test_local_executor_uses_native_budgets_without_launcher_state(self) -> None:
        self.assertIn("executor.cpus = params.max_cpus as int", self.config)
        self.assertIn("executor.memory = params.max_memory", self.config)
        self.assertNotIn("taskset", self.launcher)
        self.assertNotIn("NXF_OPTS", self.launcher)
        self.assertNotIn("VISUM_LOCAL_CPU_CAP", self.workflow)
        self.assertIn('exec "$nextflow_bin" run "$ROOT_DIR/visum_nextflow.nf" "$@"', self.launcher)

    def test_utility_processes_inherit_one_cpu(self) -> None:
        process_block = re.search(
            r"process\s*\{(?P<body>.*?)\n\}", self.config, re.DOTALL
        )
        self.assertIsNotNone(process_block)
        self.assertRegex(process_block.group("body"), r"\bcpus\s*=\s*1\b")
        self.assertNotIn("params.threads", process_block.group("body"))

    def test_all_twenty_wall_time_defaults_are_null(self) -> None:
        expected = {
            'utility_time', 'genomad_time', 'checkv_time', 'tesorter_time',
            'deep6_time', 'deepmicroclass2_time', 'virsorter2_time', 'ct3_time',
            'vitap_time', 'vitap_standardizer_time', 'vicat_nonviral_build_time',
            'vicat_build_time', 'vicat_orf_time', 'vicat_time',
            'vicat_standardizer_time', 'gianthunter_time', 'virbot_time',
            'vcontact3_time', 'vcontact3_standardizer_time', 'harmonizer_time',
        }
        defaults = dict(re.findall(r'(?m)^\s*(\w+_time)\s*=\s*([^\n]+)', self.config))
        self.assertEqual(set(defaults), expected)
        self.assertTrue(all(value.strip() == 'null' for value in defaults.values()), defaults)

    def test_profiles_cannot_reintroduce_shipped_time_limits(self) -> None:
        declarations = re.findall(r'(?m)^\s*(?:process\.)?time\s*=\s*([^\n]+)', self.config)
        self.assertEqual(declarations, ['{ params.utility_time }'])

    def test_every_module_time_directive_uses_an_overridable_parameter(self) -> None:
        defaults = set(re.findall(r'(?m)^\s*(\w+_time)\s*=', self.config))
        declarations = []
        for path in MODULES.glob('*.nf'):
            for directive in re.findall(r'(?m)^\s*time\s+([^\n]+)', path.read_text(encoding='utf-8')):
                with self.subTest(module=path.name, directive=directive):
                    match = re.fullmatch(r'(?:params\.(\w+)|\{\s*params\.(\w+)\s*\})', directive.strip())
                    self.assertIsNotNone(match, 'A module bypasses the configurable unlimited default')
                    self.assertIn(match.group(1) or match.group(2), defaults)
                    declarations.append(path.name)
        self.assertGreaterEqual(len(declarations), 25)

    def test_vicat_setup_validation_and_build_share_the_explicit_time_override(self) -> None:
        for name, parameter in (
            ('vicat_database.nf', 'vicat_build_time'),
            ('vicat_nonviral_database.nf', 'vicat_nonviral_build_time'),
        ):
            with self.subTest(module=name):
                self.assertIn('time { params.' + parameter + ' }', self.read_module(name))

    def test_harmonizer_has_a_fixed_configurable_32_gb_request(self) -> None:
        self.assertRegex(self.config, r"harmonizer_memory\s*=\s*'32 GB'")
        self.assertRegex(self.config, r'harmonizer_cpus\s*=\s*1\b')
        module = self.read_module('viharmony.nf')
        self.assertIn('memory params.harmonizer_memory', module)
        self.assertIn('time params.harmonizer_time', module)
        self.assertIn('Math.min(params.harmonizer_cpus as int, params.max_cpus as int)', module)
        self.assertRegex(self.config, r"max_memory\s*=\s*'64 GB'")

    def test_each_discovery_program_uses_its_own_cpu_parameter(self) -> None:
        expected = {
            "run_genomad.nf": "params.genomad_cpus",
            "run_virsorter2.nf": "params.virsorter2_cpus",
            "run_cenotetaker3.nf": "params.ct3_cpus",
            "run_deep6.nf": "params.deep6_cpus",
            "run_deepmicroclass2.nf": "params.deepmicroclass2_cpus",
            "run_virbot.nf": "params.virbot_cpus",
            "run_gianthunter.nf": "params.gianthunter_cpus",
            "predict_vicat_orfs.nf": "params.vicat_orf_cpus",
            "run_vicat_diamond.nf": "params.vicat_cpus",
            "prepare_vicat_reference_subset.nf": "params.vicat_standardizer_cpus",
            "prepare_vicat_hits.nf": "params.vicat_standardizer_cpus",
            "standardize_vicat.nf": "params.vicat_standardizer_cpus",
            "run_vitap.nf": "params.vitap_cpus",
            "run_tesorter.nf": "params.tesorter_cpus",
        }
        for module_name, declaration in expected.items():
            with self.subTest(module=module_name):
                module = self.read_module(module_name)
                self.assertIn(declaration, module)
                self.assertIn("params.max_cpus", module)
                self.assertIn("Math.min", module)

    def test_supported_tool_thread_flags_receive_task_cpus(self) -> None:
        expected = {
            "run_genomad.nf": "--threads ${task.cpus}",
            "run_virsorter2.nf": "-j ${task.cpus}",
            "run_cenotetaker3.nf": "-t ${task.cpus}",
            "run_virbot.nf": '--threads "${task.cpus}"',
            "run_gianthunter.nf": '--threads "${task.cpus}"',
            "predict_vicat_orfs.nf": '-j "${task.cpus}"',
            "run_vicat_diamond.nf": '--threads "${task.cpus}"',
            "prepare_vicat_reference_subset.nf": '--threads "${task.cpus}"',
            "prepare_vicat_hits.nf": '--threads "${task.cpus}"',
            "standardize_vicat.nf": '--threads "${task.cpus}"',
            "run_vitap.nf": '-p "${task.cpus}"',
            "run_tesorter.nf": '-p "${task.cpus}"',
            "vicat_database.nf": '--threads "${task.cpus}"',
            "virsorter2_database.nf": "-j ${task.cpus}",
        }
        for module_name, argument in expected.items():
            with self.subTest(module=module_name):
                self.assertIn(argument, self.read_module(module_name))

    def test_deep_learning_tools_limit_native_thread_pools(self) -> None:
        for module_name in ("run_deep6.nf", "run_deepmicroclass2.nf"):
            module = self.read_module(module_name)
            with self.subTest(module=module_name):
                self.assertIn('OMP_NUM_THREADS="${task.cpus}"', module)
                self.assertIn('MKL_NUM_THREADS="${task.cpus}"', module)
                self.assertIn('OPENBLAS_NUM_THREADS="${task.cpus}"', module)
                self.assertIn('NUMEXPR_NUM_THREADS="${task.cpus}"', module)

    def test_ambiguous_legacy_thread_parameters_fail_fast(self) -> None:
        self.assertIn("if( params.threads != null )", self.workflow)
        self.assertIn("--threads is no longer used", self.workflow)
        for old_name, new_name in (
            ("vicat_build_threads", "vicat_build_cpus"),
            ("vicat_orf_threads", "vicat_orf_cpus"),
            ("vicat_threads", "vicat_cpus"),
        ):
            self.assertIn(f"params.{old_name} != null", self.workflow)
            self.assertIn(f"--{new_name}", self.workflow)

    def test_analysis_cpu_requests_are_positive_and_capped_by_total_budget(self) -> None:
        self.assertIn("parsePositiveIntegerParameter(params.max_cpus", self.workflow)
        for parameter in (
            "genomad_cpus",
            "virsorter2_cpus",
            "ct3_cpus",
            "deep6_cpus",
            "deepmicroclass2_cpus",
            "virbot_cpus",
            "gianthunter_cpus",
            "vicat_orf_cpus",
            "vicat_cpus",
            "vicat_standardizer_cpus",
            "vitap_cpus",
            "tesorter_cpus",
        ):
            self.assertIn(f"'--{parameter}': params.{parameter}", self.workflow)


if __name__ == "__main__":
    unittest.main()
