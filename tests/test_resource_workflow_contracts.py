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

    def read_module(self, name: str) -> str:
        return (MODULES / name).read_text(encoding="utf-8")

    def test_local_executor_has_total_cpu_and_memory_budgets(self) -> None:
        self.assertIn("max_cpus      = 8", self.config)
        self.assertIn("max_memory    = '48 GB'", self.config)
        self.assertIn("executor.$local.cpus = params.max_cpus", self.config)
        self.assertIn("executor.$local.memory = params.max_memory", self.config)

    def test_utility_processes_inherit_one_cpu(self) -> None:
        process_block = re.search(
            r"process\s*\{(?P<body>.*?)\n\}", self.config, re.DOTALL
        )
        self.assertIsNotNone(process_block)
        self.assertRegex(process_block.group("body"), r"\bcpus\s*=\s*1\b")
        self.assertNotIn("params.threads", process_block.group("body"))

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
        ):
            self.assertIn(f"'--{parameter}': params.{parameter}", self.workflow)


if __name__ == "__main__":
    unittest.main()
