import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
class Vcontact3DatabaseWorkflowContractTests(unittest.TestCase):
    def test_setup_is_gated_and_wired(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        self.assertIn("params.run_vcontact3", workflow)
        self.assertNotIn("include { PREPARE_VCONTACT3_DATABASE }", workflow)
        setup = (ROOT / 'workflows/database_setup.nf').read_text()
        for expected in (
            "include { PREPARE_VCONTACT3_DATABASE }",
            "params.vcontact3_db",
            "params.vcontact3_auto_download",
            "params.vcontact3_update_database",
            "params.vcontact3_cleanup_archive",
            "PREPARE_VCONTACT3_DATABASE(ch_vcontact3_database_request)",
        ):
            self.assertIn(expected, setup)

    def test_module_uses_official_latest_command_and_managed_reuse(self) -> None:
        module = (ROOT / "modules/local/vcontact3_database.nf").read_text(
            encoding="utf-8"
        )
        for expected in (
            "vcontact3 prepare_databases",
            "--get-version latest",
            "python -m json.tool",
            "vcontact3 prepare_databases",
            "flock 9",
            "skipped-user-supplied-database",
            "skipped-existing-managed-database",
            "downloaded-latest-and-installed",
            ".visum_db_complete",
            'releases/$VERSION',
            'ln -sfn "releases/$VERSION" "$DB_ROOT/current"',
        ):
            self.assertIn(expected, module)

    def test_environment_pins_compatible_official_release(self) -> None:
        environment = (ROOT / "envs/vcontact3.yml").read_text(encoding="utf-8")
        self.assertIn("d57ae81e0be5ab57a8a14692a4cb8ab255ce2d52", environment)
        self.assertIn("mmseqs2", environment)
        self.assertIn("pandas=2.1.4", environment)
        self.assertIn("pyarrow=15.0.2", environment)

    def test_setup_checks_parquet_before_any_database_preparation(self) -> None:
        module = (ROOT / 'modules/local/vcontact3_database.nf').read_text()
        runtime = module.split('    require_compatible_runtime() {', 1)[1].split('\n    }', 1)[0]
        self.assertIn('python -B "!{projectDir}/bin/check_vcontact3_parquet.py"', runtime)
        self.assertLess(module.index('    require_compatible_runtime\n'),
                        module.index('    vcontact3 prepare_databases'))

if __name__ == "__main__":
    unittest.main()
