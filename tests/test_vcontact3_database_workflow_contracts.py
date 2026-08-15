import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
class Vcontact3DatabaseWorkflowContractTests(unittest.TestCase):
    def test_setup_is_gated_and_wired(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        for expected in (
            "include { PREPARE_VCONTACT3_DATABASE }",
            "params.run_vcontact3",
            "params.vcontact3_db",
            "params.vcontact3_auto_download",
            "params.vcontact3_update_database",
            "params.vcontact3_cleanup_archive",
            "PREPARE_VCONTACT3_DATABASE(ch_vcontact3_database_request)",
        ):
            self.assertIn(expected, workflow)

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

if __name__ == "__main__":
    unittest.main()
