import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class VitapDatabaseWorkflowContractTests(unittest.TestCase):
    def test_setup_is_explicitly_gated_and_wired(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        self.assertIn("include { PREPARE_VITAP_DATABASE }", workflow)
        self.assertIn("params.run_vitap", workflow)
        self.assertIn("params.vitap_db", workflow)
        self.assertIn("params.vitap_vmr", workflow)
        self.assertIn("params.vitap_update_database", workflow)
        self.assertIn("PREPARE_VITAP_DATABASE(ch_vitap_database_request)", workflow)

    def test_module_validates_runtime_files_and_diamond_databases(self) -> None:
        module = (ROOT / "modules/local/vitap_database.nf").read_text(
            encoding="utf-8"
        )
        for required in (
            "VMR_genome_*.gff",
            "uniref90.dmnd",
            "uniref90.accession2taxid",
            "taxdmp/nodes.dmp",
            "taxdmp/names.dmp",
            "Realm Kingdom Phylum Class Order Family Genus Species",
            "diamond dbinfo",
        ):
            self.assertIn(required, module)

    def test_managed_build_is_locked_reusable_and_release_aware(self) -> None:
        module = (ROOT / "modules/local/vitap_database.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("flock 9", module)
        self.assertIn("https://ictv.global/vmr/current", (ROOT / "visum.config").read_text())
        self.assertIn("skipped-existing-managed-database", module)
        self.assertIn("skipped-current-release-already-installed", module)
        self.assertIn("skipped-user-supplied-database", module)
        self.assertIn(".visum_vmr_metadata.tsv", module)
        self.assertIn(".visum_db_complete", module)

    def test_environment_pins_current_vitap_release(self) -> None:
        environment = (ROOT / "envs/vitap.yml").read_text(encoding="utf-8")
        self.assertIn("vitap=1.12", environment)
        self.assertIn("openpyxl", environment)


if __name__ == "__main__":
    unittest.main()
