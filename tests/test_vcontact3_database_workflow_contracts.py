import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "bin" / "validate_vcontact3_database.py"


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
            'file("${projectDir}/bin/validate_vcontact3_database.py")',
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
            'python "$VALIDATOR"',
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

    def test_validator_accepts_complete_manifest_and_rejects_missing_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            manifest = {
                "232": {
                    "RefSeq": {
                        "domains": {
                            domain: {
                                "Genomes Report": f"v232/{domain}/genomes.parquet",
                                "gene2genome": f"v232/{domain}/gene2genome.parquet",
                                "proteins": f"v232/{domain}/proteins.faa",
                                "nucleotides": f"v232/{domain}/nucleotides.fna",
                                "identities": {
                                    "95": f"v232/{domain}/identity95.parquet"
                                },
                            }
                            for domain in ("prokaryotes", "eukaryotes")
                        }
                    },
                    "VOGDB": {
                        "mmseqs_db": "v232/VOGDB/vog_profiles",
                        "rank_exclusivity": "v232/VOGDB/rank_exclusivity.parquet",
                        "host_domain": "v232/VOGDB/host_domain.parquet",
                    },
                }
            }
            (root / "232.json").write_text(json.dumps(manifest), encoding="utf-8")

            referenced_files = []
            for domain in ("prokaryotes", "eukaryotes"):
                domain_data = manifest["232"]["RefSeq"]["domains"][domain]
                referenced_files.extend(
                    domain_data[key]
                    for key in ("Genomes Report", "gene2genome", "proteins", "nucleotides")
                )
                referenced_files.append(domain_data["identities"]["95"])
            referenced_files.extend(
                manifest["232"]["VOGDB"][key]
                for key in ("rank_exclusivity", "host_domain")
            )
            referenced_files.append("v232/VOGDB/vog_profiles.dbtype")

            for relative_path in referenced_files:
                target = root / relative_path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text("data\n", encoding="utf-8")

            accepted = subprocess.run(
                [sys.executable, str(VALIDATOR), str(root), "--field", "version"],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(accepted.returncode, 0, accepted.stderr)
            self.assertEqual(accepted.stdout.strip(), "232")

            (root / "v232/eukaryotes/proteins.faa").unlink()
            rejected = subprocess.run(
                [sys.executable, str(VALIDATOR), str(root)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("missing or empty", rejected.stderr)


if __name__ == "__main__":
    unittest.main()
