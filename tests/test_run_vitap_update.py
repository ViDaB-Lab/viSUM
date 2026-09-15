import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from run_vitap_update import configure_diamond_threads, patch_updater_source


class VitapUpdateCompatibilityTests(unittest.TestCase):
    def test_patch_deduplicates_downloads_and_preserves_coordinate_regions(self):
        upstream = '''def upd(args):
    downloaded_ids = {
        f[:-6] for f in os.listdir(output_folder) if f.endswith(".fasta")
    }
    with ThreadPoolExecutor(max_workers=3) as executor:
        progress_bar = tqdm(total=total_rows, desc="Downloading genomes")
        progress_bar.close()
    # Processing integrated viral sequences
    for row in rows:
        pass
    print("[INFO] All files successfully downloaded and processed")
'''
        patched = patch_updater_source(upstream)

        self.assertIn("unique_download_rows", patched)
        self.assertIn("coordinate_regions", patched)
        self.assertIn('segment_id = f"{virus_id}.segment{segment_index}"', patched)
        self.assertIn("expected_fasta_names", patched)
        self.assertIn("stale_fasta_paths", patched)
        self.assertIn("os.remove(input_fasta)", patched)
        self.assertNotIn("for row in rows:\n        pass", patched)

    def test_patch_fails_loudly_if_upstream_layout_changes(self):
        with self.assertRaisesRegex(RuntimeError, "start marker not found"):
            patch_updater_source("def upd(args): pass")

    def test_diamond_commands_receive_requested_threads(self):
        calls = []

        class FakeSubprocess:
            @staticmethod
            def run(command, *args, **kwargs):
                calls.append(command)
                return 0

        class FakeModule:
            subprocess = FakeSubprocess

        configure_diamond_threads(FakeModule, 20)
        FakeModule.subprocess.run(
            ["diamond", "blastp", "--query", "query.faa"]
        )
        FakeModule.subprocess.run(
            ["diamond", "makedb", "--in", "reference.faa"]
        )
        FakeModule.subprocess.run(
            ["diamond", "blastp", "--threads", "6", "--query", "query.faa"]
        )
        FakeModule.subprocess.run(["wget", "example.invalid"])

        self.assertEqual(
            calls[0][:4], ["diamond", "blastp", "--threads", "20"]
        )
        self.assertEqual(
            calls[1][:4], ["diamond", "makedb", "--threads", "20"]
        )
        self.assertEqual(calls[2][2:4], ["--threads", "6"])
        self.assertEqual(calls[3], ["wget", "example.invalid"])


if __name__ == "__main__":
    unittest.main()
