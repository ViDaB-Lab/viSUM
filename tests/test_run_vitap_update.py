import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from run_vitap_update import patch_updater_source


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
        self.assertIn("os.remove(input_fasta)", patched)
        self.assertNotIn("for row in rows:\n        pass", patched)

    def test_patch_fails_loudly_if_upstream_layout_changes(self):
        with self.assertRaisesRegex(RuntimeError, "start marker not found"):
            patch_updater_source("def upd(args): pass")


if __name__ == "__main__":
    unittest.main()
