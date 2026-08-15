from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from run_vitap_assignment import patch_assignment_source


class VitapAssignmentRunnerTests(unittest.TestCase):
    def test_replaces_the_upstream_uniref90_thread_literal(self) -> None:
        source = (
            'first = ["-p", threads]\n'
            'second = ["--threads", "10", "--quiet"]\n'
        )

        patched = patch_assignment_source(source)

        self.assertIn('["--threads", threads, "--quiet"]', patched)
        self.assertNotIn('["--threads", "10", "--quiet"]', patched)
        self.assertIn('["-p", threads]', patched)

    def test_rejects_an_unrecognized_upstream_source(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "upstream implementation may have changed"):
            patch_assignment_source('subprocess.run(["diamond", "blastp"])')


if __name__ == "__main__":
    unittest.main()
