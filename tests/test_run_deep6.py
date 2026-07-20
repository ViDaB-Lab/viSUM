import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "bin"))

import run_deep6


class Deep6RunnerTests(unittest.TestCase):
    def test_500_model_does_not_match_1500_model(self):
        with tempfile.TemporaryDirectory() as temp_directory:
            model_directory = Path(temp_directory)
            model_500 = model_directory / "model_500_test.h5"
            model_1500 = model_directory / "model_1500_test.h5"
            model_500.write_bytes(b"500")
            model_1500.write_bytes(b"1500")

            selected = run_deep6.select_model(model_directory, "500")

            self.assertEqual(selected, model_500)

    def test_duplicate_exact_models_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp_directory:
            model_directory = Path(temp_directory)
            (model_directory / "model_250_first.h5").write_bytes(b"first")
            (model_directory / "model_250_second.h5").write_bytes(b"second")

            with self.assertRaisesRegex(ValueError, "found 2"):
                run_deep6.select_model(model_directory, "250")

    def test_missing_model_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_directory:
            model_directory = Path(temp_directory)

            with self.assertRaisesRegex(ValueError, "found 0"):
                run_deep6.select_model(model_directory, "1000")


if __name__ == "__main__":
    unittest.main()
