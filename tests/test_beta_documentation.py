"""Keep the beta's taxonomy defaults and public entry points consistent."""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class BetaDocumentationTests(unittest.TestCase):
    def test_taxonomy_refinement_is_enabled_by_default(self):
        config = (ROOT / "visum.config").read_text(encoding="utf-8")
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        for tool in ("vitap", "vcontact3"):
            self.assertRegex(config, rf"(?m)^\s*run_{tool}\s*=\s*true\s*$")
            self.assertRegex(workflow, rf"--run_{tool} BOOL[^\n]*\[true\]")
        self.assertRegex(config, r"(?m)^\s*run_vicat\s*=\s*false\s*$")

    def test_readme_has_four_plain_language_modules(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        diagram = readme.split("```mermaid\n", 1)[1].split("```", 1)[0]
        for n in range(1, 5):
            self.assertIn(f"Module {n}", diagram)
        for internal in ("DISCOVERY_GATE", "PROJECT_VICAT", "discovery gate", "project viCAT"):
            self.assertNotIn(internal, diagram)
        self.assertIn("INPUT -->|DNA| DNA", diagram)
        self.assertIn("INPUT -->|RNA| RNA", diagram)

    def test_local_documentation_links_resolve(self):
        for name in ("README.md", "docs/README.md", "docs/usage.md", "docs/benchmarks/retained-negative-context.md"):
            source = ROOT / name
            for link in re.findall(r"\]\(([^)]+)\)", source.read_text(encoding="utf-8")):
                target = link.split("#", 1)[0]
                if not target or "://" in target:
                    continue
                self.assertTrue((source.parent / target).exists(), f"Broken link in {name}: {link}")


if __name__ == "__main__":
    unittest.main()
