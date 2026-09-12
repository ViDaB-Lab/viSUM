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

    def test_readme_cites_each_external_pipeline_program(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        roles = readme.split("## Programs and roles\n", 1)[1].split("\n## ", 1)[0]
        citations = readme.split("## Who to cite\n", 1)[1]
        programs = re.findall(r"^\| ([^|]+?) \| (?:Yes|No) \|", roles, re.MULTILINE)
        self.assertTrue(programs)
        for program in [*programs, "Nextflow"]:
            if program == "viCAT":
                self.assertIn("viCAT and viHARMONY are components of viSUM", citations)
                continue
            self.assertRegex(
                citations,
                rf"(?m)^\| {re.escape(program)} \| .*https://doi\.org/[^)]+",
            )

    def test_readme_negative_exclusions_are_explicit_and_consistent(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        benchmark = readme.split("## Preliminary performance\n", 1)[1].split("\n## ", 1)[0]
        self.assertNotIn("geNomad", benchmark)
        for excluded in (48, 55):
            retained, total = 101 - excluded, 2390 - excluded
            self.assertIn(f"{retained:,}/{total:,} ({100 * retained / total:.2f}%)", benchmark)
        self.assertIn("conditional, evidence-adjusted", benchmark)
        self.assertIn("27/1,600 (1.69%)", benchmark)
        review = (ROOT / "docs/benchmarks/retained-negative-context.md").read_text(encoding="utf-8")
        cohorts = review.split("The excluded source IDs are:\n", 1)[1].split("The 55-input scenario", 1)[0]
        for prefix, digits, expected in (("CELL_NEG_", 6, 5), ("MITO_NEG_", 6, 5),
                                         ("PLASMID_NEG_", 6, 38), ("NEG_TE_LTR_", 4, 7)):
            ids = re.findall(rf"{prefix}\d{{{digits}}}", cohorts)
            self.assertEqual(len(ids), expected)
            self.assertEqual(len(set(ids)), expected)
        self.assertIn("direct_hmm", review)
        self.assertIn("viral_like_mobile_element_with_viral_support", review)


if __name__ == "__main__":
    unittest.main()
