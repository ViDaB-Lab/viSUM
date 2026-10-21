from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class ViharmonyWorkflowContractTests(unittest.TestCase):
    def test_workflow_runs_harmonizer_after_refinement_and_taxonomy_tools(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        self.assertIn("include { VIHARMONY }", workflow)
        self.assertGreater(workflow.index("VIHARMONY(ch_harmony_inputs)"), workflow.index("STANDARDIZE_VCONTACT3(RUN_VCONTACT3.out.results)"))
        self.assertIn("ch_harmony_evidence_for_sample", workflow)
        self.assertIn("assets/empty_discovery_evidence.tsv", workflow)
        self.assertIn("ch_refined_for_harmony", workflow)
        self.assertIn("PROJECT_VICAT_REFINED.out.projection", workflow)

    def test_harmonizer_requires_each_enabled_per_sample_evidence_artifact(self) -> None:
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")

        self.assertIn("expectedHarmonyToolsByType", workflow)
        self.assertIn("expectedHarmonyToolsByType[type].size() +", workflow)
        self.assertIn("extraHarmonyEvidenceCount", workflow)
        self.assertIn("'viharmony', expectedTools, actualTools, evidenceFiles", workflow)
        self.assertIn("ch_harmony_inputs = ch_refined_for_harmony", workflow)

    def test_module_declares_agreed_final_outputs(self) -> None:
        module = (ROOT / "modules" / "local" / "viharmony.nf").read_text(encoding="utf-8")
        for suffix in (
            ".final.normalized.fasta", ".final.original_ids.fasta",
            ".final_metadata.tsv", ".review_queue.tsv",
            ".sequence_disposition.tsv", ".sequence_map.tsv",
            ".harmonizer_manifest.json", ".database_candidates.fasta",
            ".database_candidates.tsv",
            ".all_candidates.fasta",
            ".provisional.normalized.fasta",
            ".provisional.original_ids.fasta", ".provisional_metadata.tsv",
        ):
            self.assertIn(suffix, module)
        self.assertIn("--audit-mode", module)
        self.assertIn("--vcontact3-min-taxonomy-length", module)
        self.assertNotIn('.review_candidates.fasta', module)
        self.assertIn('path("audit/${prefix}.*_audit.tsv.gz"),\n          optional: true', module)
        self.assertIn("viHARMONY decision policy v1.4", module)
        self.assertIn("TE-routing-v1", module)

    def test_config_uses_one_canonical_msl_parameter(self) -> None:
        config = (ROOT / "visum.config").read_text(encoding="utf-8")
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        self.assertIn("harmonizer_audit", config)
        self.assertIn('ictv_csv = "${projectDir}/assets/ICTV_VMR_MSL41.csv"', config)
        self.assertNotIn("harmonizer_ictv_csv", config)
        self.assertIn("def ictvCsv = file(params.ictv_csv)", workflow)
        self.assertNotIn("harmonizerIctvCsv", workflow)
        self.assertIn("vcontact3_min_taxonomy_length", config)

    def test_workers_respect_nextflow_cpu_allocation(self) -> None:
        config = (ROOT / "visum.config").read_text(encoding="utf-8")
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        module = (ROOT / "modules/local/viharmony.nf").read_text(encoding="utf-8")
        self.assertIn("harmonizer_cpus   = 1", config)
        self.assertIn("Math.min(params.harmonizer_cpus as int, params.max_cpus as int)", module)
        self.assertIn('--workers "${task.cpus}"', module)
        self.assertIn("'--harmonizer_cpus': params.harmonizer_cpus", workflow)
        self.assertIn("--harmonizer_cpus INT", workflow)

    def test_reporting_cap_defaults_and_native_command_are_wired(self) -> None:
        config = (ROOT / "visum.config").read_text(encoding="utf-8")
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        module = (ROOT / "modules/local/viharmony.nf").read_text(encoding="utf-8")
        self.assertIn("harmonizer_max_rank = 'family'", config)
        self.assertIn("--harmonizer_max_rank RANK", workflow)
        self.assertIn("Invalid --harmonizer_max_rank", workflow)
        self.assertIn('def maxRank = params.harmonizer_max_rank.toString().trim().toLowerCase()', module)
        self.assertIn('--max-rank "${maxRank}"', module)
        self.assertIn('audit/${prefix}.*_audit.tsv.gz', module)

    def test_published_results_refresh_on_reporting_changes_and_cached_restore(self) -> None:
        module = (ROOT / "modules/local/viharmony.nf").read_text(encoding="utf-8")
        self.assertIn('"${params.outdir}/${prefix}_results/viharmony"', module)
        self.assertIn("mode: 'copy', overwrite: true, failOnError: true", module)

    def test_reconciliation_bundle_and_manifest_are_native_task_inputs(self) -> None:
        config = (ROOT / "visum.config").read_text(encoding="utf-8")
        workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
        module = (ROOT / "modules/local/viharmony.nf").read_text(encoding="utf-8")
        self.assertIn('ictv_reconciliation = "${projectDir}/assets/ictv_reconciliation/', config)
        self.assertIn("path(ictv_reconciliation)", module)
        self.assertIn("path(ictv_reconciliation_manifest)", module)
        self.assertIn('--ictv-reconciliation "${ictv_reconciliation}"', module)
        self.assertIn("groupFiles, ictvCsv, ictvReconciliation, ictvReconciliationManifest, harmonizerAudit", workflow)
        self.assertIn("${prefix}.*_audit.tsv.gz", module)


if __name__ == "__main__":
    unittest.main()
