import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class CheckVDatabaseWorkflowContractTests(unittest.TestCase):
    def test_workflow_wires_optional_checkv_database_setup(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("PREPARE_CHECKV_DATABASE", workflow)
        self.assertIn("params.run_checkv", workflow)
        self.assertIn("params.checkv_db", workflow)
        self.assertIn("params.checkv_auto_download", workflow)

    def test_database_module_validates_both_database_components(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/checkv_database.nf").read_text(
            encoding="utf-8"
        )
        for required_path in (
            "genome_db/checkv_reps.faa",
            "genome_db/checkv_reps.fna",
            "genome_db/checkv_reps.tsv",
            "genome_db/checkv_reps.dmnd",
            "hmm_db/checkv_hmms.tsv",
            "hmm_db/genome_lengths.tsv",
            "hmm_db/checkv_hmms",
        ):
            self.assertIn(required_path, module)
        self.assertIn("diamond dbinfo", module)

    def test_database_install_is_staged_locked_and_reusable(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/checkv_database.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("flock 9", module)
        self.assertIn("checkv download_database", module)
        self.assertIn("skipped-existing-managed-database", module)
        self.assertIn("skipped-user-supplied-database", module)

    def test_analysis_runs_only_on_discovery_candidates(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("include { RUN_CHECKV }", workflow)
        self.assertIn("DISCOVERY_GATE.out.candidates", workflow)
        self.assertIn("RUN_CHECKV(", workflow)

    def test_checkv_standardization_and_provirus_refinement_are_wired(self) -> None:
        workflow = (REPOSITORY_ROOT / "visum_nextflow.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("include { STANDARDIZE_CHECKV }", workflow)
        self.assertIn("STANDARDIZE_CHECKV(ch_checkv_standardizer_input)", workflow)
        self.assertIn("STANDARDIZE_CHECKV.out.evidence", workflow)
        self.assertIn("include { REFINE_PROVIRAL_REGIONS }", workflow)
        self.assertIn("REFINE_PROVIRAL_REGIONS(ch_provirus_refinement_inputs)", workflow)
        self.assertIn("STANDARDIZE_GENOMAD.out.evidence", workflow)
        self.assertIn("STANDARDIZE_CENOTETAKER3.out.evidence", workflow)

    def test_refinement_module_publishes_fasta_mapping_and_audit(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/refine_proviral_regions.nf").read_text(
            encoding="utf-8"
        )
        for expected_output in (
            "refined_candidates.fasta",
            "provirus_region_map.tsv",
            "provirus_boundary_audit.tsv",
            "provirus_refinement_summary.tsv",
        ):
            self.assertIn(expected_output, module)

    def test_analysis_preserves_primary_reports_and_zero_candidate_runs(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/run_checkv.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn("checkv end_to_end", module)
        self.assertIn("skipped_no_discovery_candidates", module)
        self.assertIn("quality-summary coverage mismatch", module)
        for expected_output in (
            "checkv_quality_summary.tsv",
            "checkv_completeness.tsv",
            "checkv_contamination.tsv",
            "checkv_complete_genomes.tsv",
            "checkv_proviruses.fna",
            "checkv_run_metadata.tsv",
        ):
            self.assertIn(expected_output, module)

    def test_zero_candidate_headers_match_checkv_1_1_1(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/run_checkv.nf").read_text(
            encoding="utf-8"
        )
        quality_header = (
            "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count"
            "\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality"
            "\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq"
            "\twarnings"
        )
        completeness_header = (
            "contig_id\tcontig_length\tviral_length\taai_expected_length"
            "\taai_completeness\taai_confidence\taai_error\taai_num_hits"
            "\taai_top_hit\taai_id\taai_af\thmm_completeness_lower"
            "\thmm_completeness_upper\thmm_num_hits\tkmer_freq"
        )
        complete_genomes_header = (
            "contig_id\tcontig_length\tkmer_freq\tprediction_type"
            "\tconfidence_level\tconfidence_reason\trepeat_length"
            "\trepeat_count\trepeat_n_freq\trepeat_mode_base_freq\trepeat_seq"
        )
        self.assertIn(quality_header, module)
        self.assertIn(completeness_header, module)
        self.assertIn(complete_genomes_header, module)
        self.assertNotIn("completeness_method\tcomplete_genome_type", module)

    def test_run_metadata_records_resolved_database_path(self) -> None:
        module = (REPOSITORY_ROOT / "modules/local/run_checkv.nf").read_text(
            encoding="utf-8"
        )
        self.assertIn('readlink -f "${checkv_database}"', module)
        self.assertIn('"\\$RESOLVED_DATABASE_PATH" "\\$DATABASE_RELEASE"', module)


if __name__ == "__main__":
    unittest.main()
