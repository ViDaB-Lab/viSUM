import gzip
from pathlib import Path

from bin.prepare_vicat_competitive_references import normalize_cellular_metadata


ROOT = Path(__file__).resolve().parents[1]


def test_vicat_publish_patterns_do_not_reference_process_inputs() -> None:
    for name in (
        "predict_vicat_orfs.nf",
        "run_vicat_diamond.nf",
        "standardize_vicat.nf",
    ):
        text = (ROOT / "modules" / "local" / name).read_text(encoding="utf-8")
        pattern_lines = [line.strip() for line in text.splitlines() if "pattern:" in line]
        assert pattern_lines == ["pattern: '*.vicat_*'"]


def test_vicat_orf_commands_use_supported_cli_flags_and_pinned_versions() -> None:
    module = (ROOT / "modules" / "local" / "predict_vicat_orfs.nf").read_text(encoding="utf-8")
    environment = (ROOT / "envs" / "vicat.yml").read_text(encoding="utf-8")
    assert " -q " not in module
    assert "pyrodigal-gv=0.3.2" in environment
    assert "pyrodigal-rv=0.1.0" in environment


def test_vicat_diamond_uses_one_literal_command_and_lean_tuple() -> None:
    module = (ROOT / "modules" / "local" / "run_vicat_diamond.nf").read_text(encoding="utf-8")
    commands = [line.strip() for line in module.splitlines() if line.strip().startswith("diamond blastp")]
    assert len(commands) == 1
    assert "bin/run_vicat_diamond.py" not in module
    assert "DIAMOND_ARGS=(" not in module
    assert "path(orf_metadata)" not in module
    assert "path(orf_log)" not in module
    for parameter in (
        "params.vicat_diamond_sensitivity",
        "params.vicat_min_bitscore", "params.vicat_min_query_cover",
        "params.vicat_top_percent", "params.vicat_block_size",
        "params.vicat_index_chunks",
    ):
        assert parameter in commands[0]
    assert "--evalue" not in commands[0]
    assert " evalue " in commands[0]


def test_vicat_standardizer_emits_reference_audit() -> None:
    module = (ROOT / "modules" / "local" / "standardize_vicat.nf").read_text(encoding="utf-8")
    assert 'path("${prefix}.vicat_reference_audit.tsv"), emit: audit' in module
    assert '--output-audit "${prefix}.vicat_reference_audit.tsv"' in module


def test_vicat_uses_separate_orf_and_contig_taxonomy_thresholds() -> None:
    module = (ROOT / "modules" / "local" / "standardize_vicat.nf").read_text(encoding="utf-8")
    config = (ROOT / "visum.config").read_text(encoding="utf-8")
    workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
    for parameter in ("vicat_orf_taxonomy_support", "vicat_contig_taxonomy_support"):
        assert parameter in module
        assert parameter in config
        assert parameter in workflow
    assert "vicat_taxonomy_support" not in module
    assert "vicat_taxonomy_support" not in config
    assert "vicat_taxonomy_support" not in workflow


def test_vicat_competitive_database_extension_is_explicit_and_labeled() -> None:
    config = (ROOT / "visum.config").read_text(encoding="utf-8")
    workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
    module = (ROOT / "modules" / "local" / "vicat_database.nf").read_text(encoding="utf-8")
    builder = (ROOT / "bin" / "build_vicat_competitive_database.sh").read_text(encoding="utf-8")
    tier_builder = (ROOT / "bin" / "build_vicat_competitive_tiers.sh").read_text(
        encoding="utf-8"
    )
    helper = (ROOT / "bin" / "prepare_vicat_competitive_references.py").read_text(encoding="utf-8")

    for parameter in (
        "vicat_viral_representatives",
        "vicat_cellular_proteins",
        "vicat_cellular_metadata",
        "vicat_competitive_managed_db",
        "vicat_cellular_min_seq_id",
        "vicat_cellular_coverage",
        "vicat_cellular_min_protein_length",
    ):
        assert parameter in config
        assert parameter in workflow

    assert "viSUM-managed-competitive" in workflow
    assert "build_vicat_competitive_database.sh" in module
    assert "mmseqs linclust" in builder
    assert "vicat_viral_cellular.dmnd" in builder
    assert 'bash "$SCRIPT_DIR/build_vicat_competitive_database.sh"' in tier_builder
    assert 'f">{label}|{source_id}\\n"' in helper


def test_vicat_analysis_prefers_competitive_database_and_wires_manifest() -> None:
    analysis = (ROOT / "modules" / "local" / "run_vicat_diamond.nf").read_text(
        encoding="utf-8"
    )
    standardizer = (ROOT / "modules" / "local" / "standardize_vicat.nf").read_text(
        encoding="utf-8"
    )
    helper = (ROOT / "bin" / "standardize_vicat.py").read_text(encoding="utf-8")

    assert "vicat_viral_cellular.dmnd" in analysis
    assert "vicat_competitive_reference_manifest.parquet" in analysis
    assert "--reference-manifest" in standardizer
    assert "reference_class" in helper
    assert 'hit["reference_class"] == "viral"' in helper


def test_vicat_projects_precomputed_loci_after_provirus_refinement() -> None:
    workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
    module = (ROOT / "modules" / "local" / "project_vicat_refined.nf").read_text(
        encoding="utf-8"
    )
    helper = (ROOT / "bin" / "project_vicat_refined.py").read_text(encoding="utf-8")
    assert "PROJECT_VICAT_REFINED(ch_vicat_refinement_projection)" in workflow
    assert "REFINE_PROVIRAL_REGIONS.out.refined" in workflow
    assert "STANDARDIZE_VICAT.out.loci" in workflow
    assert "project_vicat_refined.py" in module
    assert "crosses_boundary" in helper
    assert '"evidence_scope": "refined_region"' in helper


def test_vicat_competitive_metadata_is_normalized_before_duckdb() -> None:
    helper = (ROOT / "bin" / "prepare_vicat_competitive_references.py").read_text(
        encoding="utf-8"
    )
    assert "def normalize_cellular_metadata" in helper
    assert 'elif "\\\\t" in first_line:' in helper
    assert "Observed columns:" in helper
    assert "cellular_metadata.normalized.tsv" in helper


def test_vicat_competitive_metadata_detects_gzip_without_suffix(tmp_path: Path) -> None:
    source = tmp_path / "staged_cellular_metadata"
    destination = tmp_path / "normalized.tsv"
    with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
        handle.write(
            "protein_id\tcellular_group\tsource_accession\torganism_name\ttaxid\n"
            "protein_1\tbacteria\tGCF_000001\tExample bacterium\t1234\n"
        )

    columns = normalize_cellular_metadata(source, destination)

    assert columns == {
        "protein_id",
        "cellular_group",
        "source_accession",
        "organism_name",
        "taxid",
    }
    assert destination.read_text(encoding="utf-8") == (
        "protein_id\tcellular_group\tsource_accession\torganism_name\ttaxid\n"
        "protein_1\tbacteria\tGCF_000001\tExample bacterium\t1234\n"
    )
