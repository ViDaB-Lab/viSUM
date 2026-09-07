import gzip
from pathlib import Path

from bin.prepare_vicat_competitive_references import normalize_cellular_metadata


ROOT = Path(__file__).resolve().parents[1]


def test_vicat_publish_patterns_do_not_reference_process_inputs() -> None:
    for name in (
        "predict_vicat_orfs.nf",
        "run_vicat_diamond.nf",
        "run_vicat_nonviral_diamond.nf",
        "standardize_vicat.nf",
    ):
        text = (ROOT / "modules" / "local" / name).read_text(encoding="utf-8")
        pattern_lines = [line.strip() for line in text.splitlines() if "pattern:" in line]
        assert len(pattern_lines) == 1
        assert "*.vicat_" in pattern_lines[0]
        if name.startswith("run_vicat_"):
            assert "saveAs:" in text
            assert "filename.contains('.raw.') ? null" in text


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


def test_vicat_standardizer_declares_and_enforces_duckdb_threads() -> None:
    module = (ROOT / "modules" / "local" / "standardize_vicat.nf").read_text(
        encoding="utf-8"
    )
    helper = (ROOT / "bin" / "standardize_vicat.py").read_text(encoding="utf-8")
    config = (ROOT / "visum.config").read_text(encoding="utf-8")
    assert "params.vicat_standardizer_cpus" in module
    assert "params.max_cpus" in module
    assert '--threads "${task.cpus}"' in module
    assert "vicat_standardizer_cpus = 8" in config
    assert 'parser.add_argument("--threads", type=int, default=1)' in helper
    assert 'connection.execute(f"SET threads = {threads}")' in helper
    assert 'args.reference_manifest, args.threads' in helper


def test_vicat_hit_preparation_is_cached_before_margin_classification() -> None:
    workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
    subset_module = (
        ROOT / "modules" / "local" / "prepare_vicat_reference_subset.nf"
    ).read_text(encoding="utf-8")
    preparation = (ROOT / "modules" / "local" / "prepare_vicat_hits.nf").read_text(
        encoding="utf-8"
    )
    standardizer = (ROOT / "modules" / "local" / "standardize_vicat.nf").read_text(
        encoding="utf-8"
    )
    subset_helper = (ROOT / "bin" / "prepare_vicat_reference_subset.py").read_text(
        encoding="utf-8"
    )
    helper = (ROOT / "bin" / "prepare_vicat_hits.py").read_text(encoding="utf-8")
    assert "include { PREPARE_VICAT_REFERENCE_SUBSET }" in workflow
    assert "include { PREPARE_VICAT_HITS }" in workflow
    assert ".collect()" in workflow
    assert "ch_vicat_reference_subset = PREPARE_VICAT_REFERENCE_SUBSET.out.subset" in workflow
    assert "PREPARE_VICAT_REFERENCE_SUBSET.out.subset.first()" not in workflow
    assert "PREPARE_VICAT_HITS(" in workflow
    assert "ch_vicat_standardizer_script = Channel.value(" in workflow
    assert "PREPARE_VICAT_HITS.out.results,\n            ch_vicat_standardizer_script" in workflow
    assert "path standardizer_script" in standardizer
    assert 'python "${standardizer_script}"' in standardizer
    assert "vicat_competitive_min_margin" not in subset_module
    assert "vicat_competitive_min_margin" not in preparation
    assert "--prepared-hits" in standardizer
    assert "ORDER BY" not in subset_helper
    assert "ORDER BY" not in helper
    assert "relevant_reference_ids" in subset_helper
    assert "FORMAT PARQUET, COMPRESSION ZSTD" in helper
    assert "--viral-reference-subset" in preparation
    assert "--viral-reference-subset-metadata" in preparation
    assert "--nonviral-metadata" in preparation
    assert "SET memory_limit" in subset_helper
    assert "SET temp_directory" in subset_helper


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


def test_vicat_guarded_dna_single_locus_rescue_is_wired_through_workflow() -> None:
    module = (ROOT / "modules" / "local" / "standardize_vicat.nf").read_text(
        encoding="utf-8"
    )
    config = (ROOT / "visum.config").read_text(encoding="utf-8")
    workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")

    assert 'vicat_dna_single_locus_rescue = \'strict\'' in config
    assert '--dna-single-locus-rescue "${params.vicat_dna_single_locus_rescue}"' in module
    assert "--vicat_dna_single_locus_rescue off|strict" in workflow
    assert "vicatDnaSingleLocusRescue in ['off', 'strict']" in workflow


def test_legacy_combined_database_builder_is_not_exposed_by_runtime_config() -> None:
    config = (ROOT / "visum.config").read_text(encoding="utf-8")
    workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
    module = (ROOT / "modules" / "local" / "vicat_database.nf").read_text(encoding="utf-8")
    builder = (ROOT / "bin" / "build_vicat_competitive_database.sh").read_text(encoding="utf-8")
    tier_builder = (ROOT / "bin" / "build_vicat_competitive_tiers.sh").read_text(
        encoding="utf-8"
    )
    helper = (ROOT / "bin" / "prepare_vicat_competitive_references.py").read_text(encoding="utf-8")

    for parameter in (
        "vicat_viral_representatives", "vicat_cellular_proteins",
        "vicat_cellular_metadata", "vicat_competitive_managed_db",
    ):
        assert parameter not in config
        assert parameter not in workflow

    assert "vicat_nonviral_db" in config
    assert "vicat_nonviral_db" in workflow
    assert "build_vicat_competitive_database.sh" in module
    assert "mmseqs linclust" in builder
    assert "vicat_viral_cellular.dmnd" in builder
    assert 'bash "$SCRIPT_DIR/build_vicat_competitive_database.sh"' in tier_builder
    assert 'f">{label}|{source_id}\\n"' in helper


def test_vicat_nonviral_builder_clusters_each_reference_class_independently() -> None:
    builder = (ROOT / "bin" / "build_vicat_nonviral_database.sh").read_text(
        encoding="utf-8"
    )
    helper = (
        ROOT / "bin" / "prepare_vicat_nonviral_cluster_metadata.py"
    ).read_text(encoding="utf-8")

    for reference_class in (
        "cellular_chromosome",
        "cellular_unplaced",
        "plasmid",
        "plastid",
        "mitochondrial",
        "shared_nonviral",
    ):
        assert reference_class in builder

    assert 'for reference_class in "${CLASSES[@]}"' in builder
    assert "mmseqs linclust" in builder
    assert '--min-seq-id "$MIN_SEQ_ID" --cov-mode 0 -c "$COVERAGE"' in builder
    assert "--min-len" not in builder
    assert "vicat_nonviral_cluster_membership.tsv.gz" in builder
    assert "vicat_nonviral_representative_metadata.parquet" in builder
    assert "vicat_nonviral.dmnd" in builder
    assert "provirus_flank_eligible" in helper
    assert "cluster_member_count" in helper


def test_vicat_nonviral_database_module_has_no_nested_groovy_triple_quotes() -> None:
    module = (
        ROOT / "modules" / "local" / "vicat_nonviral_database.nf"
    ).read_text(encoding="utf-8")
    # The process script itself uses one Groovy triple-quoted string. Python
    # embedded in its heredoc must not contain another triple-quoted string,
    # which would terminate the Groovy script at Nextflow compile time.
    assert module.count('"""') == 2
    assert "SELECT count(*) FROM read_parquet(?)" in module
    assert "parse_diamond_dbinfo.py" in module


def test_vicat_database_validators_share_version_tolerant_dbinfo_parser() -> None:
    viral = (ROOT / "modules" / "local" / "vicat_database.nf").read_text(
        encoding="utf-8"
    )
    nonviral = (
        ROOT / "modules" / "local" / "vicat_nonviral_database.nf"
    ).read_text(encoding="utf-8")

    assert viral.count("parse_diamond_dbinfo.py") == 2
    assert "parse_diamond_dbinfo.py" in nonviral
    assert "awk '$1 == \"Sequences\"" not in viral
    assert 'startswith("Sequences")' not in nonviral


def test_vicat_analysis_searches_viral_and_class_aware_nonviral_databases() -> None:
    analysis = (ROOT / "modules" / "local" / "run_vicat_diamond.nf").read_text(
        encoding="utf-8"
    )
    preparation = (
        ROOT / "modules" / "local" / "prepare_vicat_reference_subset.nf"
    ).read_text(
        encoding="utf-8"
    )
    helper = (ROOT / "bin" / "prepare_vicat_reference_subset.py").read_text(
        encoding="utf-8"
    )
    standardizer_helper = (ROOT / "bin" / "standardize_vicat.py").read_text(
        encoding="utf-8"
    )
    workflow = (ROOT / "visum_nextflow.nf").read_text(encoding="utf-8")
    nonviral_analysis = (
        ROOT / "modules" / "local" / "run_vicat_nonviral_diamond.nf"
    ).read_text(encoding="utf-8")

    assert "IMGVR5_UViG_representatives.dmnd" in analysis
    assert "vicat_viral_cellular.dmnd" not in analysis
    assert "vicat_nonviral.dmnd" in nonviral_analysis
    assert "RUN_VICAT_NONVIRAL_DIAMOND" in workflow
    assert "ch_vicat_dual_alignments" in workflow
    assert "--reference-manifest" not in preparation
    assert "reference_class" in helper
    assert 'hit["reference_class"] == "viral"' in standardizer_helper
    assert "best_nonviral_reference_class" in standardizer_helper


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
