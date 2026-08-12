from pathlib import Path


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
