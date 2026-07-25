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
