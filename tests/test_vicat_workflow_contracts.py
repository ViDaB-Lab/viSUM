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
