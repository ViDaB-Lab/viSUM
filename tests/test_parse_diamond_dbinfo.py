import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from parse_diamond_dbinfo import parse_sequence_count


@pytest.mark.parametrize(
    ("output", "expected"),
    [
        ("Sequences 178734704\n", 178734704),
        ("Database sequences 178734704\nDatabase letters 43900319257\n", 178734704),
        ("Database sequences: 178,734,704\n", 178734704),
    ],
)
def test_parses_supported_diamond_dbinfo_formats(output: str, expected: int) -> None:
    assert parse_sequence_count(output) == expected


def test_rejects_missing_sequence_count() -> None:
    with pytest.raises(ValueError, match="no sequence count"):
        parse_sequence_count("Database letters 43900319257\n")


def test_command_line_parser_reads_standard_input() -> None:
    completed = subprocess.run(
        [sys.executable, str(ROOT / "bin" / "parse_diamond_dbinfo.py")],
        input="Database sequences 42\n",
        text=True,
        capture_output=True,
        check=True,
    )
    assert completed.stdout.strip() == "42"
