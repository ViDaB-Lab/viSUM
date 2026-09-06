#!/usr/bin/env python3

"""Extract a DIAMOND database sequence count from ``diamond dbinfo`` output."""

from __future__ import annotations

import re
import sys


SEQUENCE_COUNT = re.compile(
    r"^\s*(?:Database\s+)?sequences\s*:?\s+([0-9][0-9,]*)\s*$",
    re.IGNORECASE,
)


def parse_sequence_count(text: str) -> int:
    counts = {
        int(match.group(1).replace(",", ""))
        for line in text.splitlines()
        if (match := SEQUENCE_COUNT.match(line))
    }
    if not counts:
        raise ValueError("DIAMOND dbinfo output contains no sequence count")
    if len(counts) != 1:
        raise ValueError(
            "DIAMOND dbinfo output contains conflicting sequence counts: "
            + ", ".join(str(count) for count in sorted(counts))
        )
    return counts.pop()


def main() -> None:
    try:
        print(parse_sequence_count(sys.stdin.read()))
    except ValueError as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
