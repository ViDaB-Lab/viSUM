#!/usr/bin/env python3

"""Resolve standardized provirus boundaries and extract canonical viral regions."""

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


TOOL_PRIORITY = {"genomad": 0, "checkv": 1, "cenotetaker3": 2, "vicat": 3}

REQUIRED_EVIDENCE_COLUMNS = {
    "sample_id",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "coordinates",
    "tool",
    "classification",
}

MAP_COLUMNS = [
    "sample_id",
    "input_type",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "coordinates",
    "original_length",
    "refined_length",
    "boundary_source",
    "supporting_boundary_tools",
    "boundary_status",
]

AUDIT_COLUMNS = [
    "sample_id",
    "parent_sequence_id",
    "locus_id",
    "tool",
    "call_sequence_id",
    "start",
    "end",
    "length",
    "selected",
    "selected_tool",
    "selected_start",
    "selected_end",
    "supporting_boundary_tools",
    "boundary_status",
]

SUMMARY_COLUMNS = [
    "sample_id",
    "input_type",
    "input_candidate_count",
    "unchanged_candidate_count",
    "refined_parent_count",
    "refined_region_count",
    "boundary_call_count",
    "boundary_conflict_count",
    "ct3_only_boundary_call_count",
    "ct3_only_locus_skipped_count",
    "vicat_advisory_boundary_call_count",
    "vicat_only_locus_skipped_count",
    "checkv_only_locus_skipped_count",
    "checkv_vicat_supported_locus_count",
    "vicat_support_min_overlap_fraction",
    "allow_ct3_only_refinement",
    "evidence_file_count",
]


@dataclass(frozen=True)
class BoundaryCall:
    parent_id: str
    tool: str
    sequence_id: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select geNomad, CheckV, then Cenote-Taker 3 provirus boundaries "
            "and write the refined FASTA used by taxonomy stages."
        )
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=("dna", "rna"))
    parser.add_argument("--candidate-fasta", required=True, type=Path)
    parser.add_argument("--evidence", nargs="*", default=[], type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--output-map", required=True, type=Path)
    parser.add_argument("--output-audit", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    parser.add_argument(
        "--allow-ct3-only-refinement",
        "--allow_ct3_only_refinement",
        action="store_true",
        help=(
            "Allow Cenote-Taker 3 boundaries to modify the FASTA without an "
            "overlapping geNomad or CheckV boundary. Disabled by default; CT3-only "
            "calls remain in the boundary audit while the parent stays unchanged."
        ),
    )
    parser.add_argument(
        "--vicat-support-min-overlap-fraction",
        "--vicat_support_min_overlap_fraction",
        type=float,
        default=0.5,
        help=(
            "Minimum fraction of a viCAT viral-cluster span that must overlap a "
            "CheckV region before viCAT can corroborate CheckV trimming [0.5]."
        ),
    )
    return parser.parse_args()


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        return list(reader)


def read_fasta(handle: TextIO) -> Iterator[tuple[str, str]]:
    identifier: str | None = None
    sequence_parts: list[str] = []
    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                yield identifier, "".join(sequence_parts)
            header = line[1:].strip()
            if not header:
                raise ValueError(f"Empty FASTA header at line {line_number}")
            identifier = header.split(maxsplit=1)[0]
            sequence_parts = []
        else:
            if identifier is None:
                raise ValueError(
                    f"Sequence data precedes the first FASTA header at line {line_number}"
                )
            sequence_parts.append("".join(line.split()))
    if identifier is not None:
        yield identifier, "".join(sequence_parts)


def load_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for identifier, sequence in read_fasta(handle):
            if identifier in records:
                raise ValueError(f"Duplicate candidate FASTA identifier: {identifier}")
            if not sequence:
                raise ValueError(f"Empty candidate FASTA sequence: {identifier}")
            records[identifier] = sequence
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_coordinates(value: str, sequence_id: str) -> tuple[int, int]:
    coordinate = value.strip()
    if coordinate.count("-") != 1:
        raise ValueError(f"Invalid provirus coordinates for {sequence_id}: {value}")
    start_text, end_text = coordinate.split("-", 1)
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError as error:
        raise ValueError(
            f"Non-integer provirus coordinates for {sequence_id}: {value}"
        ) from error
    if start < 1 or end < start:
        raise ValueError(f"Invalid provirus coordinates for {sequence_id}: {value}")
    return start, end


def load_boundary_calls(
    paths: list[Path], sample_id: str, fasta_records: dict[str, str]
) -> list[BoundaryCall]:
    calls: list[BoundaryCall] = []
    seen_paths: set[Path] = set()
    seen_calls: set[tuple[str, str, int, int]] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen_paths:
            raise ValueError(f"Evidence file supplied more than once: {path}")
        seen_paths.add(resolved)
        for row in read_tsv(path, REQUIRED_EVIDENCE_COLUMNS):
            if row["sample_id"] != sample_id:
                raise ValueError(
                    f"Evidence sample '{row['sample_id']}' in {path} does not match "
                    f"'{sample_id}'"
                )
            if row["classification"].strip() != "virus":
                continue
            if row["record_type"].strip() != "provirus":
                continue
            tool = row["tool"].strip().lower()
            if tool not in TOOL_PRIORITY:
                continue
            parent_id = row["parent_sequence_id"].strip()
            sequence_id = row["sequence_id"].strip()
            if not parent_id or not sequence_id:
                raise ValueError(f"Provirus evidence lacks parent/child identity in {path}")
            if parent_id not in fasta_records:
                raise ValueError(
                    f"Provirus parent is absent from candidate FASTA: {parent_id}"
                )
            start, end = parse_coordinates(row["coordinates"], sequence_id)
            if end > len(fasta_records[parent_id]):
                raise ValueError(
                    f"Provirus boundary exceeds parent {parent_id}: {start}-{end}"
                )
            key = (parent_id, tool, start, end)
            if key in seen_calls:
                raise ValueError(
                    f"Duplicate {tool} provirus boundary for {parent_id}: {start}-{end}"
                )
            seen_calls.add(key)
            calls.append(BoundaryCall(parent_id, tool, sequence_id, start, end))
    return calls


def intervals_overlap(left: BoundaryCall, right: BoundaryCall) -> bool:
    return left.start <= right.end and right.start <= left.end


def overlap_fraction_of_support(
    boundary: BoundaryCall, support: BoundaryCall
) -> float:
    overlap = max(0, min(boundary.end, support.end) - max(boundary.start, support.start) + 1)
    return overlap / support.length


def corroborating_tools(
    locus: list[BoundaryCall],
    selected: BoundaryCall | None,
    vicat_support_min_overlap_fraction: float,
) -> list[str]:
    tools = {call.tool for call in locus if call.tool != "vicat"}
    if selected is None:
        tools.update(call.tool for call in locus if call.tool == "vicat")
    elif any(
        overlap_fraction_of_support(selected, call)
        >= vicat_support_min_overlap_fraction
        for call in locus
        if call.tool == "vicat"
    ):
        tools.add("vicat")
    return sorted(tools)


def group_loci(calls: list[BoundaryCall]) -> list[list[BoundaryCall]]:
    """Group transitively overlapping boundary calls into independent loci."""
    remaining = sorted(calls, key=lambda call: (call.start, call.end, call.tool))
    loci: list[list[BoundaryCall]] = []
    while remaining:
        locus = [remaining.pop(0)]
        changed = True
        while changed:
            changed = False
            kept: list[BoundaryCall] = []
            for candidate in remaining:
                if any(intervals_overlap(candidate, member) for member in locus):
                    locus.append(candidate)
                    changed = True
                else:
                    kept.append(candidate)
            remaining = kept
        loci.append(sorted(locus, key=lambda call: (call.start, call.end, call.tool)))
    return loci


def unique_call_for_tool(locus: list[BoundaryCall], tool: str) -> BoundaryCall:
    eligible = [call for call in locus if call.tool == tool]
    if len(eligible) != 1:
        details = ", ".join(f"{call.tool}:{call.start}-{call.end}" for call in eligible)
        raise ValueError(
            "A single tool reported overlapping provirus calls for one locus; "
            f"manual disambiguation is required: {details}"
        )
    return eligible[0]


def select_boundary(
    locus: list[BoundaryCall],
    allow_ct3_only_refinement: bool,
    vicat_support_min_overlap_fraction: float,
) -> BoundaryCall | None:
    tools = {call.tool for call in locus}
    if "genomad" in tools:
        return unique_call_for_tool(locus, "genomad")
    if "checkv" in tools and "cenotetaker3" in tools:
        return unique_call_for_tool(locus, "checkv")
    if "cenotetaker3" in tools:
        return (
            unique_call_for_tool(locus, "cenotetaker3")
            if allow_ct3_only_refinement else None
        )
    if "checkv" in tools and "vicat" in tools:
        checkv = unique_call_for_tool(locus, "checkv")
        if any(
            overlap_fraction_of_support(checkv, call)
            >= vicat_support_min_overlap_fraction
            for call in locus
            if call.tool == "vicat"
        ):
            return checkv
    return None


def run(args: argparse.Namespace) -> None:
    if not 0.0 <= args.vicat_support_min_overlap_fraction <= 1.0:
        raise ValueError("viCAT support overlap fraction must be between 0 and 1")
    fasta_records = load_fasta(args.candidate_fasta)
    calls = load_boundary_calls(args.evidence, args.sample_id, fasta_records)
    calls_by_parent: dict[str, list[BoundaryCall]] = defaultdict(list)
    for call in calls:
        calls_by_parent[call.parent_id].append(call)

    refined_fasta: list[tuple[str, str]] = []
    map_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    refined_parent_count = 0
    refined_region_count = 0
    conflict_count = 0
    unchanged_count = 0
    ct3_only_boundary_call_count = 0
    ct3_only_locus_skipped_count = 0
    vicat_advisory_boundary_call_count = sum(call.tool == "vicat" for call in calls)
    vicat_only_locus_skipped_count = 0
    checkv_only_locus_skipped_count = 0
    checkv_vicat_supported_locus_count = 0

    for parent_id, sequence in fasta_records.items():
        parent_calls = calls_by_parent.get(parent_id, [])
        if not parent_calls:
            unchanged_count += 1
            refined_fasta.append((parent_id, sequence))
            map_rows.append(
                {
                    "sample_id": args.sample_id,
                    "input_type": args.input_type,
                    "sequence_id": parent_id,
                    "parent_sequence_id": "",
                    "record_type": "input_contig",
                    "coordinates": "",
                    "original_length": len(sequence),
                    "refined_length": len(sequence),
                    "boundary_source": "",
                    "supporting_boundary_tools": "",
                    "boundary_status": "unchanged_no_provirus_call",
                }
            )
            continue

        selected_loci: list[
            tuple[int, list[BoundaryCall], BoundaryCall, list[str], str]
        ] = []
        for locus_index, locus in enumerate(group_loci(parent_calls), start=1):
            selected = select_boundary(
                locus,
                args.allow_ct3_only_refinement,
                args.vicat_support_min_overlap_fraction,
            )
            locus_tools = sorted({call.tool for call in locus})
            supporting_tools = corroborating_tools(
                locus, selected, args.vicat_support_min_overlap_fraction
            )
            distinct_boundaries = {(call.start, call.end) for call in locus}
            if locus_tools == ["cenotetaker3"]:
                ct3_only_boundary_call_count += len(locus)
            if locus_tools == ["vicat"]:
                boundary_status = "not_selected_vicat_advisory_only"
                vicat_only_locus_skipped_count += 1
            elif selected is None and locus_tools == ["checkv"]:
                boundary_status = "not_selected_checkv_only_candidate"
                checkv_only_locus_skipped_count += 1
            elif selected is None and set(locus_tools) == {"checkv", "vicat"}:
                boundary_status = "not_selected_vicat_overlap_below_threshold"
                checkv_only_locus_skipped_count += 1
            elif selected is None:
                boundary_status = "not_selected_ct3_only_default"
                ct3_only_locus_skipped_count += 1
            elif set(locus_tools) == {"checkv", "vicat"}:
                boundary_status = "selected_checkv_with_vicat_support"
                checkv_vicat_supported_locus_count += 1
            elif len(locus) == 1:
                boundary_status = "selected_single_tool"
            elif len(distinct_boundaries) == 1:
                boundary_status = "selected_exact_agreement"
            else:
                boundary_status = "selected_boundary_conflict"
                conflict_count += 1

            locus_id = f"{parent_id}:locus_{locus_index}"
            for call in locus:
                audit_rows.append(
                    {
                        "sample_id": args.sample_id,
                        "parent_sequence_id": parent_id,
                        "locus_id": locus_id,
                        "tool": call.tool,
                        "call_sequence_id": call.sequence_id,
                        "start": call.start,
                        "end": call.end,
                        "length": call.length,
                        "selected": "true" if selected and call == selected else "false",
                        "selected_tool": selected.tool if selected else "",
                        "selected_start": selected.start if selected else "",
                        "selected_end": selected.end if selected else "",
                        "supporting_boundary_tools": ",".join(supporting_tools),
                        "boundary_status": boundary_status,
                    }
                )
            if selected is not None:
                selected_loci.append(
                    (locus_index, locus, selected, supporting_tools, boundary_status)
                )

        if not selected_loci:
            supporting_tools = sorted({call.tool for call in parent_calls})
            unchanged_count += 1
            refined_fasta.append((parent_id, sequence))
            map_rows.append(
                {
                    "sample_id": args.sample_id,
                    "input_type": args.input_type,
                    "sequence_id": parent_id,
                    "parent_sequence_id": "",
                    "record_type": "input_contig",
                    "coordinates": "",
                    "original_length": len(sequence),
                    "refined_length": len(sequence),
                    "boundary_source": "",
                    "supporting_boundary_tools": ",".join(supporting_tools),
                    "boundary_status": (
                        "unchanged_vicat_advisory_only"
                        if supporting_tools == ["vicat"]
                        else "unchanged_checkv_candidate_not_corroborated"
                        if "checkv" in supporting_tools
                        else "unchanged_ct3_only_boundary_not_allowed"
                    ),
                }
            )
            continue

        refined_parent_count += 1
        for _, _, selected, supporting_tools, boundary_status in selected_loci:
            refined_sequence = sequence[selected.start - 1 : selected.end]
            output_id = f"{parent_id}|viral_region_{selected.start}_{selected.end}"
            refined_fasta.append((output_id, refined_sequence))
            refined_region_count += 1
            map_rows.append(
                {
                    "sample_id": args.sample_id,
                    "input_type": args.input_type,
                    "sequence_id": output_id,
                    "parent_sequence_id": parent_id,
                    "record_type": (
                        "provirus" if args.input_type == "dna" else "viral_region"
                    ),
                    "coordinates": f"{selected.start}-{selected.end}",
                    "original_length": len(sequence),
                    "refined_length": len(refined_sequence),
                    "boundary_source": selected.tool,
                    "supporting_boundary_tools": ",".join(supporting_tools),
                    "boundary_status": boundary_status,
                }
            )

    if len({identifier for identifier, _ in refined_fasta}) != len(refined_fasta):
        raise ValueError("Refined FASTA identifiers are not unique")

    write_fasta(args.output_fasta, refined_fasta)
    write_tsv(args.output_map, MAP_COLUMNS, map_rows)
    write_tsv(args.output_audit, AUDIT_COLUMNS, audit_rows)
    write_tsv(
        args.output_summary,
        SUMMARY_COLUMNS,
        [
            {
                "sample_id": args.sample_id,
                "input_type": args.input_type,
                "input_candidate_count": len(fasta_records),
                "unchanged_candidate_count": unchanged_count,
                "refined_parent_count": refined_parent_count,
                "refined_region_count": refined_region_count,
                "boundary_call_count": len(calls),
                "boundary_conflict_count": conflict_count,
                "ct3_only_boundary_call_count": ct3_only_boundary_call_count,
                "ct3_only_locus_skipped_count": ct3_only_locus_skipped_count,
                "vicat_advisory_boundary_call_count": vicat_advisory_boundary_call_count,
                "vicat_only_locus_skipped_count": vicat_only_locus_skipped_count,
                "checkv_only_locus_skipped_count": checkv_only_locus_skipped_count,
                "checkv_vicat_supported_locus_count": checkv_vicat_supported_locus_count,
                "vicat_support_min_overlap_fraction": args.vicat_support_min_overlap_fraction,
                "allow_ct3_only_refinement": str(
                    args.allow_ct3_only_refinement
                ).lower(),
                "evidence_file_count": len(args.evidence),
            }
        ],
    )
    print(
        f"Provirus refinement sample={args.sample_id} candidates={len(fasta_records)} "
        f"refined_parents={refined_parent_count} regions={refined_region_count}"
    )


def main() -> None:
    args = parse_args()
    try:
        run(args)
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
