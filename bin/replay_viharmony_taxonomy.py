#!/usr/bin/env python3

"""Replay viHARMONY taxonomy policy from cached evidence-audit tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path

import run_viharmony


OUTPUT_COLUMNS = [
    "input_type", "split", "original_contig_name", "final_sequence_id",
    "strict_taxonomy", "strict_taxonomy_rank", "strict_species_correct",
    "exploratory_taxonomy", "exploratory_taxonomy_rank",
    "exploratory_species_correct",
    "species_assignment_status", "species_assignment_basis",
    "strict_taxonomy_stop_reason", "strict_taxonomy_conflict_rank",
    "strict_taxonomy_conflicting_taxa",
    *[
        f"{scope}_{rank}"
        for rank in run_viharmony.RANKS[1:]
        for scope in ("truth", "strict", "exploratory")
    ],
]

if len(OUTPUT_COLUMNS) != len(set(OUTPUT_COLUMNS)):
    raise RuntimeError("Taxonomy replay output schema contains duplicate columns")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def heldout_group(row: dict[str, str]) -> str:
    for column in ("family", "order", "class", "parent_virus_id"):
        value = row.get(column, "").strip()
        if value and value.lower() != "unclassified":
            return f"{column}:{value}"
    return f"sequence:{row['benchmark_sequence_id']}"


def split_for(row: dict[str, str], modulus: int) -> str:
    digest = hashlib.sha256(heldout_group(row).encode("utf-8")).digest()
    return "heldout" if int.from_bytes(digest[:8], "big") % modulus == 0 else "development"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def replay(
    results_root: Path,
    truth_path: Path,
    msl_path: Path,
    output_directory: Path,
    holdout_modulus: int,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    if holdout_modulus < 2:
        raise ValueError("Holdout modulus must be at least 2")
    truth_rows = read_tsv(truth_path)
    truth_by_id = {
        row["benchmark_sequence_id"]: row for row in truth_rows
    }
    if len(truth_by_id) != len(truth_rows):
        raise ValueError("Truth table contains duplicate benchmark_sequence_id values")
    msl, _ = run_viharmony.load_msl(msl_path)

    replay_rows: list[dict[str, object]] = []
    for metadata_path in sorted(results_root.rglob("*.final_metadata.tsv")):
        stem = metadata_path.name.removesuffix(".final_metadata.tsv")
        audit_path = metadata_path.with_name(f"{stem}.evidence_audit.tsv.gz")
        if not audit_path.exists():
            raise ValueError(f"Missing cached evidence audit: {audit_path}")
        evidence_by_final: dict[str, list[dict[str, str]]] = defaultdict(list)
        for evidence in read_gzip_tsv(audit_path):
            if evidence.get("evidence_application") == "used_for_final_decision":
                evidence_by_final[evidence["final_sequence_id"]].append(evidence)

        for metadata in read_tsv(metadata_path):
            original_id = metadata["original_contig_name"]
            truth = truth_by_id.get(original_id)
            if truth is None or truth.get("truth_label", "").lower() != "viral":
                continue
            final_id = metadata["final_sequence_id"]
            decision = run_viharmony.decide_taxonomy(
                evidence_by_final.get(final_id, []), msl
            )
            strict_species = decision.strict.get("species", "")
            exploratory_species = decision.exploratory.get("species", "")
            truth_species = truth.get("species", "")
            replay_row: dict[str, object] = {
                "input_type": truth.get("input_type", ""),
                "split": split_for(truth, holdout_modulus),
                "original_contig_name": original_id,
                "final_sequence_id": final_id,
                "truth_family": truth.get("family", ""),
                "truth_genus": truth.get("genus", ""),
                "truth_species": truth_species,
                "strict_taxonomy": run_viharmony.taxonomy_string(decision.strict),
                "strict_taxonomy_rank": decision.strict_rank,
                "strict_species": strict_species or "NA",
                "strict_species_correct": str(strict_species == truth_species).lower() if strict_species else "NA",
                "exploratory_taxonomy": run_viharmony.taxonomy_string(decision.exploratory),
                "exploratory_taxonomy_rank": decision.exploratory_rank,
                "exploratory_species": exploratory_species or "NA",
                "exploratory_species_correct": str(exploratory_species == truth_species).lower() if exploratory_species else "NA",
                "species_assignment_status": decision.species_status,
                "species_assignment_basis": decision.species_basis,
                "strict_taxonomy_stop_reason": decision.strict_stop_reason,
                "strict_taxonomy_conflict_rank": decision.conflict_rank,
                "strict_taxonomy_conflicting_taxa": ",".join(decision.conflicting_taxa) or "NA",
            }
            for rank in run_viharmony.RANKS[1:]:
                truth_taxon = truth.get(rank, "")
                if truth_taxon.lower() == "unclassified":
                    truth_taxon = ""
                replay_row[f"truth_{rank}"] = truth_taxon or "NA"
                replay_row[f"strict_{rank}"] = decision.strict.get(rank, "") or "NA"
                replay_row[f"exploratory_{rank}"] = (
                    decision.exploratory.get(rank, "") or "NA"
                )
            replay_rows.append(replay_row)

    summary_rows: list[dict[str, object]] = []
    for input_type in sorted({str(row["input_type"]) for row in replay_rows}):
        for split in ("all", "development", "heldout"):
            selected = [
                row for row in replay_rows
                if row["input_type"] == input_type
                and (split == "all" or row["split"] == split)
            ]
            for mode in ("strict", "exploratory"):
                called = [row for row in selected if row[f"{mode}_species"] != "NA"]
                correct = [
                    row for row in called
                    if row[f"{mode}_species_correct"] == "true"
                ]
                summary_rows.append({
                    "input_type": input_type,
                    "split": split,
                    "mode": mode,
                    "eligible": len(selected),
                    "called": len(called),
                    "correct": len(correct),
                    "incorrect": len(called) - len(correct),
                    "coverage": f"{len(called) / len(selected):.6f}" if selected else "NA",
                    "precision": f"{len(correct) / len(called):.6f}" if called else "NA",
                })

    rank_summary_rows: list[dict[str, object]] = []
    for input_type in sorted({str(row["input_type"]) for row in replay_rows}):
        for split in ("all", "development", "heldout"):
            selected = [
                row for row in replay_rows
                if row["input_type"] == input_type
                and (split == "all" or row["split"] == split)
            ]
            for rank in run_viharmony.RANKS[1:]:
                eligible = [row for row in selected if row[f"truth_{rank}"] != "NA"]
                for mode in ("strict", "exploratory"):
                    called = [
                        row for row in eligible if row[f"{mode}_{rank}"] != "NA"
                    ]
                    correct = [
                        row for row in called
                        if row[f"{mode}_{rank}"] == row[f"truth_{rank}"]
                    ]
                    rank_summary_rows.append({
                        "input_type": input_type,
                        "split": split,
                        "rank": rank,
                        "mode": mode,
                        "eligible": len(eligible),
                        "called": len(called),
                        "correct": len(correct),
                        "incorrect": len(called) - len(correct),
                        "coverage": f"{len(called) / len(eligible):.6f}" if eligible else "NA",
                        "precision": f"{len(correct) / len(called):.6f}" if called else "NA",
                    })

    output_directory.mkdir(parents=True, exist_ok=True)
    write_tsv(output_directory / "taxonomy_policy_replay.tsv", OUTPUT_COLUMNS, replay_rows)
    write_tsv(
        output_directory / "taxonomy_policy_summary.tsv",
        ["input_type", "split", "mode", "eligible", "called", "correct", "incorrect", "coverage", "precision"],
        summary_rows,
    )
    write_tsv(
        output_directory / "taxonomy_rank_summary.tsv",
        ["input_type", "split", "rank", "mode", "eligible", "called", "correct", "incorrect", "coverage", "precision"],
        rank_summary_rows,
    )
    return replay_rows, summary_rows, rank_summary_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    parser.add_argument("--ictv-msl", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--holdout-modulus", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        rows, summary, rank_summary = replay(
            args.results_root, args.truth, args.ictv_msl,
            args.output_directory, args.holdout_modulus,
        )
    except (OSError, ValueError, csv.Error) as error:
        raise SystemExit(f"ERROR: {error}") from error
    print(
        f"Taxonomy replay complete: rows={len(rows)} "
        f"species_summaries={len(summary)} rank_summaries={len(rank_summary)}"
    )


if __name__ == "__main__":
    main()
