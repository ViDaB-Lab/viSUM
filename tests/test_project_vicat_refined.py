import csv
import sys
from argparse import Namespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import project_vicat_refined
from evidence_schema import TAXONOMY_COLUMNS
from standardize_vicat import LOCUS_COLUMNS


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_projects_only_contained_loci_and_recalculates_region_taxonomy(tmp_path: Path) -> None:
    loci = tmp_path / "loci.tsv"
    base = {column: "" for column in LOCUS_COLUMNS}
    taxonomy = {column: f"{column[:3]}Taxon" for column in TAXONOMY_COLUMNS}
    rows = []
    for locus_id, coordinates, classification, score in [
        ("cell_left", "1-200", "cellular_supported", ""),
        ("viral_one", "301-500", "viral_supported", "120"),
        ("viral_two", "551-750", "viral_supported", "110"),
        ("crossing", "740-900", "cellular_supported", ""),
    ]:
        row = dict(base)
        row.update({
            "sample_id": "sample", "sequence_id": "parent", "locus_id": locus_id,
            "coordinates": coordinates, "locus_classification": classification,
            "best_bitscore": score, "best_reference_id": "VIRAL|repA" if score else "",
            "classification_rank": "species" if score else "", "orf_callers": "pyrodigal-gv",
            "competitive_mode": "true",
        })
        if score:
            row.update(taxonomy)
        rows.append(row)
    write_tsv(loci, LOCUS_COLUMNS, rows)

    regions = tmp_path / "regions.tsv"
    write_tsv(regions, [
        "sample_id", "sequence_id", "parent_sequence_id", "record_type",
        "coordinates", "refined_length",
    ], [{
        "sample_id": "sample", "sequence_id": "parent|viral_region_301_800",
        "parent_sequence_id": "parent", "record_type": "provirus",
        "coordinates": "301-800", "refined_length": "500",
    }])
    evidence, projection = tmp_path / "evidence.tsv", tmp_path / "projection.tsv"
    project_vicat_refined.run(Namespace(
        sample_id="sample", input_type="dna", loci=loci, region_map=regions,
        contig_taxonomy_support=0.6, cluster_min_viral_loci=2,
        cluster_max_neutral_gap=1, output_evidence=evidence,
        output_projection=projection,
    ))

    evidence_row = read_tsv(evidence)[0]
    assert evidence_row["sequence_id"] == "parent|viral_region_301_800"
    assert evidence_row["parent_sequence_id"] == "parent"
    assert evidence_row["classification"] == "virus"
    assert evidence_row["evidence_scope"] == "refined_region"
    assert evidence_row["viral_supported_loci"] == "2"
    assert evidence_row["classification_rank"] == "species"
    projected = {row["source_locus_id"]: row for row in read_tsv(projection)}
    assert projected["viral_one"]["projected_coordinates"] == "1-200"
    assert projected["viral_two"]["projected_coordinates"] == "251-450"
    assert projected["crossing"]["projection_status"] == "crosses_boundary"
    assert projected["crossing"]["projected_coordinates"] == ""


def test_empty_refined_sample_writes_header_only_outputs(tmp_path: Path) -> None:
    loci = tmp_path / "loci.tsv"
    write_tsv(loci, LOCUS_COLUMNS, [])
    regions = tmp_path / "regions.tsv"
    write_tsv(
        regions,
        [
            "sample_id", "sequence_id", "parent_sequence_id", "record_type",
            "coordinates", "refined_length",
        ],
        [],
    )
    evidence = tmp_path / "evidence.tsv"
    projection = tmp_path / "projection.tsv"

    project_vicat_refined.run(Namespace(
        sample_id="sample", input_type="dna", loci=loci, region_map=regions,
        contig_taxonomy_support=0.6, cluster_min_viral_loci=2,
        cluster_max_neutral_gap=1, output_evidence=evidence,
        output_projection=projection,
    ))

    assert read_tsv(evidence) == []
    assert read_tsv(projection) == []
    assert evidence.read_text(encoding="utf-8").startswith("sample_id\t")
    assert projection.read_text(encoding="utf-8").startswith("sample_id\t")
