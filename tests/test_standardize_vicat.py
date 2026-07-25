import csv
import os
import subprocess
import sys
from pathlib import Path

import duckdb


ROOT = Path(__file__).resolve().parents[1]
TAXONOMY = [
    "d__Domain", "r__Realm", "k__Kingdom", "p__Phylum", "c__Class",
    "o__Order", "f__Family", "g__Genus", "s__Species",
]


def write_tsv(path: Path, columns: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def make_lookup(path: Path) -> None:
    connection = duckdb.connect()
    connection.execute(
        f"""COPY (SELECT * FROM (VALUES
        ('repA', 2, 2, 1.0, false, '', 'geNomad',
         'd__Viruses','r__R','k__K','p__P','c__C','o__O','f__A','g__A','s__A'),
        ('repB', 3, 3, 1.0, false, '', 'geNomad',
         'd__Viruses','r__R','k__K','p__P','c__C','o__O','f__B','g__B','s__B'),
        ('repC', 1, 1, 1.0, false, '', 'geNomad',
         'd__Viruses','r__R','k__K','p__P','c__C','o__O','f__A','g__A','s__A')
        ) t(representative_protein_id, member_votu_count, classified_votu_count,
            taxonomy_coverage, taxonomy_conflict, taxonomy_conflict_rank,
            taxonomy_methods, "d__Domain", "r__Realm", "k__Kingdom", "p__Phylum",
            "c__Class", "o__Order", "f__Family", "g__Genus", "s__Species"))
        TO '{str(path).replace("'", "''")}' (FORMAT PARQUET)"""
    )
    connection.close()


def test_locus_voting_and_zero_hit_contig(tmp_path: Path) -> None:
    orf_map = tmp_path / "orfs.tsv"
    write_tsv(
        orf_map,
        ["sample_id", "orf_id", "sequence_id", "start", "end", "strand", "protein_length", "translation_table", "callers"],
        [
            ["sample", "orfA", "sample_c000001", 1, 300, "+", 100, 11, "pyrodigal-gv"],
            ["sample", "orfB", "sample_c000001", 1, 300, "+", 100, 11, "pyrodigal-rv"],
            ["sample", "orfC", "sample_c000001", 400, 700, "+", 100, 11, "pyrodigal-gv"],
            ["sample", "orfD", "sample_c000002", 1, 300, "+", 100, 11, "pyrodigal-gv"],
        ],
    )
    header_map = tmp_path / "headers.tsv"
    write_tsv(header_map, ["sequence_id", "length"], [["sample_c000001", 1000], ["sample_c000002", 500]])
    diamond = tmp_path / "diamond.tsv"
    columns = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp"]
    write_tsv(
        diamond, columns,
        [
            ["orfA", "repA", 80, 100, 100, 120, 1, 100, 1, 100, "1e-30", 100, 100, 83],
            ["orfB", "repB", 79, 100, 100, 120, 1, 100, 1, 100, "1e-29", 96, 100, 83],
            ["orfC", "repC", 75, 90, 100, 110, 1, 90, 1, 90, "1e-20", 80, 90, 82],
        ],
    )
    lookup = tmp_path / "lookup.parquet"
    make_lookup(lookup)
    loci = tmp_path / "loci.tsv"
    evidence = tmp_path / "evidence.tsv"
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(ROOT / "bin") + os.pathsep + environment.get("PYTHONPATH", "")
    subprocess.run(
        [
            sys.executable, str(ROOT / "bin" / "standardize_vicat.py"),
            "--sample-id", "sample", "--input-type", "rna",
            "--orf-map", str(orf_map), "--diamond", str(diamond),
            "--taxonomy-lookup", str(lookup), "--header-map", str(header_map),
            "--taxonomy-support", "0.60", "--locus-overlap", "0.80",
            "--output-loci", str(loci), "--output-evidence", str(evidence),
        ],
        check=True, env=environment,
    )
    with loci.open(encoding="utf-8", newline="") as handle:
        locus_rows = list(csv.DictReader(handle, delimiter="\t"))
    with evidence.open(encoding="utf-8", newline="") as handle:
        evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(locus_rows) == 3
    assert locus_rows[0]["caller_taxonomy_conflict"] == "true"
    assert len(evidence_rows) == 1
    assert evidence_rows[0]["sequence_id"] == "sample_c000001"
    assert evidence_rows[0]["classification_rank"] == "order"
    assert evidence_rows[0]["taxonomy_conflict"] == "true"
    assert evidence_rows[0]["score"] == "1"
