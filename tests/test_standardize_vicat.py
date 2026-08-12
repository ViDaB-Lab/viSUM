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
         'd__Viruses','r__R','k__K','p__P','c__C','o__O','f__A','g__A','s__A'),
        ('repRealmConflict', 2, 2, 1.0, true, 'realm', 'geNomad',
         'd__Viruses','r__unclassified','k__unclassified','p__unclassified',
         'c__unclassified','o__unclassified','f__unclassified','g__unclassified','s__unclassified'),
        ('repFamilyConflict', 2, 2, 1.0, true, 'family', 'geNomad',
         'd__Viruses','r__R','k__K','p__P','c__C','o__O',
         'f__unclassified','g__unclassified','s__unclassified')
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
            ["orfA", "repC", 78, 98, 100, 110, 1, 98, 1, 98, "1e-28", 97, 98, 89],
            ["orfB", "repB", 79, 100, 100, 120, 1, 100, 1, 100, "1e-29", 96, 100, 83],
            ["orfC", "repC", 75, 90, 100, 110, 1, 90, 1, 90, "1e-20", 80, 90, 82],
        ],
    )
    lookup = tmp_path / "lookup.parquet"
    make_lookup(lookup)
    loci = tmp_path / "loci.tsv"
    audit = tmp_path / "audit.tsv"
    evidence = tmp_path / "evidence.tsv"
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(ROOT / "bin") + os.pathsep + environment.get("PYTHONPATH", "")
    subprocess.run(
        [
            sys.executable, str(ROOT / "bin" / "standardize_vicat.py"),
            "--sample-id", "sample", "--input-type", "rna",
            "--orf-map", str(orf_map), "--diamond", str(diamond),
            "--taxonomy-lookup", str(lookup), "--header-map", str(header_map),
            "--orf-taxonomy-support", "0.60",
            "--contig-taxonomy-support", "0.60", "--locus-overlap", "0.80",
            "--output-loci", str(loci), "--output-audit", str(audit),
            "--output-evidence", str(evidence),
        ],
        check=True, env=environment,
    )
    with loci.open(encoding="utf-8", newline="") as handle:
        locus_rows = list(csv.DictReader(handle, delimiter="\t"))
    with evidence.open(encoding="utf-8", newline="") as handle:
        evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
    with audit.open(encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(locus_rows) == 3
    assert locus_rows[0]["caller_taxonomy_conflict"] == "true"
    assert len(evidence_rows) == 1
    assert evidence_rows[0]["sequence_id"] == "sample_c000001"
    # The conflicted locus abstains below order; the independent resolved locus
    # remains eligible to supply family through species taxonomy.
    assert evidence_rows[0]["classification_rank"] == "species"
    assert evidence_rows[0]["taxonomy_conflict"] == "true"
    assert evidence_rows[0]["score"] == "1"
    assert len(audit_rows) == 4
    assert {row["reference_id"] for row in audit_rows} == {"repA", "repB", "repC"}
    assert all(row["orf_taxonomy_support_threshold"] == "0.6" for row in audit_rows)
    rep_a = next(row for row in audit_rows if row["reference_id"] == "repA")
    assert rep_a["selected_orf_for_locus"] == "true"
    assert rep_a["selected_best_reference"] == "true"
    assert rep_a["lineage_vote_representative"] == "true"
    assert rep_a["g__Genus"] == "g__A"
    duplicate_lineage = next(
        row for row in audit_rows if row["orf_id"] == "orfA" and row["reference_id"] == "repC"
    )
    assert duplicate_lineage["lineage_vote_representative"] == "false"
    assert duplicate_lineage["lineage_vote_weight"] == ""


def test_conflicted_references_remain_viral_hits_and_abstain_below_safe_rank(tmp_path: Path) -> None:
    orf_map = tmp_path / "orfs.tsv"
    write_tsv(
        orf_map,
        ["sample_id", "orf_id", "sequence_id", "start", "end", "strand", "protein_length", "translation_table", "callers"],
        [
            ["sample", "orfRealmConflict", "sample_c000001", 1, 300, "+", 100, 11, "pyrodigal-gv"],
            ["sample", "orfResolved", "sample_c000001", 400, 700, "+", 100, 11, "pyrodigal-gv"],
            ["sample", "orfFamilyConflict", "sample_c000001", 800, 1100, "+", 100, 11, "pyrodigal-gv"],
        ],
    )
    header_map = tmp_path / "headers.tsv"
    write_tsv(header_map, ["sequence_id", "length"], [["sample_c000001", 1200]])
    diamond = tmp_path / "diamond.tsv"
    columns = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp"]
    write_tsv(
        diamond,
        columns,
        [
            ["orfRealmConflict", "repRealmConflict", 90, 100, 100, 100, 1, 100, 1, 100, "1e-40", 120, 100, 100],
            ["orfResolved", "repA", 85, 100, 100, 100, 1, 100, 1, 100, "1e-30", 100, 100, 100],
            ["orfFamilyConflict", "repFamilyConflict", 88, 100, 100, 100, 1, 100, 1, 100, "1e-35", 110, 100, 100],
        ],
    )
    lookup = tmp_path / "lookup.parquet"
    make_lookup(lookup)
    loci = tmp_path / "loci.tsv"
    audit = tmp_path / "audit.tsv"
    evidence = tmp_path / "evidence.tsv"
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(ROOT / "bin") + os.pathsep + environment.get("PYTHONPATH", "")
    subprocess.run(
        [
            sys.executable, str(ROOT / "bin" / "standardize_vicat.py"),
            "--sample-id", "sample", "--input-type", "rna",
            "--orf-map", str(orf_map), "--diamond", str(diamond),
            "--taxonomy-lookup", str(lookup), "--header-map", str(header_map),
            "--orf-taxonomy-support", "0.60",
            "--contig-taxonomy-support", "0.60", "--locus-overlap", "0.80",
            "--output-loci", str(loci), "--output-audit", str(audit),
            "--output-evidence", str(evidence),
        ],
        check=True, env=environment,
    )

    with loci.open(encoding="utf-8", newline="") as handle:
        locus_rows = {row["locus_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with evidence.open(encoding="utf-8", newline="") as handle:
        evidence_row = next(csv.DictReader(handle, delimiter="\t"))
    with audit.open(encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert len(locus_rows) == 3
    realm_conflict = next(row for row in locus_rows.values() if row["coordinates"] == "1-300")
    assert realm_conflict["has_qualified_hit"] == "true"
    assert realm_conflict["classification_rank"] == "domain"
    assert realm_conflict["r__Realm"] == "r__unclassified"
    realm_audit = next(row for row in audit_rows if row["reference_id"] == "repRealmConflict")
    assert realm_audit["reference_taxonomy_conflict"] == "true"
    assert realm_audit["reference_taxonomy_conflict_rank"] == "realm"

    family_conflict = next(row for row in locus_rows.values() if row["coordinates"] == "800-1100")
    assert family_conflict["classification_rank"] == "order"
    assert family_conflict["o__Order"] == "o__O"
    assert family_conflict["f__Family"] == "f__unclassified"

    # All three loci remain viral-homology hits. At each rank, truncated
    # references abstain, so the safely resolved locus can carry lower ranks.
    assert evidence_row["hit_loci"] == "3"
    assert evidence_row["score"] == "1"
    assert evidence_row["classification"] == "virus"
    assert evidence_row["classification_rank"] == "species"
    assert evidence_row["r__Realm"] == "r__R"
    assert evidence_row["f__Family"] == "f__A"
    assert evidence_row["taxonomy_eligible_loci"] == "1"
    assert evidence_row["taxonomy_supporting_loci"] == "1"
