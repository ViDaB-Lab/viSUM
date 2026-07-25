import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_prediction(directory: Path, name: str, records: list[tuple]) -> tuple[Path, Path]:
    faa = directory / f"{name}.faa"
    gff = directory / f"{name}.gff"
    with faa.open("w", encoding="utf-8") as proteins, gff.open("w", encoding="utf-8") as features:
        features.write("##gff-version 3\n")
        for feature_id, parent, start, end, strand, translation in records:
            proteins.write(f">{feature_id} # {start} # {end} # 1 # ID={feature_id};partial=00\n{translation}\n")
            features.write(
                f"{parent}\tpyrodigal\tCDS\t{start}\t{end}\t.\t{strand}\t0\t"
                f"ID={feature_id};transl_table=11\n"
            )
    return faa, gff


def test_combines_exact_calls_but_retains_alternatives(tmp_path: Path) -> None:
    gv_faa, gv_gff = write_prediction(
        tmp_path, "gv", [("1_1", "sample_c000001", 1, 90, "+", "MPEPTIDE")]
    )
    rv_faa, rv_gff = write_prediction(
        tmp_path,
        "rv",
        [
            ("1_1", "sample_c000001", 1, 90, "+", "MPEPTIDE"),
            ("1_2", "sample_c000001", 4, 90, "+", "MALT"),
        ],
    )
    output_faa = tmp_path / "combined.faa"
    output_map = tmp_path / "combined.tsv"
    subprocess.run(
        [
            sys.executable, str(ROOT / "bin" / "combine_vicat_orfs.py"),
            "--sample-id", "sample", "--input-type", "rna",
            "--gv-proteins", str(gv_faa), "--gv-gff", str(gv_gff),
            "--rv-proteins", str(rv_faa), "--rv-gff", str(rv_gff),
            "--output-proteins", str(output_faa), "--output-map", str(output_map),
        ],
        check=True,
    )
    with output_map.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 2
    assert rows[0]["callers"] == "pyrodigal-gv,pyrodigal-rv"
    assert rows[1]["callers"] == "pyrodigal-rv"
    assert output_faa.read_text(encoding="utf-8").count(">") == 2


def test_accepts_pyrodigal_gv_short_header_id_and_full_gff_id(tmp_path: Path) -> None:
    faa = tmp_path / "gv.faa"
    gff = tmp_path / "gv.gff"
    faa.write_text(
        ">sample_c000001_1 # 89 # 586 # -1 # ID=1_1;partial=00;start_type=ATG\nMPEPTIDE\n",
        encoding="utf-8",
    )
    gff.write_text(
        "##gff-version 3\n"
        "sample_c000001\tpyrodigal\tCDS\t89\t586\t.\t-\t0\t"
        "ID=sample_c000001_1;partial=00;transl_table=11\n",
        encoding="utf-8",
    )
    output_faa = tmp_path / "combined.faa"
    output_map = tmp_path / "combined.tsv"
    subprocess.run(
        [
            sys.executable, str(ROOT / "bin" / "combine_vicat_orfs.py"),
            "--sample-id", "sample", "--input-type", "dna",
            "--gv-proteins", str(faa), "--gv-gff", str(gff),
            "--output-proteins", str(output_faa), "--output-map", str(output_map),
        ],
        check=True,
    )
    with output_map.open(encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["sequence_id"] == "sample_c000001"
    assert row["start"] == "89"
    assert row["end"] == "586"
    assert row["strand"] == "-"
