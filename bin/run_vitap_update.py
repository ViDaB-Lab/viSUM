#!/usr/bin/env python3
"""Run VITAP upd with a compatibility fix for repeated coordinate accessions.

VITAP 1.12 deletes the downloaded parent FASTA after extracting the first
coordinate-qualified region. Current ICTV VMR releases can describe multiple
segments on the same parent accession (for example, AH009795), so the second
region fails and the shared ``*_segment.fasta`` name would overwrite regions.

This wrapper patches the installed updater in memory. It does not modify the
Conda environment or VITAP installation.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace


DOWNLOAD_START = "    downloaded_ids = {"
DOWNLOAD_END = "        progress_bar.close()"
REGION_START = "    # Processing integrated viral sequences"
REGION_END = '    print("[INFO] All files successfully downloaded and processed")'


DOWNLOAD_REPLACEMENT = '''    downloaded_ids = {
        f[:-6] for f in os.listdir(output_folder) if f.endswith(".fasta")
    }
    # A VMR can list several coordinate-defined regions from one accession.
    # Download each parent accession once; VITAP 1.12 otherwise schedules
    # concurrent writes to the same FASTA path.
    unique_download_rows = {}
    for row in rows:
        unique_download_rows.setdefault(row[0], row)
    with ThreadPoolExecutor(max_workers=3) as executor:
        progress_bar = tqdm(
            total=len(unique_download_rows), desc="Downloading genomes"
        )
        futures = [
            executor.submit(
                download_and_process_genome,
                row,
                output_folder,
                downloaded_ids,
                progress_bar,
                counter_lock,
            )
            for row in unique_download_rows.values()
        ]
        for future in as_completed(futures):
            future.result()
        progress_bar.close()'''


REGION_REPLACEMENT = '''    # Processing integrated viral sequences.
    # VITAP 1.12 removes the parent after every region and writes every region
    # to the same filename. Group regions by parent, retain the parent until all
    # regions are extracted, and use IDs whose prefix still maps to the VMR
    # accession in VITAP's downstream ``split(".", 1)[0]`` normalization.
    coordinate_regions = defaultdict(list)
    for row in rows:
        virus_id = row[0]
        start_end_sites = row[-1]
        if start_end_sites != "full_length":
            coordinate_regions[virus_id].append(start_end_sites)

    for virus_id, regions in coordinate_regions.items():
        input_fasta = os.path.join(output_folder, f"{virus_id}.fasta")
        if not os.path.isfile(input_fasta) or os.path.getsize(input_fasta) == 0:
            raise FileNotFoundError(
                f"Downloaded parent FASTA is unavailable for coordinate "
                f"extraction: {input_fasta}"
            )

        records = list(SeqIO.parse(input_fasta, "fasta"))
        if len(records) != 1:
            raise ValueError(
                f"Expected one record in {input_fasta}; found {len(records)}"
            )
        parent_record = records[0]
        parent_length = len(parent_record.seq)

        # Remove a partial file left by an earlier unpatched VITAP 1.12 run.
        legacy_segment = os.path.join(output_folder, f"{virus_id}_segment.fasta")
        if os.path.exists(legacy_segment):
            os.remove(legacy_segment)

        for segment_index, start_end_sites in enumerate(regions, start=1):
            start, end = map(int, start_end_sites.split("~"))
            if start < 1 or end < start or end > parent_length:
                raise ValueError(
                    f"Invalid VMR coordinates for {virus_id}: {start}..{end} "
                    f"(parent length {parent_length})"
                )
            segment_id = f"{virus_id}.segment{segment_index}"
            output_fasta = os.path.join(output_folder, f"{segment_id}.fasta")
            segment_record = SeqRecord(
                parent_record.seq[start - 1:end],
                id=segment_id,
                description=(
                    f"parent={virus_id} coordinates={start}..{end}"
                ),
            )
            SeqIO.write([segment_record], output_fasta, "fasta")

        os.remove(input_fasta)
        for fai in glob.glob(os.path.join(output_folder, "*.fai")):
            os.remove(fai)

    # The managed VMR_Genome directory is a download cache shared across VMR
    # releases. Remove retired accessions before VITAP merges every cached
    # FASTA into the new release database.
    expected_fasta_names = {
        f"{row[0]}.fasta"
        for row in rows
        if row[-1] == "full_length"
    }
    for virus_id, regions in coordinate_regions.items():
        expected_fasta_names.update(
            f"{virus_id}.segment{index}.fasta"
            for index in range(1, len(regions) + 1)
        )
    cached_fasta_paths = glob.glob(os.path.join(output_folder, "*.fasta"))
    cached_fasta_names = {os.path.basename(path) for path in cached_fasta_paths}
    missing_fasta_names = sorted(expected_fasta_names - cached_fasta_names)
    if missing_fasta_names:
        examples = ", ".join(missing_fasta_names[:10])
        raise FileNotFoundError(
            f"VITAP genome download is incomplete; missing "
            f"{len(missing_fasta_names)} expected FASTA file(s): {examples}"
        )
    stale_fasta_paths = [
        path
        for path in cached_fasta_paths
        if os.path.basename(path) not in expected_fasta_names
    ]
    for stale_path in stale_fasta_paths:
        os.remove(stale_path)
    if stale_fasta_paths:
        print(
            f"[viSUM] Removed {len(stale_fasta_paths)} stale genome FASTA "
            f"file(s) not present in the current VMR."
        )
    print("[INFO] All files successfully downloaded and processed")'''


def replace_section(source: str, start: str, end: str, replacement: str) -> str:
    start_index = source.find(start)
    if start_index < 0:
        raise RuntimeError(f"VITAP compatibility start marker not found: {start}")
    end_index = source.find(end, start_index)
    if end_index < 0:
        raise RuntimeError(f"VITAP compatibility end marker not found: {end}")
    end_index += len(end)
    return source[:start_index] + replacement + source[end_index:]


def patch_updater_source(source: str) -> str:
    patched = replace_section(
        source, DOWNLOAD_START, DOWNLOAD_END, DOWNLOAD_REPLACEMENT
    )
    return replace_section(
        patched, REGION_START, REGION_END, REGION_REPLACEMENT
    )


def load_patched_updater() -> ModuleType:
    vitap_cli = shutil.which("VITAP")
    cli_sibling = (
        Path(vitap_cli).resolve().with_name("VITAP_upd.py")
        if vitap_cli
        else None
    )
    if cli_sibling is not None and cli_sibling.is_file():
        source_path = cli_sibling
    else:
        spec = importlib.util.find_spec("VITAP_upd")
        if spec is None or spec.origin is None:
            raise RuntimeError("Could not locate the installed VITAP_upd module.")
        source_path = Path(spec.origin)
    # Bioconda installs VITAP_upd.py and its helper module beside the VITAP
    # executable rather than in site-packages.
    sys.path.insert(0, str(source_path.parent))
    source = source_path.read_text(encoding="utf-8")
    patched_source = patch_updater_source(source)
    module = ModuleType("VITAP_upd_visum_compat")
    module.__file__ = str(source_path)
    exec(compile(patched_source, str(source_path), "exec"), module.__dict__)
    return module


def configure_diamond_threads(module: ModuleType, threads: int) -> None:
    """Inject the Nextflow CPU allocation into every DIAMOND subprocess."""
    if threads < 1:
        raise ValueError("--threads must be at least 1")
    original_run = module.subprocess.run

    def controlled_run(command, *args, **kwargs):
        if (
            isinstance(command, (list, tuple))
            and len(command) >= 2
            and Path(str(command[0])).name == "diamond"
            and "--threads" not in command
            and "-p" not in command
        ):
            command = [
                command[0],
                command[1],
                "--threads",
                str(threads),
                *command[2:],
            ]
            print(f"[viSUM] DIAMOND command limited to {threads} thread(s).")
        return original_run(command, *args, **kwargs)

    module.subprocess.run = controlled_run


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vmr", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--threads", required=True, type=int)
    args = parser.parse_args()

    updater = load_patched_updater()
    configure_diamond_threads(updater, args.threads)
    updater.upd(SimpleNamespace(vmr=args.vmr, out=args.out, db=args.db))


if __name__ == "__main__":
    main()
