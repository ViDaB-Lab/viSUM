#!/usr/bin/env python3
"""Compare trimming authority using saved discovery evidence, not stale taxonomy.

Reproduces saved refinement with CT3-only trimming enabled before evaluating
the disabled policy. Reports sequence changes and affected OLD primary parents;
does not claim new final sensitivity/FPR without downstream evidence regeneration.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import refine_proviral_regions as refiner


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def signatures(rows):
    result = {}
    for row in rows:
        parent = row["parent_sequence_id"] or row["sequence_id"]
        result.setdefault(parent, []).append((row["sequence_id"], row["record_type"], row["coordinates"], row["refined_length"]))
    return {k: sorted(v) for k, v in result.items()}


def execute(folder, output, fasta_root):
    sample = folder.name.removesuffix("_results")
    saved = folder / "refinement"
    old_map = read(saved / f"{sample}.provirus_region_map.tsv")
    parents = set(signatures(old_map))
    dest = output / sample
    dest.mkdir()
    candidate = folder / "discovery_gate" / f"{sample}.discovery_candidates.fasta"
    if not candidate.exists():
        if fasta_root is None:
            raise FileNotFoundError(f"Missing candidates; supply --fasta-root: {candidate}")
        normalized = fasta_root / folder.name / "prep" / f"{sample}.normalized.fasta"
        sequences = refiner.load_fasta(normalized)
        if parents - sequences.keys():
            raise ValueError(f"Missing normalized parents: {sample}")
        candidate = dest / "reconstructed_discovery_candidates.fasta"
        with candidate.open("w", encoding="utf-8", newline="\n") as handle:
            for parent in sorted(parents):
                handle.write(f">{parent}\n{sequences[parent]}\n")
    sequences = refiner.load_fasta(candidate)
    if set(sequences) != parents:
        raise ValueError(f"Discovery input/old map mismatch: {sample}")
    evidence = [folder / tool / f"{sample}.{tool}_evidence.tsv"
                for tool in ("genomad", "virsorter2", "cenotetaker3", "checkv", "vicat")]
    evidence.append(folder / "vicat" / f"{sample}.vicat_provirus_evidence.tsv")
    for path in evidence:
        if not path.exists():
            raise FileNotFoundError(path)
    outputs = {}
    for enabled in (True, False):
        variant = "ct3_on" if enabled else "ct3_off"
        target = dest / variant
        target.mkdir()
        args = SimpleNamespace(sample_id=sample, input_type=old_map[0]["input_type"] if old_map else ("rna" if "RNA" in sample else "dna"),
            candidate_fasta=candidate, evidence=evidence, allow_ct3_only_refinement=enabled,
            vicat_support_min_overlap_fraction=0.5, output_fasta=target/"refined.fasta",
            output_map=target/"map.tsv", output_audit=target/"audit.tsv", output_summary=target/"summary.tsv")
        refiner.run(args)
        outputs[variant] = read(args.output_map)
    # Compare semantic table contents, not platform-specific newline bytes.
    if outputs["ct3_on"] != old_map:
        raise ValueError(f"Baseline map mismatch: {sample}; do not interpret this replay")
    baseline_fasta = saved / f"{sample}.refined_candidates.fasta"
    if not baseline_fasta.exists() and fasta_root:
        baseline_fasta = fasta_root / folder.name / "refinement" / f"{sample}.refined_candidates.fasta"
    if not baseline_fasta.exists():
        raise FileNotFoundError(f"Missing baseline FASTA: {baseline_fasta}")
    if refiner.load_fasta(baseline_fasta) != refiner.load_fasta(dest/"ct3_on/refined.fasta"):
        raise ValueError(f"Baseline sequence mismatch: {sample}")
    new_fasta = refiner.load_fasta(dest/"ct3_off/refined.fasta")
    before, after = signatures(old_map), signatures(outputs["ct3_off"])
    if before.keys() != after.keys():
        raise ValueError(f"Parent lost or introduced: {sample}")
    changed = sorted(k for k in before if before[k] != after[k])
    for row in outputs["ct3_off"]:
        if row["boundary_source"] == "cenotetaker3":
            raise ValueError("CT3 still selected with authority disabled")
        parent = row["parent_sequence_id"] or row["sequence_id"]
        expected = sequences[parent]
        if row["coordinates"]:
            a, b = map(int, row["coordinates"].split("-"))
            expected = expected[a-1:b]
        if new_fasta[row["sequence_id"]] != expected:
            raise ValueError("New output is not its declared parent slice")
    retained = read(folder/"viharmony"/f"{sample}.database_candidates.tsv")
    primary = {r["normalized_name"] for r in retained}
    detail = [{"parent": p, "old": before[p], "new": after[p], "was_primary": p in primary} for p in changed]
    (dest/"changes.json").write_text(json.dumps(detail, indent=2), encoding="utf-8")
    result = dict(sample=sample, baseline_map_and_fasta_match=True, discovery_parents=len(parents),
                  changed_parents=len(changed), changed_old_primary_parents=len(primary.intersection(changed)),
                  old_refined_parents=sum(any(v[1] != "input_contig" for v in rows) for rows in before.values()),
                  new_refined_parents=sum(any(v[1] != "input_contig" for v in rows) for rows in after.values()),
                  final_metrics_status="requires_changed_sequence_downstream_rerun")
    (dest/"input_hashes.json").write_text(json.dumps({str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [candidate, *evidence]}, indent=2), encoding="utf-8")
    print(json.dumps(result), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--fasta-root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    folders = sorted(args.results_root.glob("*_results"))
    if not folders:
        parser.error("No sample result folders")
    args.output_dir.mkdir(parents=True, exist_ok=False)
    results = [execute(folder, args.output_dir, args.fasta_root) for folder in folders]
    (args.output_dir/"summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
