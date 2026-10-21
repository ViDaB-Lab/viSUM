#!/usr/bin/env python3
"""Build a standalone offline prototype from the pinned official ICTV export.

No downloads, update checks, pipeline integration, or scientific rule changes.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import gzip
import hashlib
import io
import json
from pathlib import Path
import re

from ictv_reconciliation import RANKS, SCHEMA, sha256


SOURCE_COMMIT = "e83bc2c73a746165a209ba88d76ef8a8fb8635d2"
SOURCE_URL = "https://github.com/ICTV-Virus-Knowledgebase/ICTVdatabase"
SOURCE_BLOBS = {
    "taxonomy_node_export.utf8.txt": "5c8646ed0857bb199a4299dfbfdfab815b25a2e1",
    "taxonomy_node_merge_split.utf8.txt": "072100c76f5eca9fd645eeb55b3023a28c535df0",
    "taxonomy_toc.utf8.txt": "52f54ab2f2ba83ca2f0a4b8361fbee27dfb9c053",
    "taxonomy_level.utf8.txt": "74cb8dcf2677bb901aee0bab41c6bacca52d7d2d",
    "LICENSE.txt": "3f00874ad9f5fb00fed152078981b960961fc4e2",
}


def read_export(path: Path):
    # ICTV's published reader uses backslash escaping and normalizes its \N
    # null sentinel before parsing. Do not silently drop malformed columns.
    text = path.read_text(encoding="utf-8-sig")
    text = re.sub(r"(^|\t)\\N(?=\t|\r?\n|$)", r"\1", text, flags=re.MULTILINE)
    reader = csv.DictReader(io.StringIO(text), delimiter="\t", escapechar="\\")
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed ICTV export row in {path.name}")
        yield row


def verified_sources(source_dir: Path):
    metadata = {}
    for name, expected in SOURCE_BLOBS.items():
        data = (source_dir / name).read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if blob != expected:
            raise ValueError(f"Source does not match pinned official Git blob: {name}")
        metadata[name] = {"sha256": hashlib.sha256(data).hexdigest(), "git_blob": blob, "bytes": len(data)}
    return metadata


def forward_history_link(link):
    """Do not use reversals or taxon recreations as equivalence mappings."""
    return int(link["rev_count"]) == 0 and int(link["is_recreated"]) == 0


def build_bundle(source_dir: Path, target_release: str):
    sources = verified_sources(source_dir)
    releases = list(read_export(source_dir / "taxonomy_toc.utf8.txt"))
    chosen = [r for r in releases if r["version_tag"] == target_release]
    if len(chosen) != 1:
        raise ValueError("Target release is absent or ambiguous in pinned official release metadata")
    target_number = int(chosen[0]["msl_release_num"])
    raw = list(read_export(source_dir / "taxonomy_node_export.utf8.txt"))
    nodes = {r["taxnode_id"]: r for r in raw if int(r["msl_release_num"]) <= target_number}
    if len(nodes) != sum(int(r["msl_release_num"]) <= target_number for r in raw):
        raise ValueError("Duplicate source taxonomy identifiers")
    visible = {nid: r for nid, r in nodes.items() if r["is_hidden"] != "1" and r["is_deleted"] != "1" and r["is_typo"] != "1"}
    current = {nid: r for nid, r in visible.items() if int(r["msl_release_num"]) == target_number}
    # Include the hidden tree root only as an ancestry container.
    root_id = chosen[0]["tree_id"]
    current[root_id] = nodes[root_id]
    by_identity = defaultdict(set)
    for nid, r in current.items():
        by_identity[r["ictv_id"]].add(nid)
    identity_links = defaultdict(set)
    for identity, targets in by_identity.items():
        identity_links[identity].update(targets)
    for link in read_export(source_dir / "taxonomy_node_merge_split.utf8.txt"):
        # This official table is already a transitive closure. Do not infer
        # edges from similar names or resurrect an abolished taxon.
        if forward_history_link(link):
            identity_links[link["prev_ictv_id"]].update(by_identity.get(link["next_ictv_id"], ()))

    current_ancestors = {}
    for nid in current:
        ancestors, seen, cursor = {}, set(), nid
        while True:
            if cursor in seen or cursor not in current:
                raise ValueError("Invalid target ancestry in official export")
            seen.add(cursor)
            r = current[cursor]
            if r["rank"] in RANKS:
                ancestors[r["rank"]] = cursor
            if r["parent_id"] == cursor:
                break
            cursor = r["parent_id"]
        current_ancestors[nid] = ancestors

    # A realm can be officially renamed while several kingdoms move out of
    # it. An old broad assignment then has multiple possible present-day
    # destinations. Derive this from historical child-taxon relocations, not a
    # hand-maintained Monodnaviria exception.
    membership = defaultdict(set)
    for nid, r in visible.items():
        if int(r["msl_release_num"]) == target_number:
            continue
        destinations = identity_links.get(r["ictv_id"], ())
        cursor, seen = r["parent_id"], set()
        while cursor in nodes and cursor not in seen:
            seen.add(cursor)
            ancestor = nodes[cursor]
            rank = ancestor["rank"]
            if rank in RANKS:
                direct = identity_links.get(ancestor["ictv_id"], ())
                if direct and all(current[t]["rank"] == rank for t in direct):
                    membership[cursor].update(current_ancestors[t][rank] for t in destinations if rank in current_ancestors[t])
                # Only immediate formal child relocations expand a group.
                # Projecting every remote species move onto every ancestor
                # incorrectly makes otherwise unchanged kingdoms ambiguous.
                break
            if ancestor["parent_id"] == cursor:
                break
            cursor = ancestor["parent_id"]

    proposals_by_identity = defaultdict(set)
    for r in nodes.values():
        if r["out_filename"] and int(r["msl_release_num"]) < target_number:
            proposals_by_identity[r["ictv_id"]].add(r["out_filename"])
    aliases = defaultdict(list)
    for nid, r in visible.items():
        if r["rank"] not in RANKS:
            continue
        targets = set(identity_links.get(r["ictv_id"], ())) | membership.get(nid, set())
        aliases[r["rank"] + "\t" + r["name"]].append({
            "source_node_id": nid,
            "ictv_id": r["ictv_id"],
            "release": int(r["msl_release_num"]),
            "targets": sorted(targets),
            "proposals": sorted(proposals_by_identity.get(r["ictv_id"], ())),
        })
    target_nodes = [{"id": nid, "parent": r["parent_id"], "rank": r["rank"], "name": r["name"], "ictv_id": r["ictv_id"]} for nid, r in sorted(current.items())]
    bundle = {
        "schema": SCHEMA,
        "target_release": target_release,
        "target_msl_number": target_number,
        "source_releases": sorted({int(r["msl_release_num"]) for r in visible.values()}),
        "source_commit": SOURCE_COMMIT,
        "source_url": SOURCE_URL,
        "attribution": "Derived from ICTVdatabase, International Committee on Taxonomy of Viruses; CC BY 4.0. Fields projected and historical relationships indexed by viSUM prototype; not endorsed by ICTV.",
        "nodes": target_nodes,
        "aliases": dict(sorted(aliases.items())),
    }
    return bundle, sources


def write_bundle(bundle: dict, sources: dict, output: Path):
    if output.exists() or output.with_name(output.name + ".manifest.json").exists():
        raise ValueError("Refusing to overwrite an existing versioned bundle")
    output.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(bundle, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    output.write_bytes(gzip.compress(data, mtime=0))
    manifest = {
        "schema": SCHEMA,
        "target_release": bundle["target_release"],
        "source_commit": SOURCE_COMMIT,
        "source_url": SOURCE_URL,
        "sources": sources,
        "bundle_sha256": sha256(output),
        "target_node_count": len(bundle["nodes"]),
        "historical_name_rank_count": len(bundle["aliases"]),
        "attribution": bundle["attribution"],
    }
    output.with_name(output.name + ".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-dir", required=True, type=Path)
    p.add_argument("--target-release", default="MSL41.v1")
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    bundle, sources = build_bundle(args.source_dir, args.target_release)
    print(json.dumps(write_bundle(bundle, sources, args.output), indent=2))


if __name__ == "__main__":
    main()
