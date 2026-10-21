#!/usr/bin/env python3
"""Standalone, offline ICTV name reconciliation. Not used by the pipeline.

This translates nomenclature, not evidence strength or biological identity.
The caller must preserve its original prediction and make its own decisions.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
from pathlib import Path


RANKS = ("realm", "kingdom", "phylum", "class", "order", "family", "genus", "species")
SCHEMA = "visum-ictv-reconciliation-prototype-v1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@dataclass
class Resolution:
    original: dict[str, str]
    normalized: dict[str, str]
    status: str
    target_release: str
    source_release: int | None
    unresolved: list[dict]
    mappings: list[dict]


class Reconciler:
    def __init__(self, bundle: dict):
        if bundle.get("schema") != SCHEMA:
            raise ValueError("Unsupported reconciliation bundle schema")
        self.bundle = bundle
        self.target_release = bundle["target_release"]
        self.nodes = {str(n["id"]): n for n in bundle["nodes"]}
        if not self.nodes:
            raise ValueError("Empty target taxonomy")
        if len(self.nodes) != len(bundle["nodes"]):
            raise ValueError("Duplicate target node identifiers")
        self.aliases = bundle["aliases"]
        self.lineages: dict[str, dict[str, str]] = {}
        self.children: dict[str, list[str]] = {nid: [] for nid in self.nodes}
        self.current: dict[str, list[str]] = {}
        for nid, node in self.nodes.items():
            parent = str(node["parent"])
            if parent != nid:
                if parent not in self.nodes:
                    raise ValueError(f"Missing target parent: {parent}")
                self.children[parent].append(nid)
            if node["rank"] in RANKS:
                self.current.setdefault(node["rank"] + "\t" + node["name"], []).append(nid)
        # Iterative ancestry validation avoids recursion limits on corrupt data.
        for nid in self.nodes:
            trail, seen, cursor = [], set(), nid
            while cursor not in self.lineages:
                if cursor in seen:
                    raise ValueError("Cycle in target taxonomy")
                seen.add(cursor)
                node = self.nodes[cursor]
                trail.append(cursor)
                if str(node["parent"]) == cursor:
                    self.lineages[cursor] = {}
                    break
                cursor = str(node["parent"])
            for item in reversed(trail):
                node = self.nodes[item]
                lineage = dict(self.lineages.get(str(node["parent"]), {}))
                if node["rank"] in RANKS:
                    if node["rank"] in lineage:
                        raise ValueError("Repeated formal rank in target ancestry")
                    lineage[node["rank"]] = node["name"]
                self.lineages[item] = lineage
        self.bits = {nid: 1 << i for i, nid in enumerate(sorted(self.nodes))}
        self.subtrees = dict(self.bits)
        # The ancestry length includes subranks, unlike formal lineage length.
        depths = {}
        for nid in self.nodes:
            cursor, depth = nid, 0
            while str(self.nodes[cursor]["parent"]) != cursor:
                depth += 1
                cursor = str(self.nodes[cursor]["parent"])
            depths[nid] = depth
        for nid in sorted(self.nodes, key=depths.get, reverse=True):
            parent = str(self.nodes[nid]["parent"])
            if parent != nid:
                self.subtrees[parent] |= self.subtrees[nid]
        self.all_bits = (1 << len(self.nodes)) - 1
        self.by_bit = {value: nid for nid, value in self.bits.items()}
        for records in self.aliases.values():
            for record in records:
                if any(str(nid) not in self.nodes for nid in record["targets"]):
                    raise ValueError("Historical mapping points outside target taxonomy")

    @classmethod
    def load(cls, path: Path) -> "Reconciler":
        manifest_path = path.with_name(path.name + ".manifest.json")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("bundle_sha256") != sha256(path):
            raise ValueError("Reconciliation bundle checksum mismatch")
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            bundle = json.load(handle)
        if bundle["target_release"] != manifest.get("target_release"):
            raise ValueError("Bundle/manifest target release mismatch")
        return cls(bundle)

    def _options(self, rank: str, name: str, source_release: int | None):
        key = rank + "\t" + name
        records = self.aliases.get(key, [])
        if source_release is not None:
            records = [r for r in records if r["release"] == source_release]
            return sorted({str(n) for r in records for n in r["targets"]}), records
        # Without a known source release, an exact current name retains its
        # current meaning. Obsolete names require historical relationships.
        if key in self.current:
            return self.current[key], []
        return sorted({str(n) for r in records for n in r["targets"]}), records

    def resolve(self, lineage: dict[str, str], source_release: int | None = None) -> Resolution:
        strict = self._resolve_strict(lineage, source_release)
        if strict.status != "incompatible":
            return strict

        # A named higher prediction is unsupported, not a branch conflict,
        # when *every* surviving lower placement has that rank unassigned.
        # Work upward from jointly compatible lower claims; never select one
        # side of an actual assigned-ancestor conflict or an ambiguous split.
        known = []
        for rank, name in strict.original.items():
            options, records = self._options(rank, name, source_release)
            if not options:
                break  # Preserve the existing unknown-rank barrier.
            known.append((rank, name, options, records))
        # Check *all* known claim pairs before omitting any labels. Otherwise
        # a lower unassigned taxon could hide an explicit conflict between two
        # upper predictions (e.g. a realm and an order in another realm).
        for index, (rank, _, options, _) in enumerate(known):
            if any(self.nodes[n]["rank"] != rank for n in options):
                continue
            allowed = 0
            for nid in options:
                allowed |= self.subtrees[nid]
            for _, _, lower_options, _ in known[index + 1:]:
                if (all(rank in self.lineages[n] for n in lower_options)
                        and not any(allowed & self.subtrees[n] for n in lower_options)):
                    return strict
        filtered = dict(strict.original)
        omitted, omitted_mappings = [], []
        for index in range(len(known) - 2, -1, -1):
            rank, name, options, records = known[index]
            suffix = {r: n for r, n, _, _ in known[index + 1:] if r in filtered}
            lower = self._resolve_strict(suffix, source_release)
            if lower.status == "incompatible":
                return strict  # A genuine conflict survives below this rank.
            if not lower.mappings:
                continue
            viable = lower.mappings[-1]["compatible_target_node_ids"]
            rank_index = RANKS.index(rank)
            if (not viable or any(self.nodes[n]["rank"] != rank for n in options)
                    or not all(rank not in self.lineages[n]
                               and max((RANKS.index(r) for r in self.lineages[n]), default=-1) > rank_index
                               for n in viable)):
                continue
            del filtered[rank]
            omitted.append({"rank": rank, "name": name,
                            "reason": "target_higher_rank_unassigned",
                            "supporting_lower_target_node_ids": viable})
            omitted_mappings.append({"rank": rank, "name": name,
                "target_node_ids": options, "history": records,
                "compatible_target_node_ids": [],
                "disposition": "not_used_target_higher_rank_unassigned"})
        if not omitted:
            return strict
        corrected = self._resolve_strict(filtered, source_release)
        if corrected.status == "incompatible":
            return strict
        return Resolution(strict.original, corrected.normalized, "partial",
            self.target_release, source_release,
            sorted(omitted + corrected.unresolved, key=lambda item: RANKS.index(item["rank"])),
            sorted(omitted_mappings + corrected.mappings, key=lambda item: RANKS.index(item["rank"])))

    def _resolve_strict(self, lineage: dict[str, str], source_release: int | None = None) -> Resolution:
        unknown_ranks = set(lineage) - set(RANKS)
        if unknown_ranks:
            raise ValueError(f"Unsupported ranks: {sorted(unknown_ranks)}")
        if source_release is not None and source_release not in self.bundle["source_releases"]:
            raise ValueError("Source release not present in this pinned history snapshot")
        original = {r: str(lineage.get(r, "")).strip() for r in RANKS}
        original = {r: n for r, n in original.items() if n and n.lower() not in {"unclassified", "na", "none", "-"}}
        if not original:
            return Resolution(original, {}, "no_named_taxonomy", self.target_release, source_release, [], [])
        candidates, mappings, unresolved, deepest = self.all_bits, [], [], -1
        for rank in RANKS:
            if rank not in original:
                continue
            name = original[rank]
            options, records = self._options(rank, name, source_release)
            if not options:
                unresolved.append({"rank": rank, "name": name, "reason": "no_verified_current_mapping"})
                # Do not use deeper claims to guess through an unknown rank.
                break
            allowed = 0
            for nid in options:
                allowed |= self.subtrees[nid]
            intersection = candidates & allowed
            if not intersection:
                unresolved.append({"rank": rank, "name": name, "reason": "incompatible_lineage"})
                return Resolution(original, {}, "incompatible", self.target_release, source_release, unresolved, mappings)
            candidates = intersection
            # If a formal rank was demoted to an intermediate rank that viSUM
            # does not report (e.g. subfamily), stop at its closest supported
            # formal ancestor rather than inferring one of its descendants.
            mapped_rank_indices = [max((RANKS.index(r) for r in self.lineages[n]), default=-1) for n in options]
            deepest = max(deepest, max(mapped_rank_indices, default=-1))
            mappings.append({"rank": rank, "name": name, "target_node_ids": options, "history": records})
        surviving_bits = candidates
        normalized = {}
        representative = self.by_bit[candidates & -candidates]
        for rank in RANKS[:deepest + 1]:
            name = self.lineages[representative].get(rank, "")
            if not name:
                continue
            containing = 0
            for nid in self.current[rank + "\t" + name]:
                containing |= self.subtrees[nid]
            if not (surviving_bits & ~containing):
                normalized[rank] = name
        ambiguous = False
        for mapping in mappings:
            viable = [nid for nid in mapping["target_node_ids"] if surviving_bits & self.subtrees[nid]]
            mapping["compatible_target_node_ids"] = viable
            # Parent and descendant alternatives are still biologically
            # ambiguous; no member is selected merely because it exists.
            if len(viable) > 1:
                ambiguous = True
        names_changed = any(normalized.get(rank) != name for rank, name in original.items())
        status = "partial" if unresolved else "ambiguous" if ambiguous else "translated" if names_changed else "unchanged"
        return Resolution(original, normalized, status, self.target_release, source_release, unresolved, mappings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--lineage-json", required=True, help="JSON object of unprefixed rank/name pairs")
    parser.add_argument("--source-release", type=int)
    args = parser.parse_args()
    result = Reconciler.load(args.bundle).resolve(json.loads(args.lineage_json), args.source_release)
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
