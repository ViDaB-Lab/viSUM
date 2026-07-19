"""Shared column definitions for sparse viSUM tool-evidence tables."""

TAXONOMY_COLUMNS = [
    "d__Domain",
    "r__Realm",
    "k__Kingdom",
    "p__Phylum",
    "c__Class",
    "o__Order",
    "f__Family",
    "g__Genus",
    "s__Species",
]

CORE_EVIDENCE_COLUMNS = [
    "sample_id",
    "sequence_id",
    "parent_sequence_id",
    "record_type",
    "coordinates",
    "tool",
    "classification",
    "score",
    "score_type",
    "length",
    "topology",
    "n_genes",
    "n_hallmarks",
    *TAXONOMY_COLUMNS,
]

RANK_PREFIXES = ("d__", "r__", "k__", "p__", "c__", "o__", "f__", "g__", "s__")


def unclassified_taxonomy(classification: str) -> dict[str, str]:
    """Return rank columns with a viral domain only for viral evidence."""
    ranks = ["unclassified"] * len(TAXONOMY_COLUMNS)
    if classification == "virus":
        ranks[0] = "Viruses"

    return {
        column: f"{prefix}{rank}"
        for column, prefix, rank in zip(TAXONOMY_COLUMNS, RANK_PREFIXES, ranks)
    }
