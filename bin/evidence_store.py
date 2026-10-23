"""Compact, ordered TSV storage for parent-side viHARMONY inputs.

Only representation changes here. Candidate evaluation and reporting receive
ordinary dictionaries. TSV parsing preserves csv.DictReader's field ordering,
duplicate-header, missing-field and extra-field behavior, including None keys.
"""
from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
import csv
from pathlib import Path
import sys


# Share bounded categorical values, not arbitrary identifiers or free text.
_CATEGORICAL = frozenset((
    "sample_id", "record_type", "tool", "classification", "evidence_strength",
    "strength_basis", "score_type", "topology", "assignment_method",
    "vcontact3_assignment_method", "vcontact3_database_domain", "vcontact3_database_version",
    "classification_rank", "tesorter_order", "tesorter_superfamily", "tesorter_evidence_category",
    "d__Domain", "r__Realm", "k__Kingdom", "p__Phylum", "c__Class", "o__Order",
    "f__Family", "g__Genus", "s__Species",
))
_ROW_BITS = 32
_ROW_LIMIT = 1 << _ROW_BITS
_ROW_MASK = _ROW_LIMIT - 1


class _Schema:
    __slots__ = ("names", "positions")

    def __init__(self, names):
        self.names = tuple(names)
        self.positions = {name: position for position, name in enumerate(self.names)}


class EvidenceTable:
    """One shared schema plus tuple values, with lightweight read-only row views."""

    def __init__(self, path: Path, *, source_file: str | None = None):
        self.path = Path(path)
        self._rows: list[tuple] = []
        with self.path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            fields = next(reader, None)
            if fields is None:
                raise ValueError(f"TSV has no header: {path}")
            self.fieldnames = list(fields)
            names = tuple(dict.fromkeys(fields))
            append_source = source_file is not None and "_source_file" not in names
            source_names = ("_source_file",) if append_source else ()
            self._normal = _Schema((*names, *source_names))
            self._extra = _Schema((*names, None, *source_names))
            self._normal_width = len(self._normal.names)
            width = len(fields)
            pooled = [(index, key) for index, key in enumerate(names) if key in _CATEGORICAL]

            if len(names) != len(fields):
                # Duplicate headers require DictReader's last-value/restval rules.
                # The file is already positioned after its header; do not reread it.
                duplicate_reader = csv.DictReader(handle, fieldnames=fields, delimiter="\t")
                for row in duplicate_reader:
                    values = [row[key] for key in names]
                    if None in row:
                        values.append(row[None])
                    self._append(values, pooled, names, source_file, append_source)
            else:
                for values in reader:
                    # DictReader skips empty CSV records, not tab-only records.
                    if not values:
                        continue
                    if len(values) < width:
                        values.extend([None] * (width - len(values)))
                    elif len(values) > width:
                        values = values[:width] + [values[width:]]
                    self._append(values, pooled, names, source_file, append_source)

    def _append(self, values, pooled, names, source_file, append_source):
        if len(self._rows) >= _ROW_LIMIT:
            raise OverflowError("Compact evidence supports fewer than 2**32 rows per table")
        for index, _key in pooled:
            value = values[index]
            if isinstance(value, str) and value:
                values[index] = sys.intern(value)
        if source_file is not None:
            if append_source:
                values.append(source_file)
            else:
                values[self._normal.positions["_source_file"]] = source_file
        self._rows.append(tuple(values))

    def __len__(self):
        return len(self._rows)

    def __iter__(self) -> Iterator["EvidenceRow"]:
        for position in range(len(self._rows)):
            yield EvidenceRow(self, position)

    def row(self, position: int) -> "EvidenceRow":
        if position < 0 or position >= len(self._rows):
            raise IndexError(position)
        return EvidenceRow(self, position)

    def _schema(self, values):
        return self._extra if len(values) > self._normal_width else self._normal

    def materialize(self, position: int) -> dict:
        values = self._rows[position]
        return dict(zip(self._schema(values).names, values))


class EvidenceRow(Mapping):
    """A row reference; never sent to workers or to biological decision code."""

    __slots__ = ("table", "position")

    def __init__(self, table: EvidenceTable, position: int):
        self.table, self.position = table, position

    def __getitem__(self, key):
        values = self.table._rows[self.position]
        return values[self.table._schema(values).positions[key]]

    def get(self, key, default=None):
        values = self.table._rows[self.position]
        position = self.table._schema(values).positions.get(key)
        return default if position is None else values[position]

    def __iter__(self):
        return iter(self.table._schema(self.table._rows[self.position]).names)

    def __len__(self):
        return len(self.table._rows[self.position])


class EvidenceIndex(Mapping):
    """Ordered row references; each lookup returns a plain list of plain dicts."""

    def __init__(self):
        self._tables: list[EvidenceTable] = []
        self._table_ids: dict[int, int] = {}
        self._references: dict[str, list[int]] = {}

    def add(self, key: str, row: EvidenceRow) -> None:
        if not isinstance(row, EvidenceRow):
            raise TypeError("EvidenceIndex requires EvidenceRow references")
        if not 0 <= row.position < _ROW_LIMIT:
            raise OverflowError("Compact evidence row reference exceeds its supported range")
        identity = id(row.table)
        table_id = self._table_ids.get(identity)
        if table_id is None:
            table_id = len(self._tables)
            self._table_ids[identity] = table_id
            self._tables.append(row.table)
        token = (table_id << _ROW_BITS) | row.position
        self._references.setdefault(key, []).append(token)

    @classmethod
    def from_pairs(cls, pairs: Iterable[tuple[str, EvidenceRow]]) -> "EvidenceIndex":
        result = cls()
        for key, row in pairs:
            result.add(key, row)
        return result

    def __getitem__(self, key):
        references = self._references[key]
        # Preserve row multiplicity and within-lookup identity for repeated refs.
        materialized = {}
        result = []
        for token in references:
            if token not in materialized:
                materialized[token] = self._tables[token >> _ROW_BITS].materialize(token & _ROW_MASK)
            result.append(materialized[token])
        return result

    def get(self, key, default=None):
        return self[key] if key in self._references else default

    def __iter__(self):
        return iter(self._references)

    def __len__(self):
        return len(self._references)
