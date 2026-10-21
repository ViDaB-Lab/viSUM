"""Resolve presentation paths without normalizing scientific output contents.

Historical frozen scripts still write flat files. Current scripts use the new
bundle. Tests must read either layout without changing their value assertions.
"""
from pathlib import Path


def artifact(path):
    path=Path(path)
    if path.exists():return path
    for folder in ('audit','audit/legacy','TE_candidates','quality','groups','review'):
        candidate=path.parent/folder/path.name
        if candidate.exists():return candidate
    # Retain the original nonexistent path for tests asserting absent outputs.
    return path


def artifacts(prefix):
    return {p.name:p for p in prefix.parent.rglob(prefix.name+'.*') if p.is_file()}
