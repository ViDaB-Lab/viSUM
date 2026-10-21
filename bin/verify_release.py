#!/usr/bin/env python3
"""Verify this unpacked source package; no installation, network or analysis."""
import json
from pathlib import Path
from verify_alpha_candidate import verify

if __name__=='__main__':
    print(json.dumps(verify(Path(__file__).resolve().parents[1]),indent=2))
