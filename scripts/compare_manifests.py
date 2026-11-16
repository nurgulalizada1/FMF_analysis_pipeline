#!/usr/bin/env python3
"""Compare the deterministic sections of two run manifests.

Exits 0 if the scientific outputs, intermediates, settings and stage results
recorded by two runs are identical, and 1 otherwise, listing what differs.
The ``execution`` section (timestamps, durations, paths, versions) is ignored
by design: it is not expected to be reproducible.

Usage
-----
    python scripts/compare_manifests.py FIRST.json SECOND.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    a, b = (json.loads(Path(p).read_text())["deterministic"] for p in argv)
    problems = []
    for section in ("artifacts", "intermediates"):
        ka, kb = set(a[section]), set(b[section])
        problems += [f"{section}: only in first: {k}" for k in sorted(ka - kb)]
        problems += [f"{section}: only in second: {k}" for k in sorted(kb - ka)]
        problems += [f"{section}: content differs: {k}" for k in sorted(ka & kb)
                     if a[section][k]["sha256"] != b[section][k]["sha256"]]
    for key in ("settings", "stages", "artifacts_sha256"):
        if a[key] != b[key]:
            problems.append(f"{key} differ")
    if problems:
        print("\n".join(problems))
        return 1
    print(f"identical: {len(a['artifacts'])} artifacts, {len(a['intermediates'])} intermediates, "
          f"combined digest {a['artifacts_sha256'][:16]}…")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
