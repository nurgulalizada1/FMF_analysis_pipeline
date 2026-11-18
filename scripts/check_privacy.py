#!/usr/bin/env python3
"""Scan the repository, staged changes or generated outputs for data that should not be published.

    python scripts/check_privacy.py                  # tracked + untracked-not-ignored files
    python scripts/check_privacy.py --tracked        # tracked files only
    python scripts/check_privacy.py --staged         # what the next commit would add (pre-commit)
    python scripts/check_privacy.py --path results/  # generated outputs, notebooks, any folder
    python scripts/check_privacy.py --denylist ~/private/denylist.txt

The private denylist is read from --denylist or the MEFV_PRIVATE_DENYLIST environment
variable. It must live outside the repository. Without one only the built-in rules run, and
the report says so.

Output: file paths, rule identifiers and counts. Matched values are never printed.
Exit status: 0 no findings, 1 findings, 2 the scan could not run or the denylist is unusable.
A clean result is not a guarantee that no sensitive data are present (docs/privacy_policy.md).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mefv import privacy  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--tracked", action="store_true", help="tracked files only")
    mode.add_argument("--staged", action="store_true", help="staged additions and changes")
    mode.add_argument("--path", nargs="+", metavar="P", help="scan these files or folders")
    ap.add_argument("--denylist", help=f"private denylist file (default: ${privacy.DENYLIST_ENV})")
    ap.add_argument("--require-denylist", action="store_true",
                    help="exit 2 if no private denylist is available")
    ap.add_argument("--json", metavar="FILE", help="also write the findings as JSON (no values)")
    ap.add_argument("--root", default=".", help="repository root or a folder inside it")
    a = ap.parse_args(argv)

    try:
        root = privacy.repo_root(Path(a.root))
    except privacy.ScanError as e:
        print(f"check_privacy: {e}", file=sys.stderr)
        return 2

    dl_path = a.denylist or os.environ.get(privacy.DENYLIST_ENV)
    try:
        denylist = privacy.load_denylist(dl_path, root)
    except privacy.DenylistError as e:
        print(f"check_privacy: {e}", file=sys.stderr)
        return 2
    if privacy.tracked_denylist_files(root, dl_path):
        print("check_privacy: the private denylist is tracked by Git; untrack it", file=sys.stderr)
        return 2
    if denylist is None and dl_path:
        # A denylist that was asked for but cannot be found must not turn into a scan that
        # quietly runs without it.
        print("check_privacy: the configured private denylist was not found", file=sys.stderr)
        return 2
    if denylist is None and a.require_denylist:
        print("check_privacy: no private denylist available (--require-denylist)", file=sys.stderr)
        return 2

    try:
        if a.staged:
            items = privacy.staged_items(root)
        elif a.tracked:
            items = privacy.tracked_items(root)
        elif a.path:
            items = privacy.path_items(a.path, root)
        else:
            items = privacy.working_tree_items(root)
    except privacy.ScanError as e:
        print(f"check_privacy: scan could not be completed: {e}", file=sys.stderr)
        return 2

    result = privacy.scan_items(items, denylist, path_rules=not a.path)  # --path: content rules only
    print(result.render())
    if a.json:
        Path(a.json).write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n")
    return 0 if result.clean else 1


if __name__ == "__main__":
    raise SystemExit(main())
