#!/usr/bin/env python3
"""Compare pipeline output tables against an immutable reference baseline.

The baseline is a directory of previously archived output tables. This script
never writes to it. For every table it reports file existence, SHA-256 digest,
row and column counts, column names and order, and cell-level differences, and
writes a machine-readable manifest.

Differences are not interpreted here. A difference is reported as unexplained
unless it is listed in an accepted-differences file, which records why each one
exists. An unexplained difference makes the script exit non-zero.

Accepted differences are matched one cell at a time. Each entry names the
table, the row key, the column, the exact baseline value and the exact
current value, and must state the decision, the reason and the evidence
status. A cell is accepted only if all five identifying fields match, so an
entry cannot excuse a change in another row, another column, or a different
change in the same cell. An entry that matches no difference is reported and
fails the check too, because an expectation that no longer holds means the
record is out of date. A table whose shape or column names changed fails
outright.

The accepted-differences file can also declare output tables that were renamed
since the baseline was archived, so that content is still compared even though
the file name changed. It describes differences against a private baseline and
quotes values from it, so it is a private file; a tracked template shows its
structure.

Usage
-----
    python scripts/check_regression.py --baseline DIR --current DIR \\
        [--accepted FILE] [--manifest FILE]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd


def digest(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p: Path) -> pd.DataFrame:
    sep = "\t" if p.suffix == ".tsv" else ";"
    return pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False, encoding="utf-8-sig")


#: Fields every accepted-difference entry must have.
REQUIRED_FIELDS = ("table", "key", "column", "baseline", "current",
                   "decision", "reason", "evidence_status")

#: Allowed evidence_status values. None of them claims independent review.
EVIDENCE_STATUS = ("provisional", "pending_second_review", "verified_against_primary_source")


class SpecError(ValueError):
    """The accepted-differences file is malformed."""


def load_spec(path: Path | None) -> tuple[list[dict], list[dict], dict[str, str]]:
    if path is None:
        return [], [], {}
    if not path.exists():
        raise SpecError(f"accepted-differences file not found: {path}")
    spec = json.loads(path.read_text())
    accepted = spec.get("accepted_differences", [])
    seen = set()
    for i, d in enumerate(accepted):
        missing = [f for f in REQUIRED_FIELDS if not str(d.get(f, "")).strip()
                   and f not in ("baseline", "current")]
        missing += [f for f in ("baseline", "current") if f not in d]
        if missing:
            raise SpecError(f"accepted difference #{i} lacks {missing}")
        if d["evidence_status"] not in EVIDENCE_STATUS:
            raise SpecError(f"accepted difference #{i}: evidence_status must be one of "
                            f"{list(EVIDENCE_STATUS)}")
        cell = (d["table"], d["key"], d["column"])
        if cell in seen:
            raise SpecError(f"accepted difference #{i} repeats the cell {cell}")
        seen.add(cell)
    obsolete = spec.get("obsolete_baseline_tables", [])
    for i, d in enumerate(obsolete):
        if not d.get("table") or not d.get("reason"):
            raise SpecError(f"obsolete table #{i} needs 'table' and 'reason'")
    renamed = {}
    for i, d in enumerate(spec.get("renamed_tables", [])):
        if not d.get("baseline") or not d.get("current") or not d.get("reason"):
            raise SpecError(f"renamed table #{i} needs 'baseline', 'current' and 'reason'")
        renamed[d["baseline"]] = d["current"]
    return accepted, obsolete, renamed


#: Columns tried, in order, as the row identifier. The first column whose values
#: are unique in both tables is used. A non-unique key would let one accepted
#: difference excuse changes in several rows, so it is never used.
KEY_PREFERENCE = ("variant", "patient_id", "stratum", "definition")


def row_key(a: pd.DataFrame, b: pd.DataFrame) -> str | None:
    for col in (*KEY_PREFERENCE, a.columns[0]):
        if col in a.columns and col in b.columns and a[col].is_unique and b[col].is_unique:
            return col
    return None


def compare(baseline: Path, current: Path) -> dict:
    a, b = read(current), read(baseline)
    result: dict = {
        "sha256_baseline": digest(baseline),
        "sha256_current": digest(current),
        "rows": [len(b), len(a)],
        "columns": [len(b.columns), len(a.columns)],
        "column_names_match": list(a.columns) == list(b.columns),
        "differences": [],
    }
    result["identical_bytes"] = result["sha256_baseline"] == result["sha256_current"]
    if not result["column_names_match"] or a.shape != b.shape:
        result["status"] = "structure_differs"
        return result
    key_col = row_key(a, b)
    result["key_column"] = key_col
    mask = a.values != b.values
    for i, j in zip(*mask.nonzero()):
        key = str(a[key_col].iloc[i]) if key_col else f"row:{int(i)}"
        result["differences"].append({
            "row": int(i),
            "key": key,
            "column": str(a.columns[j]),
            "baseline": str(b.iloc[i, j]),
            "current": str(a.iloc[i, j]),
        })
    result["status"] = "identical" if not result["differences"] else "differs"
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", required=True, type=Path)
    ap.add_argument("--current", required=True, type=Path)
    ap.add_argument("--accepted", type=Path, default=None,
                    help="JSON list of accepted differences")
    ap.add_argument("--manifest", type=Path, default=None)
    args = ap.parse_args()

    try:
        accepted, obsolete, renamed = load_spec(args.accepted)
    except SpecError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    expected = {(d["table"], d["key"], d["column"]): d for d in accepted}
    used: set[tuple[str, str, str]] = set()
    obsolete_tables = {d["table"] for d in obsolete}

    report: dict = {"tables": {}, "summary": {}}
    unexplained = 0
    for base_file in sorted(args.baseline.glob("*.csv")):
        cur = args.current / renamed.get(base_file.name, base_file.name)
        if not cur.exists():
            if base_file.name in obsolete_tables:
                report["tables"][base_file.name] = {"status": "obsolete_baseline_table"}
            else:
                report["tables"][base_file.name] = {"status": "missing_in_current"}
                unexplained += 1
            continue
        r = compare(base_file, cur)
        if cur.name != base_file.name:
            r["compared_as"] = cur.name
        if r["status"] == "structure_differs":
            unexplained += 1
        for d in r["differences"]:
            cell = (base_file.name, d["key"], d["column"])
            e = expected.get(cell)
            d["accepted"] = (e is not None and e["baseline"] == d["baseline"]
                             and e["current"] == d["current"])
            if d["accepted"]:
                used.add(cell)
                d["evidence_status"] = e["evidence_status"]
            else:
                unexplained += 1
                if e is not None:
                    d["note"] = (f"an entry exists for this cell, but expects "
                                 f"{e['baseline']!r} -> {e['current']!r}")
        report["tables"][base_file.name] = r

    unused = sorted(set(expected) - used)
    unexplained += len(unused)

    renamed_current = set(renamed.values())
    extra = sorted(
        p.name for p in args.current.glob("*.csv")
        if not (args.baseline / p.name).exists()
        and p.name not in renamed_current and p.name != "run_manifest.json"
    )
    identical = [k for k, v in report["tables"].items() if v.get("status") == "identical"]
    obsolete_seen = [k for k, v in report["tables"].items()
                     if v.get("status") == "obsolete_baseline_table"]
    differing = [k for k, v in report["tables"].items() if v.get("status") == "differs"]
    report["summary"] = {
        "tables_compared": len(report["tables"]),
        "identical": sorted(identical),
        "differing": sorted(differing),
        "not_in_baseline": extra,
        "obsolete_baseline_tables": sorted(obsolete_seen),
        "renamed_tables": {b: c for b, c in sorted(renamed.items())},
        "unused_accepted_entries": [list(c) for c in unused],
        "unexplained_differences": unexplained,
    }

    target = args.manifest or (args.current / "regression_manifest.json")
    target.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    s = report["summary"]
    print(f"tables compared     : {s['tables_compared']}")
    print(f"identical           : {len(s['identical'])}")
    print(f"differing           : {len(s['differing'])} {s['differing'] or ''}")
    print(f"new since baseline  : {s['not_in_baseline'] or 'none'}")
    print(f"obsolete in baseline: {s['obsolete_baseline_tables'] or 'none'}")
    for b, c in s["renamed_tables"].items():
        print(f"renamed             : {b} -> {c}")
    structure = [k for k, v in report["tables"].items()
                 if v.get("status") in ("structure_differs", "missing_in_current")]
    for name in structure:
        print(f"  [UNEXPLAINED] {name}: {report['tables'][name]['status']}")
    for c in s["unused_accepted_entries"]:
        print(f"  [UNUSED ENTRY] {c[0]} :: {c[1]} / {c[2]} matched no difference")
    print(f"unexplained         : {s['unexplained_differences']}")
    for name in s["differing"]:
        for d in report["tables"][name]["differences"]:
            tag = (f"accepted, {d['evidence_status']}" if d["accepted"] else "UNEXPLAINED")
            print(f"  [{tag}] {name} :: {d['key']} / {d['column']}: "
                  f"{d['baseline']!r} -> {d['current']!r}" + (f"  ({d['note']})" if d.get("note") else ""))
    print(f"manifest -> {target}")
    return 1 if unexplained else 0


if __name__ == "__main__":
    sys.exit(main())
