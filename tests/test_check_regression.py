"""Tests for the regression comparison tool.

The tool decides whether a change is explained. A defect here would let an
unexplained change pass silently, so its decision logic is tested directly.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "scripts" / "check_regression.py"
_spec = importlib.util.spec_from_file_location("check_regression", TOOL)
cr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cr)


def write(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep=";", index=False, encoding="utf-8-sig")


def table(values: dict[str, str]) -> pd.DataFrame:
    # The first column repeats on every row, as in a real stratified table.
    return pd.DataFrame({
        "cohort": ["all"] * len(values),
        "variant": list(values),
        "flag": list(values.values()),
    })


def entry(key: str, column: str = "flag", baseline: str = "x", current: str = "CHANGED",
          **extra) -> dict:
    return {"table": "t.csv", "key": key, "column": column, "baseline": baseline,
            "current": current, "decision": "test decision", "reason": "test reason",
            "evidence_status": "provisional", **extra}


def spec_file(tmp_path: Path, **sections) -> Path:
    p = tmp_path / "accepted.json"
    p.write_text(json.dumps(sections))
    return p


def run(base: Path, cur: Path, accepted: Path | None = None) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(TOOL), "--baseline", str(base), "--current", str(cur),
           "--manifest", str(cur / "m.json")]
    if accepted:
        cmd += ["--accepted", str(accepted)]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_row_key_skips_a_non_unique_first_column():
    a = table({"V1": "x", "V2": "y"})
    assert cr.row_key(a, a) == "variant"


def test_row_key_is_none_when_nothing_is_unique():
    a = pd.DataFrame({"g": ["a", "a"], "h": ["b", "b"]})
    assert cr.row_key(a, a) is None


def test_identical_tables_pass(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x", "V2": "y"}))
    write(tmp_path / "cur/t.csv", table({"V1": "x", "V2": "y"}))
    assert run(tmp_path / "base", tmp_path / "cur").returncode == 0


def test_an_unexplained_difference_fails(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x", "V2": "y"}))
    write(tmp_path / "cur/t.csv", table({"V1": "x", "V2": "CHANGED"}))
    r = run(tmp_path / "base", tmp_path / "cur")
    assert r.returncode != 0
    assert "UNEXPLAINED" in r.stdout


def test_an_accepted_difference_passes_only_for_its_own_row(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x", "V2": "y"}))
    write(tmp_path / "cur/t.csv", table({"V1": "CHANGED", "V2": "CHANGED"}))
    spec = spec_file(tmp_path, accepted_differences=[entry("V1")])
    r = run(tmp_path / "base", tmp_path / "cur", spec)
    # V1 is accepted, V2 is not: one shared 'cohort' value must not excuse both.
    assert r.returncode != 0
    assert "[accepted, provisional] t.csv :: V1" in r.stdout
    assert "[UNEXPLAINED] t.csv :: V2" in r.stdout


def test_an_exactly_matching_entry_passes(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x", "V2": "y"}))
    write(tmp_path / "cur/t.csv", table({"V1": "CHANGED", "V2": "y"}))
    spec = spec_file(tmp_path, accepted_differences=[entry("V1")])
    assert run(tmp_path / "base", tmp_path / "cur", spec).returncode == 0


def test_a_different_change_in_an_accepted_cell_fails(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x"}))
    write(tmp_path / "cur/t.csv", table({"V1": "SOMETHING ELSE"}))
    spec = spec_file(tmp_path, accepted_differences=[entry("V1")])
    r = run(tmp_path / "base", tmp_path / "cur", spec)
    assert r.returncode != 0 and "expects" in r.stdout


def test_an_entry_for_another_column_does_not_excuse_a_change(tmp_path):
    base = pd.DataFrame({"variant": ["V1"], "a": ["1"], "b": ["2"]})
    cur = pd.DataFrame({"variant": ["V1"], "a": ["1"], "b": ["9"]})
    write(tmp_path / "base/t.csv", base)
    write(tmp_path / "cur/t.csv", cur)
    spec = spec_file(tmp_path, accepted_differences=[entry("V1", column="a",
                                                           baseline="2", current="9")])
    assert run(tmp_path / "base", tmp_path / "cur", spec).returncode != 0


def test_an_unused_entry_fails(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x"}))
    write(tmp_path / "cur/t.csv", table({"V1": "x"}))
    spec = spec_file(tmp_path, accepted_differences=[entry("V1")])
    r = run(tmp_path / "base", tmp_path / "cur", spec)
    assert r.returncode != 0 and "UNUSED ENTRY" in r.stdout


def test_a_structural_change_fails(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x"}))
    write(tmp_path / "cur/t.csv", table({"V1": "x", "V2": "y"}))
    r = run(tmp_path / "base", tmp_path / "cur")
    assert r.returncode != 0 and "structure_differs" in r.stdout


def test_an_incomplete_entry_is_rejected(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x"}))
    write(tmp_path / "cur/t.csv", table({"V1": "CHANGED"}))
    bad = entry("V1")
    del bad["evidence_status"]
    r = run(tmp_path / "base", tmp_path / "cur", spec_file(tmp_path, accepted_differences=[bad]))
    assert r.returncode == 2 and "evidence_status" in r.stderr


def test_an_unknown_evidence_status_is_rejected(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x"}))
    write(tmp_path / "cur/t.csv", table({"V1": "CHANGED"}))
    spec = spec_file(tmp_path, accepted_differences=[entry("V1", evidence_status="approved")])
    assert run(tmp_path / "base", tmp_path / "cur", spec).returncode == 2


def test_a_missing_accepted_file_is_an_error_not_an_empty_list(tmp_path):
    write(tmp_path / "base/t.csv", table({"V1": "x"}))
    write(tmp_path / "cur/t.csv", table({"V1": "x"}))
    r = run(tmp_path / "base", tmp_path / "cur", tmp_path / "absent.json")
    assert r.returncode == 2


def test_a_missing_table_fails_unless_declared_obsolete(tmp_path):
    write(tmp_path / "base/old.csv", table({"V1": "x"}))
    (tmp_path / "cur").mkdir()
    assert run(tmp_path / "base", tmp_path / "cur").returncode != 0
    spec = spec_file(tmp_path, obsolete_baseline_tables=[{"table": "old.csv", "reason": "r"}])
    assert run(tmp_path / "base", tmp_path / "cur", spec).returncode == 0


def test_a_renamed_table_is_still_compared(tmp_path):
    write(tmp_path / "base/old_name.csv", table({"V1": "x"}))
    write(tmp_path / "cur/new_name.csv", table({"V1": "DIFFERENT"}))
    spec = spec_file(tmp_path, renamed_tables=[
        {"baseline": "old_name.csv", "current": "new_name.csv", "reason": "r"}])
    r = run(tmp_path / "base", tmp_path / "cur", spec)
    # the rename is honoured, so the content difference is found, not skipped
    assert r.returncode != 0
    assert "UNEXPLAINED" in r.stdout


def test_the_baseline_is_never_written(tmp_path):
    base = tmp_path / "base"
    write(base / "t.csv", table({"V1": "x"}))
    write(tmp_path / "cur/t.csv", table({"V1": "y"}))
    before = {p.name: p.read_bytes() for p in base.iterdir()}
    run(base, tmp_path / "cur")
    assert {p.name: p.read_bytes() for p in base.iterdir()} == before
