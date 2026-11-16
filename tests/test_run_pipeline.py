"""Unit tests for the runner's output collection and clean-run logic."""
from __future__ import annotations

import importlib.util
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("run_pipeline", ROOT / "scripts/run_pipeline.py")
rp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rp)


def test_collect_excludes_manifests_and_hidden_files(tmp_path):
    tables = tmp_path / "tables"
    tables.mkdir()
    (tables / "t1.csv").write_text("a\n")
    (tables / "run_manifest.json").write_text("{}")   # an old manifest
    (tables / ".mefv-output").write_text("marker")
    manifest = tmp_path / "run_manifest.json"
    out, stale = rp.collect([tables], manifest, since=None)
    assert set(out) == {"tables/t1.csv"} and stale == []


def test_collect_reports_files_older_than_the_run_as_stale(tmp_path):
    tables = tmp_path / "tables"
    tables.mkdir()
    old = tables / "old.csv"
    old.write_text("x")
    past = time.time() - 3600
    import os
    os.utime(old, (past, past))
    (tables / "new.csv").write_text("y")
    out, stale = rp.collect([tables], tmp_path / "m.json", since=time.time() - 60)
    assert set(out) == {"tables/new.csv"} and stale == ["tables/old.csv"]


def test_combined_digest_depends_on_names_and_contents():
    a = {"x": {"sha256": "1"}, "y": {"sha256": "2"}}
    assert rp.combined(a) == rp.combined(dict(reversed(list(a.items()))))
    assert rp.combined(a) != rp.combined({"x": {"sha256": "1"}, "y": {"sha256": "3"}})


def test_clean_refuses_a_directory_it_did_not_create(tmp_path):
    d = tmp_path / "somewhere"
    d.mkdir()
    (d / "important.txt").write_text("keep me")
    with pytest.raises(SystemExit, match="refusing"):
        rp.clean_outputs([d])
    assert (d / "important.txt").exists()


def test_clean_empties_its_own_directory_but_keeps_hidden_files(tmp_path):
    d = tmp_path / "tables"
    rp.clean_outputs([d])                      # creates it, with a marker
    (d / "result.csv").write_text("x")
    (d / "sub").mkdir()
    (d / ".gitkeep").write_text("")
    rp.clean_outputs([d])
    assert sorted(p.name for p in d.iterdir()) == [".gitkeep", rp.MARKER]
