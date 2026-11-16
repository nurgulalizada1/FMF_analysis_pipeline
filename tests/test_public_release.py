"""Public-release checks: what may be tracked, and what must stay private.

These run against the files tracked in Git, so they check the release itself
rather than a working directory that may legitimately hold private files.
"""
from __future__ import annotations

import csv
import fnmatch
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)


def tracked() -> list[str]:
    if shutil.which("git") is None or _git("rev-parse", "--git-dir").returncode != 0:
        pytest.skip("not a Git checkout")
    return _git("ls-files").stdout.splitlines()


def allowlist() -> list[str]:
    lines = (ROOT / "release_allowlist.txt").read_text().splitlines()
    return [x.strip() for x in lines if x.strip() and not x.lstrip().startswith("#")]


def test_every_tracked_file_is_on_the_allowlist():
    patterns = allowlist()
    outside = [f for f in tracked() if not any(fnmatch.fnmatchcase(f, p) for p in patterns)]
    assert outside == []


def test_no_data_or_results_are_tracked():
    for f in tracked():
        if f.startswith("data/"):
            name = Path(f).name
            assert name in ("README.md", ".gitkeep") or ".example." in name, f
        if f.startswith("results/"):
            assert Path(f).name == ".gitkeep", f


def test_no_binary_archive_or_office_file_is_tracked():
    forbidden = (".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".heic", ".xlsx", ".xls",
                 ".docx", ".doc", ".tsv", ".bundle", ".patch", ".zip", ".tgz", ".gz", ".7z",
                 ".tar", ".pkl", ".npy", ".parquet", ".db", ".sqlite")
    assert [f for f in tracked() if f.lower().endswith(forbidden)] == []


@pytest.mark.parametrize("path", [
    "study_parameters.json",
    "study_parameters.local.json",
    "data/raw/method.csv",
    "data/processed/patients.csv",
    "data/external/population_frequencies.csv",
    "data/external/source_annotations.csv",
    "data/reviewed/variant_annotations_reviewed.csv",
    "data/reviewed/annotation_decision_log.csv",
    "data/reviewed/accepted_differences.json",
    "data/synthetic/raw/synthetic_patients.csv",
    "results/tables/table2_variant_frequencies.csv",
    "results/run_manifest.json",
    "release.bundle",
    "changes.patch",
    "PRIVATE_FMF_study_files_DO_NOT_COMMIT.tgz",
    "archive.zip",
    "laboratory_report.pdf",
    "Screenshot 2026-01-01.png",
])
def test_private_paths_are_ignored(path):
    tracked()  # skips outside a Git checkout
    r = _git("check-ignore", "--no-index", "-q", path)
    assert r.returncode == 0, f"{path} is not ignored"


def test_templates_use_placeholders_only():
    """Templates must not carry real variant labels or study settings."""
    params = json.loads((ROOT / "study_parameters.example.json").read_text())
    labels = []
    for group in params["variant_groups"].values():
        labels += group
    plan = params["linkage_plan"]
    labels += plan["primary_pair"]
    for members in plan["haplotype_blocks"].values():
        labels += members
    for anchor, partners in plan["cross_locus_pairs"].items():
        labels += [anchor, *partners]
    labels += params["affinity_comparison_variants"]
    assert labels and all(x.startswith("EXAMPLE_") for x in labels), labels
    for f in [*(ROOT / "data").glob("*/*.example.csv")]:
        with open(f, encoding="utf-8-sig") as fh:
            rows = [r for r in csv.reader(fh, delimiter=";") if r and not r[0].startswith("#")]
        header, body = rows[0], rows[1:]
        if "variant" in header:
            i = header.index("variant")
            assert all(r[i].startswith("EXAMPLE_") for r in body), f.name
