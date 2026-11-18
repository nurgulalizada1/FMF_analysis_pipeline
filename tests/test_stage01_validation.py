"""Stage 01 as a program: diagnostic-only mode, blocking, acknowledgement, safe output.

Runs scripts/01_build_dataset.py alone on a synthetic export that is corrupted
in one place per test. Synthetic data only; the tests check that the process
exit codes, the files written and the text printed are what the validation
policy promises (see docs/validation.md).
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "scripts" / "01_build_dataset.py"

_spec = importlib.util.spec_from_file_location(
    "make_synthetic_data", ROOT / "scripts" / "make_synthetic_data.py"
)
syn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(syn)


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory) -> Path:
    data = tmp_path_factory.mktemp("stage01") / "data"
    syn.write(data, syn.generate(n=200, seed=4242))
    return data


def _copy(synthetic: Path, dest: Path) -> Path:
    import shutil

    shutil.copytree(synthetic, dest)
    return dest


def _run(data: Path, out: Path, *args: str, **env) -> subprocess.CompletedProcess:
    e = os.environ.copy()
    for k in list(e):
        if k.startswith("MEFV_"):
            del e[k]
    e.update(
        MEFV_RAW_AUG=str(data / "raw/synthetic_method.csv"),
        MEFV_RAW_VARIANTS=str(data / "raw/synthetic_variants.csv"),
        MEFV_RAW_PATIENTS=str(data / "raw/synthetic_patients.csv"),
        MEFV_PROCESSED=str(out),
        MEFV_STUDY_PARAMETERS=str(data / "study_parameters.json"),
        **env,
    )
    return subprocess.run(
        [sys.executable, str(STAGE), *args], env=e, capture_output=True, text=True, cwd=ROOT
    )


def _append_variant_line(data: Path, line_from: int = 3, flip: bool = False) -> None:
    """Duplicate one data line of the variants file (``flip``: with the other zygosity)."""
    path = data / "raw/synthetic_variants.csv"
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    dup = lines[line_from]
    if flip:
        f = dup.split(";")
        f[3] = {"Het": "Hom", "Hom": "Het"}[f[3]]
        dup = ";".join(f)
    path.write_text("\n".join([*lines, dup]) + "\n", encoding="utf-8-sig")


def _report(out: Path) -> dict:
    return json.loads((out / "validation_report.json").read_text())


def test_clean_export_builds_and_writes_a_passing_report(synthetic, tmp_path):
    r = _run(synthetic, tmp_path / "o")
    assert r.returncode == 0, r.stderr[-500:]
    assert (tmp_path / "o/patients.csv").exists() and (tmp_path / "o/variants_long.csv").exists()
    assert _report(tmp_path / "o")["status"] == "passed"


def test_diagnostic_only_writes_only_the_report(synthetic, tmp_path):
    r = _run(synthetic, tmp_path / "o", "--diagnostic-only")
    assert r.returncode == 0
    assert sorted(p.name for p in (tmp_path / "o").iterdir()) == ["validation_report.json"]


def test_diagnostic_only_exits_2_on_findings_that_would_stop_the_run(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data)
    r = _run(data, tmp_path / "o", "--diagnostic-only")
    assert r.returncode == 2
    assert [p.name for p in (tmp_path / "o").iterdir()] == ["validation_report.json"]
    assert _report(tmp_path / "o")["status"] == "blocked"


def test_exact_duplicate_blocks_final_outputs_by_default(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data)
    r = _run(data, tmp_path / "o")
    assert r.returncode == 2
    assert "withheld" in r.stderr and "V5" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()
    assert _report(tmp_path / "o")["status"] == "blocked"


def test_conflicting_genotypes_block_final_outputs_by_default(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    r = _run(data, tmp_path / "o")
    assert r.returncode == 2 and not (tmp_path / "o/patients.csv").exists()
    rules = {f["rule"] for f in _report(tmp_path / "o")["findings"]}
    assert "V6b" in rules and "V5" not in rules


def _lowercase_one_label(data: Path) -> None:
    """V7: one row spells a variant label with different case from the others."""
    path = data / "raw/synthetic_variants.csv"
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    f = lines[3].split(";")
    f[1] = f[1].lower()
    lines[3] = ";".join(f)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8-sig")


def _disposition(tmp_path: Path, out: Path, rules: tuple[str, ...], **override) -> Path:
    """An owner disposition for the current findings, read from a diagnostic-only report."""
    rep = _report(out)
    for k in ("input_sha256", "config_sha256", "ruleset_sha256"):
        assert len(rep[k]) == 64                             # in the report so the owner can bind
    counts = {f["rule"]: f["records"] for f in rep["findings"] if f["severity"] == "blocking"}
    entries = [{"rule": r, "records": counts[r], "decision": "use_as_recorded",
                "decided_by": "LEAKNAME", "decided_on": "2026-01-02", "basis": "LEAKBASIS",
                **override} for r in rules]
    p = tmp_path / "private" / "disposition.json"
    p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps({"input_sha256": rep["input_sha256"],
                             "config_sha256": rep["config_sha256"],
                             "ruleset_sha256": rep["ruleset_sha256"], "dispositions": entries}))
    return p


CASES = {
    "duplicate": (lambda d: _append_variant_line(d), ("V5", "V6a")),
    "conflict": (lambda d: _append_variant_line(d, flip=True), ("V6b",)),
    "label_case": (_lowercase_one_label, ("V7",)),
}


@pytest.mark.parametrize("case", CASES)
def test_unacknowledged_state_stops_and_writes_no_analytical_output(synthetic, tmp_path, case):
    data = _copy(synthetic, tmp_path / "d")
    CASES[case][0](data)
    r = _run(data, tmp_path / "o")
    assert r.returncode == 2 and "withheld" in r.stderr
    assert sorted(p.name for p in (tmp_path / "o").iterdir()) == ["validation_report.json"]
    rep = _report(tmp_path / "o")
    assert rep["status"] == "blocked" and rep["outputs_released"] is False


@pytest.mark.parametrize("case", CASES)
def test_acknowledged_state_still_stops_and_writes_no_analytical_output(synthetic, tmp_path, case):
    """Acknowledging a finding must not release final outputs."""
    data = _copy(synthetic, tmp_path / "d")
    CASES[case][0](data)
    rules = CASES[case][1]
    r = _run(data, tmp_path / "o", "--accept-reviewed-findings", *rules)
    assert r.returncode == 2, "acknowledgement alone released the run"
    assert "Acknowledging does not release" in r.stderr
    assert sorted(p.name for p in (tmp_path / "o").iterdir()) == ["validation_report.json"]
    rep = _report(tmp_path / "o")
    assert rep["status"] == "acknowledged_unresolved" and rep["outputs_released"] is False
    assert rep["acknowledged_unresolved_rules"] == sorted(rules)
    assert rep["owner_dispositions"] == []


def test_acknowledgement_from_the_environment_does_not_release_either(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data)
    r = _run(data, tmp_path / "o", MEFV_ACCEPT_REVIEWED_FINDINGS="V5,V6a")
    assert r.returncode == 2
    assert not (tmp_path / "o/patients.csv").exists()


@pytest.mark.parametrize("case", CASES)
def test_an_explicit_owner_disposition_releases_and_is_recorded(synthetic, tmp_path, case):
    data = _copy(synthetic, tmp_path / "d")
    CASES[case][0](data)
    rules = CASES[case][1]
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", rules)
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 0, r.stderr[-500:]
    assert (tmp_path / "o/patients.csv").exists() and (tmp_path / "o/variants_long.csv").exists()
    rep = _report(tmp_path / "o")
    assert rep["status"] == "owner_disposition_recorded" and rep["outputs_released"] is True
    assert sorted(d["rule"] for d in rep["owner_dispositions"]) == sorted(rules)
    assert len(rep["owner_disposition_file_sha256"]) == 64
    # the conflict is still in the report, and nothing owner-written is
    assert {f["rule"] for f in rep["findings"]} >= set(rules)
    blob = r.stdout + r.stderr + json.dumps(rep)
    assert "LEAKNAME" not in blob and "LEAKBASIS" not in blob


def test_the_disposition_also_works_from_the_environment(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    r = _run(data, tmp_path / "o", MEFV_OWNER_DISPOSITION=str(disp))
    assert r.returncode == 0, r.stderr[-500:]


def test_a_disposition_with_a_stale_record_count_does_not_release(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    doc = json.loads(disp.read_text())
    doc["dispositions"][0]["records"] = 9999               # the owner reviewed a different count
    disp.write_text(json.dumps(doc))
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 2 and "different record count" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()


def test_the_export_changing_after_the_disposition_voids_it(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    _append_variant_line(data, line_from=5, flip=True)            # the export changes afterwards
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 2 and "different input files" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()


def _edit(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8-sig")
    assert old in text
    path.write_text(text.replace(old, new, 1), encoding="utf-8-sig")


def test_a_changed_input_with_the_same_finding_counts_does_not_inherit_the_disposition(
    synthetic, tmp_path
):
    """Counts alone are not a binding: any change to the reviewed files voids the disposition."""
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    ok = _run(data, tmp_path / "o0", "--owner-disposition", str(disp))
    assert ok.returncode == 0, ok.stderr[-500:]                    # the reviewed input releases
    # change one unrelated cell: the findings and their counts are exactly the same
    _edit(data / "raw/synthetic_patients.csv", ";male;", ";female;")
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 2 and "different input files" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()
    rep = _report(tmp_path / "o")
    assert rep["status"] == "blocked" and rep["outputs_released"] is False
    after = _report(tmp_path / "diag")["input_sha256"]
    assert rep["input_sha256"] != after


def test_a_reformatted_parameter_file_with_the_same_values_keeps_the_disposition(
    synthetic, tmp_path
):
    """The settings are bound by meaning, not by bytes: re-indenting or re-ordering is harmless."""
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    params = data / "study_parameters.json"
    doc = json.loads(params.read_text())
    params.write_text(json.dumps(dict(reversed(list(doc.items()))), indent=3))
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 0, r.stderr[-500:]


def test_a_changed_study_setting_voids_the_disposition(synthetic, tmp_path):
    """A different value in the study parameters is a different analysis: the disposition made
    under the old settings does not apply, though the data and the counts are unchanged."""
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    params = data / "study_parameters.json"
    doc = json.loads(params.read_text())
    doc["restored_individuals"] = {"coverage_confirmed": True}
    params.write_text(json.dumps(doc, indent=1))
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 2 and "different input files, settings or validation" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()
    rep = _report(tmp_path / "o")
    assert rep["status"] == "blocked" and rep["outputs_released"] is False
    assert rep["config_sha256"] != _report(tmp_path / "diag")["config_sha256"]


def test_a_disposition_made_under_other_validation_rules_does_not_release(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    doc = json.loads(disp.read_text())
    doc["ruleset_sha256"] = "7" * 64
    disp.write_text(json.dumps(doc))
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 2 and "different input files, settings or validation" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()


@pytest.mark.parametrize("key", ["input_sha256", "config_sha256", "ruleset_sha256"])
def test_a_disposition_without_a_binding_digest_is_refused(synthetic, tmp_path, key):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    doc = json.loads(disp.read_text())
    del doc[key]
    disp.write_text(json.dumps(doc))
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 2 and key in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()


def test_a_disposition_with_an_impossible_date_or_a_repeated_key_is_refused(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",), decided_on="2026-02-30")
    r = _run(data, tmp_path / "o", "--owner-disposition", str(disp))
    assert r.returncode == 2 and "real date" in r.stderr
    good = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    good.write_text(good.read_text().replace('"rule"', '"rule": "V6b", "rule"', 1))
    r = _run(data, tmp_path / "o2", "--owner-disposition", str(good))
    assert r.returncode == 2 and "repeated" in r.stderr
    assert not (tmp_path / "o2/patients.csv").exists()


def test_a_disposition_for_the_wrong_rule_does_not_release(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    wrong = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    entries = json.loads(wrong.read_text())
    entries["dispositions"][0]["rule"] = "V7"
    wrong.write_text(json.dumps(entries))
    r = _run(data, tmp_path / "o", "--owner-disposition", str(wrong))
    assert r.returncode == 2 and not (tmp_path / "o/patients.csv").exists()


def test_a_disposition_inside_the_repository_is_refused(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data)
    r = _run(data, tmp_path / "o", "--owner-disposition", str(ROOT / "owner_disposition.json"))
    assert r.returncode == 2 and "outside the repository" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()


def test_missing_or_invalid_disposition_file_fails_closed(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data)
    r = _run(data, tmp_path / "o1", "--owner-disposition", str(tmp_path / "absent.json"))
    assert r.returncode == 2 and not (tmp_path / "o1/patients.csv").exists()
    bad = tmp_path / "bad.json"
    bad.write_text("{}")
    r = _run(data, tmp_path / "o2", "--owner-disposition", str(bad))
    assert r.returncode == 2 and not (tmp_path / "o2/patients.csv").exists()


def test_diagnostic_only_never_applies_a_disposition(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    _run(data, tmp_path / "diag", "--diagnostic-only")
    disp = _disposition(tmp_path, tmp_path / "diag", ("V6b",))
    r = _run(data, tmp_path / "o", "--diagnostic-only", "--owner-disposition", str(disp))
    assert r.returncode == 2
    assert _report(tmp_path / "o")["status"] == "blocked"


def test_a_clean_export_needs_no_disposition_and_reports_none(synthetic, tmp_path):
    r = _run(synthetic, tmp_path / "o")
    rep = _report(tmp_path / "o")
    assert r.returncode == 0 and rep["outputs_released"] is True and rep["owner_dispositions"] == []


def test_a_fatal_rule_cannot_be_acknowledged(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    path = data / "raw/synthetic_variants.csv"
    text = path.read_text(encoding="utf-8-sig").replace(";Het;", ";het;", 1)
    path.write_text(text, encoding="utf-8-sig")
    r = _run(data, tmp_path / "o", "--accept-reviewed-findings", "V4")
    assert r.returncode == 2
    assert "cannot acknowledge" in r.stderr
    assert not (tmp_path / "o/patients.csv").exists()


def test_output_never_contains_identifiers_or_labels(synthetic, tmp_path):
    data = _copy(synthetic, tmp_path / "d")
    _append_variant_line(data, flip=True)
    path = data / "raw/synthetic_patients.csv"
    path.write_text(path.read_text(encoding="utf-8-sig") + "SYNTH_0001;20x1;male;1;yes;R1;x\n",
                    encoding="utf-8-sig")
    r = _run(data, tmp_path / "o")
    blob = r.stdout + r.stderr + (tmp_path / "o/validation_report.json").read_text()
    assert r.returncode == 2
    for token in ("SYNTH_", "SYN_", "20x1"):
        assert token not in blob, f"leaked {token}"


def test_a_restored_patient_is_counted_not_listed(synthetic, tmp_path):
    r = _run(synthetic, tmp_path / "o")
    assert "restored patients : <5" in r.stdout
    assert "SYNTH_" not in r.stdout


def test_stage_01_has_no_assert_statements_left():
    import ast

    tree = ast.parse(STAGE.read_text())
    assert not [n for n in ast.walk(tree) if isinstance(n, ast.Assert)]
