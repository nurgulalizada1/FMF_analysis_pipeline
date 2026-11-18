"""End-to-end tests: all nine stages on synthetic data.

Unit tests cover the library. These run the actual stage scripts in order,
from a synthetic raw export to the final tables, in a temporary directory, so
the script-level analysis logic is exercised as a whole. They need no study
data and run from a clean clone.

The pipeline is executed twice into the same output directory - the way a
user re-runs ``make demo`` - and the second run checks both determinism and
that a repeated run neither gains nor loses outputs.

The bootstrap behind the D' confidence intervals is the slowest step; these
tests use 200 replicates instead of the default 2000. The code path is the
same, and the setting is recorded in the run manifest.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from mefv import annotation as A
from mefv.utils import read_csv

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_pipeline.py"

_spec = importlib.util.spec_from_file_location(
    "make_synthetic_data", ROOT / "scripts" / "make_synthetic_data.py")
syn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(syn)

pytestmark = pytest.mark.integration

N = 200
SEED = 4242
REPLICATES = "200"

EXPECTED_TABLES = (
    "table1_cohort.csv", "table2_variant_frequencies.csv", "table3_allele_burden.csv",
    "table4_primary_pair.csv", "table4b_crosstab_primary.csv", "table5_complex_alleles.csv",
    "table6_gnomad_comparison.csv", "table7_population_affinity.csv",
    "table8_functional_concordance.csv", "table9_panel_performance.csv",
    "supp_hwe.csv", "supp_pairwise_associations.csv", "supp_haplotype_blocks.csv",
    "supp_high_allele_burden.csv",
)


def environment(data: Path, out: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.update({
        "MEFV_RAW_AUG": str(data / "raw/synthetic_method.csv"),
        "MEFV_RAW_VARIANTS": str(data / "raw/synthetic_variants.csv"),
        "MEFV_RAW_PATIENTS": str(data / "raw/synthetic_patients.csv"),
        "MEFV_PROCESSED": str(out / "processed"),
        "MEFV_SOURCE_ANNOTATIONS": str(data / "external/source_annotations.csv"),
        "MEFV_POPULATION_FREQUENCIES": str(data / "external/population_frequencies.csv"),
        "MEFV_ANNOTATION_OVERLAY": str(data / "reviewed/variant_annotations_reviewed.csv"),
        "MEFV_ANNOTATION_DECISION_LOG": str(data / "reviewed/annotation_decision_log.csv"),
        "MEFV_STUDY_PARAMETERS": str(data / "study_parameters.json"),
        "MEFV_ANNOTATIONS_GENERATED": str(out / "external/generated.csv"),
        "MEFV_DICT": str(out / "external/variant_dictionary.csv"),
        "MEFV_TABLES": str(out / "tables"),
        "MEFV_FIGURES": str(out / "figures"),
        "MEFV_ACMG": str(out / "classification"),
        "MEFV_RUN_MANIFEST": str(out / "run_manifest.json"),
        "MEFV_BOOTSTRAP_REPLICATES": REPLICATES,
        # Record every file the stages open, to show that no private default
        # location is read (see tests/open_audit/sitecustomize.py).
        "MEFV_OPEN_AUDIT_LOG": str(out.parent / f"opened_{out.name}.log"),
        "MEFV_OPEN_AUDIT_ROOT": str(ROOT),
        "PYTHONPATH": os.pathsep.join(
            [str(ROOT / "tests" / "open_audit"), env.get("PYTHONPATH", "")]),
    })
    # Make sure no developer's local study configuration can leak into the test.
    env.pop("SOURCE_DATE_EPOCH", None)
    return env


def run(data: Path, out: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(RUNNER), *extra], env=environment(data, out),
                          capture_output=True, text=True, cwd=ROOT)


def manifest(out: Path) -> dict:
    return json.loads((out / "run_manifest.json").read_text())


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    base = tmp_path_factory.mktemp("integration")
    data = base / "data"
    generated = syn.generate(n=N, seed=SEED)
    syn.write(data, generated)
    out = base / "run"
    r1 = run(data, out, "--clean")
    m1 = manifest(out) if r1.returncode == 0 else None
    snapshot = {sub: _digests(out / sub)
                for sub in ("tables", "figures", "processed", "external", "classification")}
    # Second run: same directory, no --clean - the repeated-execution case.
    r2 = run(data, out)
    m2 = manifest(out) if r2.returncode == 0 else None
    return {"data": data, "generated": generated, "out": out, "first": out,
            "r1": r1, "r2": r2, "m1": m1, "m2": m2, "snapshot": snapshot}


# --------------------------------------------------------------- execution --
def test_every_stage_succeeds(pipeline):
    r = pipeline["r1"]
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    stages = pipeline["m1"]["deterministic"]["stages"]
    assert [s["stage"][:2] for s in stages] == [f"{i:02d}" for i in range(1, 10)]
    assert all(s["returncode"] == 0 for s in stages)


def test_every_expected_table_is_produced(pipeline):
    tables = pipeline["first"] / "tables"
    missing = [t for t in EXPECTED_TABLES if not (tables / t).exists()]
    assert not missing


def test_figures_are_produced(pipeline):
    figures = sorted((pipeline["first"] / "figures").glob("*.png"))
    assert len(figures) >= 8


# ------------------------------------------------------------------ stage 01 --
def test_stage01_restores_corrects_and_resolves(pipeline):
    g = pipeline["generated"]
    pat = read_csv(pipeline["first"] / "processed/patients.csv")
    var = read_csv(pipeline["first"] / "processed/variants_long.csv")
    assert len(pat) == N
    assert pat.patient_id.notna().all() and (pat.patient_id != "").all()
    assert g["restored_id"] in set(pat.patient_id)
    assert syn.TYPO[0] not in set(var.variant)
    row = pat.set_index("patient_id").loc[g["override_id"]]
    assert row.assay_source == "laboratory_confirmed"
    # an override of a non-blank but wrong export value is labelled as an override too
    corrected = pat.set_index("patient_id").loc[g["corrected_id"]]
    assert corrected.assay_source == "laboratory_confirmed"
    assert corrected.assay == g["parameters"]["assay_overrides"][g["corrected_id"]]
    assert corrected.assay == syn.SYNTHETIC_DESIGN["assays"]["extended"]
    inference = syn.SYNTHETIC_DESIGN["assay_inference"]
    assert set(pat.assay_source) == {
        "recorded", inference["source_label"], "laboratory_confirmed"}
    inferred = pat[pat.assay_source == inference["source_label"]]
    assert (inferred.year_requested.astype(int) <= inference["until_year"]).all()
    assert (inferred.assay == inference["assay"]).all()
    restored = pat.set_index("patient_id").loc[g["restored_id"]]
    # method-only restored individuals have no laboratory testing record: they stay
    # undetermined (blank exons_tested, undetermined detection) until the owner confirms
    assert restored.exons_tested == ""
    assert restored.variant_detected == "undetermined"
    assert restored.coverage_source == "restored_from_method_coverage_unconfirmed"


# ------------------------------------------------------------------ stage 02 --
def test_stage02_dictionary_matches_the_schema(pipeline):
    path = pipeline["first"] / "external/variant_dictionary.csv"
    d = pd.read_csv(path, sep=";", dtype=str, keep_default_na=False, encoding="utf-8-sig")
    assert len(d.columns) == 17
    A.validate(d)


# ------------------------------------------------------------- analysis -----
def test_outputs_contain_only_synthetic_variant_labels(pipeline):
    """Guards against study variant names hard-coded in analysis code."""
    t2 = read_csv(pipeline["first"] / "tables/table2_variant_frequencies.csv")
    assert t2.variant.str.startswith("SYN_").all(), sorted(set(t2.variant))


def test_primary_pair_is_read_from_parameters_and_shows_linkage(pipeline):
    t4 = read_csv(pipeline["first"] / "tables/table4_primary_pair.csv")
    a, b = syn.LD_PAIR
    assert set(t4.variant_pair) == {f"{a} / {b}"}
    primary = t4[t4.stratum.str.contains("primary")]
    assert float(primary.Dprime.iloc[0]) > 0.5


def test_the_invented_complex_allele_is_detected(pipeline):
    t5 = read_csv(pipeline["first"] / "tables/table5_complex_alleles.csv")
    pairs = {frozenset((r.variant_A, r.variant_B)) for r in t5.itertuples()}
    assert frozenset(syn.COMPLEX_PAIR) in pairs


def test_periods_and_assays_come_from_the_synthetic_design(pipeline):
    design = syn.SYNTHETIC_DESIGN
    t1 = read_csv(pipeline["first"] / "tables/table1_cohort.csv")
    period = t1[t1.item == "Period"].iloc[0]
    assert period.era_early.split("–")[0] == str(design["periods"]["primary"]["first_year"])
    t4 = read_csv(pipeline["first"] / "tables/table4_primary_pair.csv")
    labels = set(t4.stratum)
    assert f"{design['periods']['primary']['label']} (primary)" in labels
    assert {design["assays"]["panel"], design["assays"]["extended"]} <= labels


def test_panel_restriction_excludes_extended_only_regions(pipeline):
    t2 = read_csv(pipeline["first"] / "tables/table2_variant_frequencies.csv")
    assert not set(t2.variant) & set(syn.EXTENDED_ONLY)
    assert set(t2.exon) <= set(syn.SYNTHETIC_DESIGN["panel_regions"])


def test_allele_burden_counts_only_panel_eligible_individuals(pipeline):
    """Burden needs every panel region tested; undetermined individuals are excluded."""
    t3 = read_csv(pipeline["first"] / "tables/table3_allele_burden.csv")
    n = pd.to_numeric(t3.n, errors="coerce").dropna()
    pat = read_csv(pipeline["first"] / "processed/patients.csv")
    panel = set(syn.SYNTHETIC_DESIGN["panel_regions"])
    eligible = sum(1 for v in pat.exons_tested if v and panel <= set(v.split(",")))
    assert eligible < N                       # the restored individual is not eligible
    assert (n <= eligible).all() and n.max() > 0


def test_table2_denominators_are_internally_consistent(pipeline):
    t2 = read_csv(pipeline["first"] / "tables/table2_variant_frequencies.csv")
    num = {c: pd.to_numeric(t2[c]) for c in (
        "n_cohort", "n_in_scope", "n_tested", "n_not_tested", "n_undetermined", "n_eligible",
        "allele_denominator", "heterozygous", "homozygous", "carriers", "alleles")}
    assert (num["n_cohort"] == N).all()
    assert (num["n_tested"] + num["n_not_tested"] + num["n_undetermined"]
            == num["n_in_scope"]).all()
    assert (num["n_eligible"] == num["n_tested"]).all()
    assert (num["allele_denominator"] == 2 * num["n_eligible"]).all()
    assert (num["carriers"] == num["heterozygous"] + num["homozygous"]).all()
    assert (num["alleles"] == num["heterozygous"] + 2 * num["homozygous"]).all()
    ok = t2.eligibility_status == "evaluable"
    af = pd.to_numeric(t2.allele_freq)
    assert ((af[ok] - num["alleles"][ok] / num["allele_denominator"][ok]).abs() < 5e-5).all()
    assert af[~ok].isna().all()               # never zero when not evaluable


def test_restored_individual_is_undetermined_in_every_denominator(pipeline):
    """Method-only restored individuals have no testing record until the owner confirms."""
    pat = read_csv(pipeline["first"] / "processed/patients.csv")
    restored = pat[pat.coverage_source == "restored_from_method_coverage_unconfirmed"]
    assert len(restored) == 1 and restored.exons_tested.iloc[0] == ""
    t2 = read_csv(pipeline["first"] / "tables/table2_variant_frequencies.csv")
    panel = t2[t2.exon.isin(syn.SYNTHETIC_DESIGN["panel_regions"])]
    assert (pd.to_numeric(panel.n_undetermined) >= 1).all()
    assert (pd.to_numeric(panel.n_tested) <= N - 1).all()
    t1 = read_csv(pipeline["first"] / "tables/table1_cohort.csv")
    und = t1[t1.item.str.contains("undetermined", case=False)]
    assert (und.total.astype(str).str.split().str[0] == "1").any()


def test_validation_report_is_recorded_and_passed_on_clean_synthetic_data(pipeline):
    rep = json.loads((pipeline["first"] / "processed/validation_report.json").read_text())
    assert rep["status"] == "passed" and rep["acknowledged_rules"] == []
    recorded = pipeline["m1"]["deterministic"]["intermediates"]
    assert any("validation_report.json" in k for k in recorded)


def test_period_restricted_variants_use_the_primary_period_scope(pipeline):
    scope = read_csv(pipeline["first"] / "tables/supp_period_scope.csv")
    assert len(scope) > 0 and set(scope.restriction_metadata) == {"variant_groups.synonymous"}
    t2 = read_csv(pipeline["first"] / "tables/table2_variant_frequencies.csv")
    restricted = t2[t2.scope != "whole cohort"]
    assert (pd.to_numeric(restricted.n_in_scope) < N).all()


def test_bootstrap_outputs_report_replicate_accounting(pipeline):
    t4 = read_csv(pipeline["first"] / "tables/table4_primary_pair.csv")
    cols = {c for c in t4.columns if "valid" in c or "dropped" in c}
    assert cols, sorted(t4.columns)


def test_reference_table_keeps_missing_distinct_from_zero(pipeline):
    t6 = read_csv(pipeline["first"] / "tables/table6_gnomad_comparison.csv")
    for col in ("x_vs_highest_status", "x_vs_mid_status", "x_vs_global_status"):
        assert set(t6[col]) <= {"computed", "undefined_reference_zero",
                                "not_evaluable_reference_missing",
                                "not_evaluable_cohort_frequency_undefined"}
    assert not t6.isin(["inf", "Inf", "-inf", "infinity"]).any().any()   # no infinite enrichment


def test_acmg_provenance_is_written_beside_the_evidence_package(pipeline):
    prov = pd.read_csv(pipeline["first"] / "classification/acmg_rule_provenance.tsv", sep="\t",
                       dtype=str, keep_default_na=False, encoding="utf-8-sig")
    assert {"BA1_threshold", "BS1_threshold", "classification"} <= set(prov.item)
    assert prov.loc[prov.item == "classification", "value"].iloc[0] == "none"


def test_no_output_contains_privacy_scanner_findings(pipeline):
    from mefv import privacy
    items = privacy.path_items([pipeline["first"] / d for d in ("tables", "classification")],
                               ROOT)
    res = privacy.scan_items(items, None, path_rules=False)
    assert res.findings == [], res.render()


# -------------------------------------------------------------- determinism --
def _digests(root: Path) -> dict[str, str]:
    out = {}
    for f in sorted(root.rglob("*")):
        if f.is_file() and not f.name.startswith("."):
            out[str(f.relative_to(root))] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


def test_manifest_is_outside_the_outputs_and_not_self_collected(pipeline):
    m = pipeline["m1"]
    keys = set(m["deterministic"]["artifacts"])
    assert not any(k.endswith("run_manifest.json") for k in keys)
    assert not (pipeline["out"] / "tables" / "run_manifest.json").exists()
    assert set(m) == {"manifest_version", "deterministic", "execution"}
    assert "generated_at" in m["execution"] and "generated_at" not in m["deterministic"]


def test_repeated_run_is_byte_identical_and_stable_in_count(pipeline):
    assert pipeline["r2"].returncode == 0, pipeline["r2"].stderr[-2000:]
    m1, m2 = pipeline["m1"], pipeline["m2"]
    # the whole deterministic section, not just the digests, is identical
    assert m1["deterministic"] == m2["deterministic"]
    assert len(m2["deterministic"]["artifacts"]) == len(m1["deterministic"]["artifacts"])
    assert m2["execution"]["stale_files_not_collected"] == []
    for sub, before in pipeline["snapshot"].items():
        after = _digests(pipeline["out"] / sub)
        assert after.keys() == before.keys(), sub
        assert [k for k in before if before[k] != after[k]] == [], sub


def test_every_artifact_in_the_directories_is_in_the_manifest(pipeline):
    on_disk = {f"{sub}/{k}" for sub in ("tables", "figures", "classification")
               for k in _digests(pipeline["out"] / sub)}
    assert on_disk == set(pipeline["m2"]["deterministic"]["artifacts"])


def test_no_private_default_location_is_read(pipeline):
    """The demo must not depend on, or read, any private file in the checkout."""
    log = pipeline["out"].parent / "opened_run.log"
    opened = set(log.read_text().splitlines()) if log.exists() else set()
    assert opened, "the open audit recorded nothing; it is not active"
    private = [ROOT / "study_parameters.json", ROOT / "data"]
    hits = sorted(p for p in opened
                  if any(Path(p).is_relative_to(q) for q in private))
    assert hits == []


# -------------------------------------------------------------- failure modes -
def test_missing_study_parameters_stop_the_pipeline(pipeline, tmp_path):
    env = environment(pipeline["data"], tmp_path)
    env["MEFV_STUDY_PARAMETERS"] = str(tmp_path / "absent.json")
    r = subprocess.run([sys.executable, str(RUNNER), "--only", "01"], env=env,
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode != 0
    assert "study parameters not found" in (r.stdout + r.stderr)


def test_missing_annotation_overlay_stops_the_pipeline(pipeline, tmp_path):
    env = environment(pipeline["data"], tmp_path)
    env["MEFV_PROCESSED"] = str(pipeline["first"] / "processed")
    env["MEFV_ANNOTATION_OVERLAY"] = str(tmp_path / "absent.csv")
    r = subprocess.run([sys.executable, str(RUNNER), "--only", "02"], env=env,
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode != 0
    assert "reviewed annotation overlay not found" in (r.stdout + r.stderr)


# ----------------------------------------- blocking findings, end to end ------
def _conflicting_export(pipeline, tmp_path):
    import shutil
    data = tmp_path / "data"
    shutil.copytree(pipeline["data"], data)
    path = data / "raw/synthetic_variants.csv"
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    f = lines[3].split(";")
    f[3] = {"Het": "Hom", "Hom": "Het"}[f[3]]
    path.write_text("\n".join([*lines, ";".join(f)]) + "\n", encoding="utf-8-sig")
    return data


def test_acknowledged_blocking_finding_stops_the_whole_pipeline(pipeline, tmp_path):
    data = _conflicting_export(pipeline, tmp_path)
    out = tmp_path / "run"
    r = run(data, out, "--clean", "--accept-reviewed-findings", "V6b")
    assert r.returncode != 0
    assert "Acknowledging does not release" in r.stderr
    m = manifest(out)["deterministic"]                    # a failed run is recorded as failed
    assert [(x["stage"][:2], x["returncode"] != 0) for x in m["stages"]] == [("01", True)]
    assert m["artifacts"] == {}
    visible = [f for sub in ("tables", "figures", "classification") for f in (out / sub).glob("*")
               if not f.name.startswith(".")]
    assert [f.name for f in visible] == ["RUN_FAILED_outputs_are_not_current.txt"] * 3
    assert not (out / "processed/patients.csv").exists()


def _disposition_file(tmp_path, data, rule):
    """Review the export with a diagnostic stage 01 run and write the matching owner disposition
    (in a folder outside any output directory), exactly as an owner would."""
    diag = run(data, tmp_path / "diag", "--only", "01")           # blocked: writes only the report
    assert diag.returncode != 0
    rep = json.loads((tmp_path / "diag/processed/validation_report.json").read_text())
    n = next(f["records"] for f in rep["findings"] if f["rule"] == rule)   # as shown: maybe "<5"
    disp = tmp_path / "private" / "disposition.json"
    disp.parent.mkdir(exist_ok=True)
    disp.write_text(json.dumps({
        "input_sha256": rep["input_sha256"], "config_sha256": rep["config_sha256"],
        "ruleset_sha256": rep["ruleset_sha256"],
        "dispositions": [{"rule": rule, "records": n, "decision": "use_as_recorded",
                          "decided_by": "synthetic owner", "decided_on": "2026-01-02",
                          "basis": "synthetic test"}]}))
    return disp


def _label_case_export(pipeline, tmp_path):
    """One variant label spelled in two cases (V7). Each individual still has one genotype
    per label, so no patient-level genotype is ambiguous."""
    import shutil
    data = tmp_path / "data"
    shutil.copytree(pipeline["data"], data)
    path = data / "raw/synthetic_variants.csv"
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    f = lines[3].split(";")
    f[1] = f[1].lower()
    lines[3] = ";".join(f)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8-sig")
    return data


def test_owner_disposition_runs_the_pipeline_and_flags_the_outputs(pipeline, tmp_path):
    data = _label_case_export(pipeline, tmp_path)
    disp = _disposition_file(tmp_path, data, "V7")
    out = tmp_path / "run"
    r = run(data, out, "--clean", "--owner-disposition", str(disp))
    assert r.returncode == 0, r.stderr[-1500:]
    t1 = read_csv(out / "tables/table1_cohort.csv")
    assert (t1.total == "owner_disposition_recorded").any()
    assert t1.item.str.contains("Owner disposition V7").any()
    final = json.loads((out / "processed/validation_report.json").read_text())
    assert final["status"] == "owner_disposition_recorded" and final["outputs_released"] is True
    assert manifest(out)["deterministic"]["run_complete"] is True
    # the validation report, with the disposition digest, is part of the manifest
    recorded = manifest(out)["deterministic"]["intermediates"]
    assert any("validation_report.json" in k for k in recorded)


def test_a_dispositioned_genotype_conflict_still_blocks_patient_level_genotype_stages(
    pipeline, tmp_path
):
    """'use_as_recorded' lets stage 01 release the export and the stages that count rows. It
    does not say which of two genotypes an individual has, so the stage that needs one genotype
    per individual and variant (linkage) refuses to run, names the rules, and the run is marked
    failed. Nothing is chosen, merged or dropped on the owner's behalf."""
    data = _conflicting_export(pipeline, tmp_path)
    disp = _disposition_file(tmp_path, data, "V6b")
    out = tmp_path / "run"
    r = run(data, out, "--clean", "--owner-disposition", str(disp))
    assert r.returncode != 0
    assert "V5/V6a/V6b" in r.stderr and "stage 04_linkage.py failed" in r.stderr
    m = manifest(out)["deterministic"]
    assert [x["returncode"] == 0 for x in m["stages"]] == [True, True, True, False]
    assert m["run_complete"] is False
    assert (out / "tables" / "RUN_FAILED_outputs_are_not_current.txt").exists()
    final = json.loads((out / "processed/validation_report.json").read_text())
    assert final["status"] == "owner_disposition_recorded"   # stage 01 itself was released
    t1 = read_csv(out / "tables/table1_cohort.csv")          # the counting stage ran, as recorded
    assert t1.item.str.contains("Owner disposition V6b").any()


# ------------------------------- nothing stale is current after a blocked run ---
def _copy_successful_run(pipeline, tmp_path):
    import shutil
    out = tmp_path / "run"
    shutil.copytree(pipeline["out"], out)
    return out


def test_a_blocked_rerun_leaves_no_output_that_can_pass_for_current(pipeline, tmp_path):
    out = _copy_successful_run(pipeline, tmp_path)
    assert (out / "processed/patients.csv").exists()            # the earlier run is complete
    bad = _conflicting_export(pipeline, tmp_path)
    r = run(bad, out)                                            # re-run, no --clean
    assert r.returncode != 0
    # the superseded processed tables are withdrawn, and the report says nothing is released
    assert not (out / "processed/patients.csv").exists()
    assert not (out / "processed/variants_long.csv").exists()
    rep = json.loads((out / "processed/validation_report.json").read_text())
    assert rep["outputs_released"] is False and rep["status"] == "blocked"
    # old tables and figures are still on disk, but every output directory says they are not current
    from importlib import util
    spec = util.spec_from_file_location("run_pipeline", RUNNER)
    rp = util.module_from_spec(spec)
    spec.loader.exec_module(rp)
    for sub in ("tables", "figures", "classification"):
        assert (out / sub / rp.FAILED_NAME).exists(), sub
    assert manifest(out)["deterministic"]["run_complete"] is False
    # no later stage can be run on what is left
    for stage in ("03", "04", "05", "06", "08", "09"):
        again = run(bad, out, "--only", stage)
        assert again.returncode != 0, f"stage {stage} ran after a blocked run"
    # a good export and a full run restore a current, complete state and remove the notes
    r2 = run(pipeline["data"], out, "--clean")
    assert r2.returncode == 0, r2.stderr[-800:]
    assert not list(out.rglob(rp.FAILED_NAME))
    assert manifest(out)["deterministic"]["run_complete"] is True


def test_later_stages_refuse_tables_that_stage_01_did_not_release(pipeline, tmp_path):
    """The tables of a successful run are still on disk; a diagnostic-only stage 01 run withdraws
    their release, so every stage that reads them must refuse."""
    out = _copy_successful_run(pipeline, tmp_path)
    r = subprocess.run([sys.executable, str(ROOT / "scripts/01_build_dataset.py"),
                        "--diagnostic-only"], env=environment(pipeline["data"], out),
                       capture_output=True, text=True, cwd=ROOT)
    assert r.returncode == 0, r.stderr[-500:]                    # the export itself is clean
    rep = json.loads((out / "processed/validation_report.json").read_text())
    assert rep["diagnostic_only"] is True and rep["outputs_released"] is False
    assert (out / "processed/patients.csv").exists()             # still there ...
    for stage in ("03", "04", "05", "06", "08", "09"):           # ... and refused
        again = run(pipeline["data"], out, "--only", stage)
        assert again.returncode != 0, f"stage {stage} accepted unreleased tables"
        assert "not released" in again.stderr, stage


@pytest.mark.parametrize("damage", ["missing", "empty_object", "malformed", "not_true"])
def test_later_stages_refuse_tables_without_a_valid_release_record(pipeline, tmp_path, damage):
    """Fail closed: the processed tables are on disk, but a missing, malformed or non-authorising
    validation report must stop every stage that reads them."""
    out = _copy_successful_run(pipeline, tmp_path)
    report = out / "processed/validation_report.json"
    assert json.loads(report.read_text())["outputs_released"] is True
    if damage == "missing":
        report.unlink()
    else:
        report.write_text({"empty_object": "{}", "malformed": "{oops",
                           "not_true": '{"outputs_released": "true"}'}[damage])
    assert (out / "processed/patients.csv").exists()
    for stage in ("03", "04", "05", "06", "08", "09"):
        again = run(pipeline["data"], out, "--only", stage)
        assert again.returncode != 0, f"stage {stage} accepted tables with a {damage} report"
        assert "not released" in again.stderr, (stage, damage)


def test_a_missing_disposition_file_withdraws_earlier_outputs(pipeline, tmp_path):
    out = _copy_successful_run(pipeline, tmp_path)
    r = run(pipeline["data"], out, "--owner-disposition", str(tmp_path / "absent.json"))
    assert r.returncode != 0
    assert not (out / "processed/patients.csv").exists()
    assert not (out / "processed/validation_report.json").exists()


# ------------------- failure markers and partial runs ---
def _failed_name():
    from importlib import util
    spec = util.spec_from_file_location("run_pipeline", RUNNER)
    rp = util.module_from_spec(spec)
    spec.loader.exec_module(rp)
    return rp.FAILED_NAME


def _failed_state(pipeline, tmp_path):
    """A complete earlier run, then a blocked full run: leftovers plus failure markers."""
    out = _copy_successful_run(pipeline, tmp_path)
    bad = _conflicting_export(pipeline, tmp_path)
    assert run(bad, out).returncode != 0
    name = _failed_name()
    for sub in ("tables", "figures", "classification", "processed"):
        assert (out / sub / name).exists(), sub
    return out, bad, name


def test_a_successful_partial_run_does_not_clear_failure_markers(pipeline, tmp_path):
    out, _bad, name = _failed_state(pipeline, tmp_path)
    r = run(pipeline["data"], out, "--only", "01")                # good export, stage 01 succeeds
    assert r.returncode == 0, r.stderr[-800:]
    for sub in ("tables", "figures", "classification", "processed"):
        assert (out / sub / name).exists(), f"{sub}: a partial run removed the failure marker"
    m = manifest(out)
    assert m["deterministic"]["run_complete"] is False
    assert m["execution"]["partial_run"] is True


def test_a_failed_partial_run_keeps_the_markers_and_is_not_complete(pipeline, tmp_path):
    out, bad, name = _failed_state(pipeline, tmp_path)
    r = run(bad, out, "--only", "01")                             # blocked again
    assert r.returncode != 0
    for sub in ("tables", "figures", "classification", "processed"):
        assert (out / sub / name).exists(), sub
    assert not (out / "processed/patients.csv").exists()
    m = manifest(out)["deterministic"]
    assert m["run_complete"] is False
    assert [x["returncode"] != 0 for x in m["stages"]] == [True]


def test_a_partial_run_starting_late_keeps_the_markers(pipeline, tmp_path):
    out, _bad, name = _failed_state(pipeline, tmp_path)
    r = run(pipeline["data"], out, "--start", "09")
    assert r.returncode != 0                      # stage 01 never re-released the tables
    assert all((out / sub / name).exists() for sub in ("tables", "figures", "classification"))
    assert manifest(out)["deterministic"]["run_complete"] is False


def test_a_full_successful_run_clears_failure_markers(pipeline, tmp_path):
    out, _bad, name = _failed_state(pipeline, tmp_path)
    r = run(pipeline["data"], out)                                # full run, no --clean
    assert r.returncode == 0, r.stderr[-800:]
    assert not list(out.rglob(name))
    m = manifest(out)
    assert m["deterministic"]["run_complete"] is True
    assert m["execution"]["partial_run"] is False


# ----------------------------------------------- degenerate inputs, end to end ---
def _edit_rows(path, fn):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    head, rows = lines[:2], lines[2:]
    path.write_text("\n".join([*head, *fn(rows)]) + "\n", encoding="utf-8-sig")


def _degenerate(pipeline, tmp_path, kind):
    import shutil
    data = tmp_path / "data"
    shutil.copytree(pipeline["data"], data)
    if kind == "no_variant_rows":
        _edit_rows(data / "raw/synthetic_variants.csv", lambda rows: [])

        def clear(rows):
            out = []
            for r in rows:
                f = r.split(";")
                f[3], f[4] = "0", "no"
                out.append(";".join(f))
            return out

        _edit_rows(data / "raw/synthetic_patients.csv", clear)
    elif kind == "nobody_eligible":
        def blank_exons(rows):
            out = []
            for r in rows:
                f = r.split(";")
                f[5] = ""
                out.append(";".join(f))
            return out

        _edit_rows(data / "raw/synthetic_patients.csv", blank_exons)
    elif kind == "one_variant_row":
        _edit_rows(data / "raw/synthetic_variants.csv", lambda rows: rows[:1])
    return data


@pytest.mark.parametrize("kind", ["no_variant_rows", "nobody_eligible", "one_variant_row"])
def test_degenerate_inputs_finish_and_report_undefined_values_as_missing(
    pipeline, tmp_path, kind
):
    """Every stage must cope with empty groups and empty tables: no crash, no figure that fails
    to draw, no zero or p = 1 standing in for something that cannot be computed."""
    data = _degenerate(pipeline, tmp_path, kind)
    out = tmp_path / "run"
    r = run(data, out, "--clean")
    assert r.returncode == 0, r.stderr[-1500:]
    assert manifest(out)["deterministic"]["run_complete"] is True
    assert len(list((out / "figures").glob("*.png"))) >= 8
    t3 = read_csv(out / "tables/table3_allele_burden.csv")
    empty = t3[t3.total == "0"]
    assert (empty.percent == "").all() and (empty.ci_95 == "").all()     # no percentage of nothing
    high = read_csv(out / "tables/supp_high_allele_burden.csv")
    assert list(high.columns) == ["patient_id", "year", "assay", "total_alleles", "genotype"]
    if kind == "nobody_eligible":
        t2 = read_csv(out / "tables/table2_variant_frequencies.csv")
        assert (t2.allele_freq == "nan").all() or (t2.allele_freq == "").all()
        assert set(t2.eligibility_status) == {"not_evaluable_no_eligible_individuals"}
        assert (read_csv(out / "tables/table9_panel_performance.csv")).empty
    if kind == "no_variant_rows":
        assert (read_csv(out / "tables/table9_panel_performance.csv")).empty


def test_non_convergence_and_bootstrap_counts_are_reported_in_every_ld_table(pipeline):
    tables = pipeline["out"] / "tables"
    t4 = read_csv(tables / "table4_primary_pair.csv")
    assert {"em_converged", "boot_replicates_requested", "boot_replicates_valid",
            "boot_replicates_dropped"} <= set(t4.columns)
    assert (t4.boot_replicates_requested.astype(int) == int(REPLICATES)).all()
    assert (t4.boot_replicates_valid.astype(int) + t4.boot_replicates_dropped.astype(int)
            == t4.boot_replicates_requested.astype(int)).all()
    for name in ("supp_pairwise_associations.csv", "supp_haplotype_blocks.csv"):
        t = read_csv(tables / name)
        assert "em_converged" in t.columns, name
        assert set(t.em_converged) <= {"True", "False", "not_run"}, name
