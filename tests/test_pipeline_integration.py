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
    assert restored.exons_tested == ",".join(syn.SYNTHETIC_DESIGN["panel_regions"])


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


def test_allele_burden_denominator_is_the_whole_cohort(pipeline):
    t3 = read_csv(pipeline["first"] / "tables/table3_allele_burden.csv")
    assert (pd.to_numeric(t3.n, errors="coerce").dropna() <= N).all()


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
