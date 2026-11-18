"""Regression tests for the finalization blockers.

Each test names the failure it guards against. All data are invented.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from mefv import DataIntegrityError, allele_dose, genotype_vector, privacy
from mefv import comparison as C
from mefv import disclosure as D
from mefv import validation as V
from mefv.eligibility import panel_eligible, record_state, region_sets
from mefv.statistics import bootstrap_dprime, ld_stats


# ------------------------------------------------------------ small counts ----
@pytest.mark.parametrize(
    ("n", "expected"), [(0, "0"), (1, "<5"), (4, "<5"), (5, "5"), (120, "120")]
)
def test_counts_below_five_are_not_shown(n, expected):
    assert D.text(n) == expected
    assert D.shown(n) == (D.SUPPRESSED if 0 < n < 5 else n)


def test_suppress_frame_and_value_counts_keep_big_counts_and_text():
    df = pd.DataFrame({"a": [1, 7, np.nan], "b": ["x", 2, 3.5], "c": [True, False, True]})
    out = D.suppress_frame(df)
    assert out["a"].tolist()[:2] == ["<5", 7] and out["b"].tolist() == ["x", "<5", 3.5]
    assert out["c"].tolist() == [True, False, True]            # booleans are not counts
    vc = D.suppressed_counts(pd.Series(["a"] * 6 + ["b"] * 2).value_counts())
    assert vc.to_dict() == {"a": 6, "b": "<5"}


def test_validation_report_text_json_and_errors_show_no_small_count():
    t = pd.DataFrame({"patient_id": ["P1", "P2"], "variant": ["V", "V"]})
    rep = V.ValidationReport(tables={"variants_rows": 3})
    rep.add("V5", records=3, individuals=2, detail={"surplus_rows": 1})
    rep.add("V6b", records=12, individuals=6, detail={"individual_variant_pairs": 6})
    d = rep.to_dict()
    by = {f["rule"]: f for f in d["findings"]}
    assert by["V5"]["records"] == "<5" and by["V5"]["individuals"] == "<5"
    assert by["V5"]["detail"] == {"surplus_rows": "<5"}
    assert by["V6b"]["records"] == 12 and by["V6b"]["individuals"] == 6
    assert d["tables"]["variants_rows"] == "<5"
    text = rep.render()
    assert "records=<5" in text and "records=12" in text and "surplus_rows=<5" in text
    with pytest.raises(DataIntegrityError) as e:
        rep.enforce()
    assert "V5: <5 record" in str(e.value) and "V6b: 12 record" in str(e.value)
    assert t is not None


def test_error_messages_from_helpers_suppress_small_counts():
    with pytest.raises(DataIntegrityError, match=r"<5 entr"):
        allele_dose(pd.Series(["Het", "x", "Hom", "y"]))
    with pytest.raises(DataIntegrityError, match=r"7 entr"):
        allele_dose(pd.Series(["x"] * 7))
    from mefv.statistics import em_haplotypes
    with pytest.raises(ValueError, match=r"<5 genotype"):
        em_haplotypes([0, 3, 1], [0, 1, 1])


# ------------------------------------------------- duplicate patient-variant ----
def _var(rows):
    return pd.DataFrame(rows, columns=["patient_id", "variant", "zygosity"])


def test_genotype_vector_one_row_per_individual_is_unchanged():
    var = _var([("P1", "V", "Het"), ("P2", "V", "Hom"), ("P3", "W", "Het")])
    assert genotype_vector(var, "V", ["P1", "P2", "P3", "P4"]) == [1, 2, 0, 0]


@pytest.mark.parametrize("second", ["Het", "Hom"], ids=["same genotype repeated", "conflicting"])
def test_genotype_vector_refuses_a_second_row_for_the_same_individual_and_variant(second):
    """The earlier code silently kept the last row, so a conflict changed the result with no
    trace and a repeated row hid a data problem."""
    var = _var([("LEAK_P1", "V", "Het"), ("LEAK_P1", "V", second), ("P2", "V", "Het")])
    with pytest.raises(DataIntegrityError) as e:
        genotype_vector(var, "V", ["LEAK_P1", "P2"])
    msg = str(e.value)
    assert "V5/V6a/V6b" in msg and "LEAK_P1" not in msg and "<5" in msg


def test_a_duplicate_on_another_variant_does_not_block_this_one():
    var = _var([("P1", "W", "Het"), ("P1", "W", "Hom"), ("P2", "V", "Het")])
    assert genotype_vector(var, "V", ["P1", "P2"]) == [0, 1]


def test_genotype_vector_refuses_a_repeated_requested_identifier():
    var = _var([("P1", "V", "Het")])
    with pytest.raises(DataIntegrityError, match="more than once"):
        genotype_vector(var, "V", ["P1", "P1"])


def test_genotype_vector_still_refuses_unrecognised_zygosity():
    with pytest.raises(DataIntegrityError):
        genotype_vector(_var([("P1", "V", "het")]), "V", ["P1"])


# ----------------------------------------------------- canonical identifiers ----
@pytest.mark.parametrize(("raw", "expected"), [
    ("P12", "P012"), (" p12 ", "P012"), ("12", "P012"), ("P0012", "P012"), (7, "P007"),
    ("P12345", "P12345"),                                  # longer than the width: kept whole
    ("P12P3", None),                                       # embedded prefix: not merged into 123
    ("PP12", None), ("P1P", None), ("1P2", None), ("P-12", None), ("P 12", None),
    ("P+12", None), ("P12.0", None), ("P١٢", None), ("P12\n3", None),
    ("", None), ("P", None), (None, None), (float("nan"), None), (pd.NA, None),
])
def test_canonical_id_never_repairs_or_truncates(raw, expected):
    assert V.canonical_id(raw, "P", 3) == expected


def test_embedded_prefix_cannot_make_two_identifiers_collide():
    """'P12P3' used to lose its inner prefix and become 123, the same as 'P123'."""
    assert V.canonical_id("P12P3", "P", 3) != V.canonical_id("P123", "P", 3)
    assert V.canonical_id("P123", "P", 3) == "P123"


def test_a_missing_exons_tested_is_unusable_not_the_region_named_nan():
    for missing in (None, float("nan"), pd.NA, np.nan):
        assert V.parse_region_set(missing) is None
    assert V.parse_region_set("R1,R2") == frozenset({"R1", "R2"})
    pat = pd.DataFrame({"patient_id": ["a", "b", "c"], "exons_tested": ["R1,R2", np.nan, None]})
    assert region_sets(pat).tolist() == [frozenset({"R1", "R2"}), None, None]
    assert panel_eligible(pat, ["R1", "R2"]).tolist() == [True, False, False]
    assert record_state(pat, "R1").tolist() == ["tested", "undetermined", "undetermined"]


# ------------------------------------------------------- None / NaN / pd.NA ----
MISSING = [None, float("nan"), np.nan, np.float32("nan"), np.float64("nan"), pd.NA]


@pytest.mark.parametrize("missing", MISSING, ids=repr)
def test_every_kind_of_missing_reference_value_is_missing_not_zero(missing):
    assert C.population_flag(missing, C.BA1_THRESHOLD) == C.NOT_EVALUABLE
    assert C.population_flag(missing, C.BS1_THRESHOLD) == C.NOT_EVALUABLE
    ratio, status = C.enrichment(0.01, missing)
    assert math.isnan(ratio) and status == "not_evaluable_reference_missing"
    ratio, status = C.enrichment(missing, 0.01)
    assert math.isnan(ratio) and status == "not_evaluable_cohort_frequency_undefined"
    group, freq, complete = C.highest_group({"a": missing, "b": 0.2}, ["a", "b"])
    assert (group, freq, complete) == ("b", 0.2, False)
    assert C.highest_group({"a": missing}, ["a"])[0] is None


def test_zero_is_a_value_and_is_not_missing():
    assert C.population_flag(0.0, C.BS1_THRESHOLD) == "no"
    assert C.enrichment(0.01, 0.0)[1] == "undefined_reference_zero"
    assert C.enrichment(0.0, 0.01) == (0.0, "computed")
    assert C.population_flag(np.float32(0.2), C.BA1_THRESHOLD) == "yes"


def test_bh_fdr_leaves_missing_values_missing():
    from mefv.statistics import bh_fdr
    out = bh_fdr([0.01, np.nan, 0.04, None, 0.03])
    assert np.isnan(out[1]) and np.isnan(out[3]) and np.isfinite(out[[0, 2, 4]]).all()


# --------------------------------------------------- empty-data edge cases ----
def test_statistics_on_no_individuals_are_undefined_not_zero_and_say_not_run():
    r = ld_stats([], [])
    assert math.isnan(r["D"]) and math.isnan(r["Dprime"]) and math.isnan(r["r2"])
    assert r["converged"] is None             # EM was not run: not "failed to converge"
    b = bootstrap_dprime([], [], 50, np.random.default_rng(1))
    assert math.isnan(b["lo"]) and math.isnan(b["hi"])
    assert (b["n_requested"], b["n_valid"], b["n_dropped"]) == (50, 0, 50)


def test_a_monomorphic_locus_gives_nan_and_a_converged_flag():
    r = ld_stats([0, 0, 0, 0], [0, 1, 2, 1])
    assert math.isnan(r["Dprime"]) and r["converged"] is True


def test_em_reports_a_non_converged_fit():
    from mefv.statistics import em_haplotypes_info
    gA, gB = [1, 1, 0, 2, 1, 0, 1, 2] * 5, [1, 0, 1, 2, 1, 1, 0, 2] * 5
    _, ok = em_haplotypes_info(gA, gB, iters=1)
    assert ok is False
    _, ok = em_haplotypes_info(gA, gB)
    assert ok is True


# ------------------------------------------------------- scanner read errors ----
def test_an_unreadable_file_is_a_scan_failure_not_a_skipped_file(tmp_path, monkeypatch):
    f = tmp_path / "a.txt"
    f.write_text("x")
    real = type(f).read_bytes

    def boom(self):
        if self.name == "a.txt":
            raise PermissionError("denied: /private/name")
        return real(self)

    monkeypatch.setattr(type(f), "read_bytes", boom)
    with pytest.raises(privacy.ScanError) as e:
        privacy.path_items([f], tmp_path)
    assert "/private/name" not in str(e.value) and "a.txt" not in str(e.value)


# ----------------------------------------------- the release gate fails closed ----
@pytest.mark.parametrize("content", [
    None,                                  # no report at all
    "{not json",                           # malformed
    "",                                    # empty file
    "[]",                                  # not an object
    "null",
    "{}",                                  # no outputs_released key
    '{"outputs_released": null}',
    '{"outputs_released": "true"}',        # a string is not the JSON value true
    '{"outputs_released": 1}',             # nor is a number
    '{"outputs_released": "yes"}',
    '{"outputs_released": false}',
    '{"released": true}',                  # wrong key
])
def test_the_release_gate_refuses_unless_the_report_says_true(tmp_path, content):
    from mefv.utils import require_released
    pat = tmp_path / "patients.csv"
    pat.write_text("x\n")
    if content is not None:
        (tmp_path / "validation_report.json").write_text(content)
    with pytest.raises(DataIntegrityError, match="not released"):
        require_released(pat)


def test_the_release_gate_accepts_only_an_explicit_true(tmp_path):
    from mefv.utils import require_released
    pat = tmp_path / "patients.csv"
    pat.write_text("x\n")
    (tmp_path / "validation_report.json").write_text('{"outputs_released": true}')
    require_released(pat)


@pytest.mark.parametrize("content", ["{not json", "{}", '{"outputs_released": false}', "[]"])
def test_allow_unreported_skips_only_a_missing_report(tmp_path, content):
    from mefv.utils import require_released
    pat = tmp_path / "patients.csv"
    pat.write_text("x\n")
    require_released(pat, allow_unreported=True)          # no report: allowed on request
    (tmp_path / "validation_report.json").write_text(content)
    with pytest.raises(DataIntegrityError):              # a report that exists is always judged
        require_released(pat, allow_unreported=True)


def test_load_refuses_tables_without_a_release_record(tmp_path):
    from mefv.utils import load
    (tmp_path / "p.csv").write_text("patient_id;year_requested;n_variants\nP1;2003;0\n")
    (tmp_path / "v.csv").write_text("patient_id;year_requested;variant\n")
    # the gate runs before the design is used, so a placeholder design is enough
    with pytest.raises(DataIntegrityError, match="not released"):
        load(tmp_path / "p.csv", tmp_path / "v.csv", object())
