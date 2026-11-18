"""Tests for the input-validation rules (src/mefv/validation.py).

All data are invented in this file. Each rule is exercised by injecting exactly
one defect into a small, clean set of tables and checking the rule identifier,
the severity and the counts that come back. A second group checks that nothing
the report, its text rendering or the exception messages contain is a value
from the data.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from mefv import DataIntegrityError, allele_dose, genotype_vector
from mefv import validation as V
from mefv.config import parse_study_design
from mefv.utils import load

PREFIX, WIDTH = "TST_", 3
DESIGN = parse_study_design(
    {
        "study_design": {
            "periods": {
                "primary": {"label": "early", "first_year": 2000, "last_year": 2010},
                "sensitivity": {"label": "late", "first_year": 2011, "last_year": 2020},
            },
            "panel_regions": ["R1", "R2"],
            "assays": {"panel": "ASSAY_T", "extended": "ASSAY_W"},
            "assay_inference": None,
        }
    }
)
FIXES = {"SYN_TYPO": "SYN_A"}


def clean() -> dict[str, pd.DataFrame]:
    """Four individuals; one has an extended-only region (R3) that was tested."""
    method = pd.DataFrame(
        {
            "Patient_ID": ["tst_1", "TST_002", "TST_003", "TST_004"],
            "YEAR": ["2003", "2005", "2012", "2015"],
            "sex": ["male", "female", "male", "female"],
            "Method": ["ASSAY_T", "ASSAY_T", "ASSAY_T", "ASSAY_W"],
        }
    )
    patients = pd.DataFrame(
        {
            "patient_id": ["TST_001", "TST_002", "TST_003", "TST_004"],
            "year_requested": ["2003", "2005", "2012", "2015"],
            "sex": ["male", "female", "male", "female"],
            "n_variants": ["2", "1", "0", "1"],
            "variant_detected": ["yes", "yes", "no", "yes"],
            "exons_tested": ["R1,R2", "R1,R2", "R1,R2", "R1,R2,R3"],
            "coverage_source": ["synthetic"] * 4,
        }
    )
    variants = pd.DataFrame(
        {
            "patient_id": ["TST_001", "TST_001", "TST_002", "TST_004"],
            "variant": ["SYN_A", "SYN_B", "SYN_A", "SYN_E"],
            "exon": ["R1", "R2", "R1", "R3"],
            "zygosity": ["Het", "Hom", "Het", "Het"],
            "year_requested": ["2003", "2003", "2005", "2015"],
        }
    )
    return {"method": method, "patients": patients, "variants": variants}


def run(t: dict[str, pd.DataFrame]) -> V.ValidationReport:
    return V.validate_raw_inputs(
        t["method"],
        t["variants"],
        t["patients"],
        design=DESIGN,
        id_prefix=PREFIX,
        id_width=WIDTH,
        nomenclature_fixes=FIXES,
    )


def row(df: pd.DataFrame, **kw) -> pd.DataFrame:
    """Append one row, filling unspecified columns from the first row."""
    new = df.iloc[[0]].copy()
    for k, v in kw.items():
        new[k] = v
    return pd.concat([df, new], ignore_index=True)


# ------------------------------------------------------------------ clean ----
def test_clean_tables_have_no_actionable_finding():
    rep = run(clean())
    assert rep.status == "passed"
    assert {f.rule_id for f in rep.findings if f.rule.severity != V.INFO} == set()
    rep.enforce()  # does not raise


# ------------------------------------------- one defect per rule: (rule, n) ---
def _m_empty_patient_row(t):
    t["patients"] = pd.concat(
        [t["patients"], pd.DataFrame([[""] * 7], columns=t["patients"].columns)], ignore_index=True
    )


def _m_blank_pid_with_content(t):
    t["patients"] = row(t["patients"], patient_id="")


def _m_noncanonical_patient_id(t):
    t["patients"].loc[0, "patient_id"] = "tst_1"


def _m_id_collision(t):
    # a second raw spelling of an existing canonical identifier
    t["patients"] = row(t["patients"], patient_id="TST_01")


def _m_dup_patient_exact(t):
    t["patients"] = row(t["patients"])


def _m_dup_patient_conflict(t):
    t["patients"] = row(t["patients"], year_requested="2004")


def _m_year_text(t):
    t["patients"].loc[1, "year_requested"] = "20x5"


def _m_year_outside(t):
    t["patients"].loc[1, "year_requested"] = "1999"


def _m_sex_case(t):
    t["patients"].loc[1, "sex"] = "Female"


def _m_detected_blank(t):
    t["patients"].loc[1, "variant_detected"] = ""


def _m_nvar_text(t):
    t["patients"].loc[1, "n_variants"] = "one"


def _m_exons_blank(t):
    t["patients"].loc[2, "exons_tested"] = ""


def _m_exons_separator(t):
    t["patients"].loc[2, "exons_tested"] = "R1; R2"


def _m_exons_chars(t):
    t["patients"].loc[2, "exons_tested"] = "R1,R2,R#"


def _m_exons_unknown(t):
    t["patients"].loc[2, "exons_tested"] = "R1,R2,RZ9"


def _m_exons_no_panel(t):
    t["patients"].loc[2, "exons_tested"] = "R1"


def _m_var_empty_row(t):
    t["variants"] = pd.concat(
        [t["variants"], pd.DataFrame([[""] * 5], columns=t["variants"].columns)], ignore_index=True
    )


def _m_var_blank_pid(t):
    t["variants"] = row(t["variants"], patient_id="")


def _m_var_noncanonical_pid(t):
    t["variants"].loc[0, "patient_id"] = "tst_1"


def _m_var_blank_label(t):
    t["variants"].loc[0, "variant"] = " "


def _m_var_blank_exon(t):
    t["variants"].loc[0, "exon"] = ""


def _m_zyg_case(t):
    t["variants"].loc[0, "zygosity"] = "het"


def _m_zyg_unknown(t):
    t["variants"].loc[0, "zygosity"] = "Homo"


def _m_zyg_blank(t):
    t["variants"].loc[0, "zygosity"] = ""


def _m_var_exact_dup(t):
    t["variants"] = row(t["variants"])


def _m_var_same_genotype_other_column(t):
    t["variants"] = row(t["variants"], year_requested="2004")


def _m_var_conflict(t):
    t["variants"] = row(t["variants"], zygosity="Hom")


def _m_label_case(t):
    t["variants"].loc[2, "variant"] = "syn_a"


def _m_nomenclature_fix(t):
    t["variants"].loc[2, "variant"] = "SYN_TYPO"


def _m_fix_creates_duplicate(t):
    # a mistyped label that the documented correction maps onto an existing one
    t["variants"] = row(t["variants"], variant="SYN_TYPO")


def _m_two_regions(t):
    t["variants"].loc[2, "exon"] = "R2"


def _m_var_year_text(t):
    t["variants"].loc[0, "year_requested"] = "x"


def _m_var_year_outside(t):
    t["variants"].loc[0, "year_requested"] = "2030"


def _m_method_unparsable(t):
    t["method"].loc[1, "Patient_ID"] = "BADID"


def _m_method_year(t):
    t["method"].loc[1, "YEAR"] = "n/a"


def _m_method_empty_row(t):
    t["method"] = pd.concat(
        [t["method"], pd.DataFrame([[""] * 4], columns=t["method"].columns)], ignore_index=True
    )


def _m_method_dup(t):
    t["method"] = row(t["method"])


def _m_orphan(t):
    t["variants"] = row(t["variants"], patient_id="TST_099")


def _m_method_only(t):
    # individual 3 exists only in the method table; one variant row belongs to it
    t["patients"] = t["patients"][t["patients"].patient_id != "TST_003"].reset_index(drop=True)
    t["variants"] = row(t["variants"], patient_id="TST_003")


def _m_no_method(t):
    t["method"] = t["method"][t["method"].Patient_ID != "TST_003"]


def _m_nvar_above(t):
    t["patients"].loc[1, "n_variants"] = "3"


def _m_nvar_below(t):
    t["patients"].loc[0, "n_variants"] = "1"


def _m_nvar_pairs(t):
    # duplicate of an existing row: declared 2 matches the 2 distinct pairs, not the 3 rows
    t["variants"] = row(t["variants"], patient_id="TST_001")
    t["variants"].loc[len(t["variants"]) - 1, "variant"] = "SYN_A"
    t["variants"].loc[len(t["variants"]) - 1, "exon"] = "R1"
    t["variants"].loc[len(t["variants"]) - 1, "zygosity"] = "Het"


def _m_flag_declared(t):
    t["patients"].loc[2, "variant_detected"] = "yes"


def _m_flag_rows(t):
    t["patients"].loc[2, "variant_detected"] = "yes"
    t["patients"].loc[2, "n_variants"] = "1"


def _m_outside_tested(t):
    t["patients"].loc[3, "exons_tested"] = "R1,R2"


def _m_variants_blank_exons_tested(t):
    t["patients"].loc[3, "exons_tested"] = ""


def _m_year_mismatch(t):
    t["method"].loc[1, "YEAR"] = "2006"


def _m_sex_mismatch(t):
    t["method"].loc[1, "sex"] = "male"


# rule id -> (mutator, records, individuals). ``None`` means "not asserted".
CASES = {
    "P1a": (_m_empty_patient_row, 1, 0),
    "P1b": (_m_blank_pid_with_content, 1, 0),
    "P2": (_m_noncanonical_patient_id, 1, 0),
    "P3": (_m_id_collision, 2, 1),
    "P4a": (_m_dup_patient_exact, 2, 1),
    "P4b": (_m_dup_patient_conflict, 2, 1),
    "P5a": (_m_year_text, 1, 0),
    "P5b": (_m_year_outside, 1, 0),
    "P6": (_m_sex_case, 1, 0),
    "P7": (_m_detected_blank, 1, 0),
    "P8": (_m_nvar_text, 1, 0),
    "E1": (_m_exons_blank, 1, 0),
    "E2": (_m_exons_separator, 1, 0),
    "E3": (_m_exons_chars, 1, 0),
    "E4a": (_m_exons_unknown, 1, 0),
    "E5": (_m_exons_no_panel, 1, 0),
    "V1a": (_m_var_empty_row, 1, 0),
    "V1b": (_m_var_blank_pid, 1, 0),
    "V10": (_m_var_noncanonical_pid, 1, 0),
    "V2": (_m_var_blank_label, 1, 0),
    "V3": (_m_var_blank_exon, 1, 0),
    "V4": (_m_zyg_case, 1, 0),
    "V5": (_m_var_exact_dup, 2, 1),
    "V6a": (_m_var_same_genotype_other_column, 2, 1),
    "V6b": (_m_var_conflict, 2, 1),
    "V7": (_m_label_case, 2, 0),
    "V8": (_m_nomenclature_fix, 1, 0),
    "V9": (_m_two_regions, 2, 0),
    "V11": (_m_var_year_text, 1, 0),
    "V12": (_m_var_year_outside, 1, 0),
    "M1": (_m_method_unparsable, 1, 0),
    "M2": (_m_method_year, 1, 0),
    "M3": (_m_method_empty_row, 1, 0),
    "X4": (_m_method_dup, 2, 1),
    "X1": (_m_orphan, 1, 1),
    "X2": (_m_method_only, 1, 1),
    "X2b": (_m_method_only, 1, 1),
    "X3": (_m_no_method, 1, 1),
    "X5a": (_m_nvar_above, 1, 1),
    "X5b": (_m_nvar_below, 1, 1),
    "X5c": (_m_nvar_pairs, 1, 0),
    "X6": (_m_flag_declared, 1, 0),
    "X7": (_m_flag_rows, 1, 0),
    "X8": (_m_outside_tested, 1, 1),
    "X9": (_m_variants_blank_exons_tested, 1, 1),
    "X10": (_m_year_mismatch, 1, 1),
    "X11": (_m_sex_mismatch, 1, 1),
}


@pytest.mark.parametrize("rule_id", sorted(CASES))
def test_rule_detects_its_defect_with_exact_counts(rule_id):
    mutate, records, individuals = CASES[rule_id]
    t = clean()
    mutate(t)
    found = {f.rule_id: f for f in run(t).findings}
    assert rule_id in found, f"{rule_id} not reported; got {sorted(found)}"
    assert (found[rule_id].records, found[rule_id].individuals) == (records, individuals)


def test_e6_reports_the_distinct_region_sets():
    f = {f.rule_id: f for f in run(clean()).findings}["E6"]
    assert f.detail == {
        "distinct_region_sets": 2,
        "individuals_with_2_regions": 3,
        "individuals_with_3_regions": 1,
    }


def test_case_variants_are_classified_by_shape_not_echoed():
    t = clean()
    t["variants"].loc[0, "zygosity"] = "het"
    t["variants"].loc[1, "zygosity"] = ""
    t["variants"].loc[2, "zygosity"] = "Homo"
    f = {f.rule_id: f for f in run(t).findings}["V4"]
    assert f.detail == {
        "case_or_whitespace_variant_of_valid": 1,
        "blank": 1,
        "unrecognised": 1,
    }


def test_nomenclature_fix_that_creates_a_duplicate_is_a_blocking_finding():
    t = clean()
    _m_fix_creates_duplicate(t)
    found = {f.rule_id: f for f in run(t).findings}
    assert "V6a" in found and "V5" not in found  # rows differ before the correction


def test_missing_column_is_reported_without_evaluating_other_rules():
    t = clean()
    t["variants"] = t["variants"].drop(columns=["zygosity"])
    rep = run(t)
    assert [f.rule_id for f in rep.findings] == ["S1"]
    assert rep.status == "failed"


def test_built_table_checks_replace_the_asserts():
    pat = pd.DataFrame({"assay": ["ASSAY_T", "not_recorded", None, "OTHER"]})
    rep = V.validate_built_tables(pat, DESIGN.assays, 5, V.ValidationReport())
    found = {f.rule_id: f for f in rep.findings}
    assert (found["A1"].records, found["A2"].records) == (2, 1)
    assert found["C1"].detail == {"expected": 5, "found": 4}
    ok = V.validate_built_tables(
        pd.DataFrame({"assay": ["ASSAY_T"] * 3}), DESIGN.assays, 3, V.ValidationReport()
    )
    assert ok.findings == []


def test_every_rule_in_the_catalogue_is_covered_by_a_test():
    tested = set(CASES) | {"E6", "S1", "A1", "A2", "C1"}
    assert tested == set(V.RULES), sorted(set(V.RULES) ^ tested)


# --------------------------------------------------------------- severity ----
@pytest.mark.parametrize(
    ("rule_id", "severity"),
    [
        ("P4a", V.FATAL), ("P4b", V.FATAL), ("V1b", V.FATAL), ("V4", V.FATAL), ("X1", V.FATAL),
        ("V5", V.BLOCKING), ("V6a", V.BLOCKING), ("V6b", V.BLOCKING), ("V7", V.BLOCKING),
        ("X5a", V.WARNING), ("X5b", V.WARNING), ("X8", V.WARNING), ("E1", V.WARNING),
        ("V8", V.INFO), ("E6", V.INFO), ("X2b", V.INFO),
    ],
)
def test_severity_policy(rule_id, severity):
    assert V.RULES[rule_id].severity == severity


def test_only_blocking_rules_can_be_acknowledged():
    assert set(V.ACKNOWLEDGEABLE) == {"V5", "V6a", "V6b", "V7"}
    assert {r.id for r in V.RULES.values() if r.severity == "blocking"} == set(V.ACKNOWLEDGEABLE)


def test_declared_and_observed_counts_are_both_reported_and_neither_wins():
    t = clean()
    _m_nvar_above(t)
    f = {f.rule_id: f for f in run(t).findings}["X5a"]
    assert f.detail == {"declared_total": 3, "observed_rows_total": 1}
    run(t).enforce()  # a warning: reported, does not stop the run


# ------------------------------------------------------------ enforcement ----
def test_fatal_findings_raise_with_rule_ids_and_counts():
    t = clean()
    _m_dup_patient_exact(t)
    with pytest.raises(DataIntegrityError, match=r"P4a: <5 record"):
        run(t).enforce()


def _rec(n):
    """The record count as the report shows it: small counts are suppressed."""
    small = isinstance(n, int) and not isinstance(n, bool) and 0 < n < 5
    return V.SUPPRESSED if small else n


def _bind(rep, tag="cd"):
    """Give the report the digests of "the run being read now" and return a binding made for
    exactly them."""
    rep.input_sha256 = tag * 32
    rep.config_sha256 = "cf" * 32
    return V.DispositionBinding(tag * 32, V.ruleset_sha256(), "cf" * 32)


def _disp(rule, _records, **kw):
    d = {"rule": rule, "records": _rec(_records), "decision": "use_as_recorded",
         "decided_by": "owner", "decided_on": "2026-01-02", "basis": "reviewed the export"}
    d.update(kw)
    return d


BLOCKING_CASES = [("V5", _m_var_exact_dup), ("V6b", _m_var_conflict), ("V7", _m_label_case)]


@pytest.mark.parametrize(("rule", "mutate"), BLOCKING_CASES)
def test_unacknowledged_blocking_finding_withholds_outputs(rule, mutate):
    t = clean()
    mutate(t)
    rep = run(t)
    with pytest.raises(DataIntegrityError, match="withheld"):
        rep.enforce()
    assert rep.status == "blocked" and rep.releases_outputs is False
    assert rule in [f.rule_id for f in rep.blocking]


@pytest.mark.parametrize(("rule", "mutate"), BLOCKING_CASES)
def test_acknowledging_does_not_release_the_run(rule, mutate):
    """Acknowledgement marks a finding as seen. It must not allow final outputs."""
    t = clean()
    mutate(t)
    rep = run(t)
    with pytest.raises(DataIntegrityError, match="not released|withheld"):
        rep.enforce({rule})
    assert rep.acknowledged == (rule,)
    assert rep.status in ("acknowledged_unresolved", "blocked")
    assert rep.releases_outputs is False
    assert rule in rep.acknowledged_unresolved
    d = rep.to_dict()
    assert d["outputs_released"] is False and d["owner_dispositions"] == []
    assert rule in d["acknowledged_unresolved_rules"]


def test_acknowledging_every_blocking_rule_present_is_still_not_a_release():
    t = clean()
    _m_var_exact_dup(t)        # V5 and V6a
    rep = run(t)
    with pytest.raises(DataIntegrityError, match="Acknowledging does not release"):
        rep.enforce({"V5", "V6a"})
    assert rep.status == "acknowledged_unresolved"


@pytest.mark.parametrize(("rule", "mutate"), BLOCKING_CASES)
def test_an_owner_disposition_for_the_exact_count_releases_and_is_recorded(rule, mutate):
    t = clean()
    mutate(t)
    rep = run(t)
    n = next(f.records for f in rep.findings if f.rule_id == rule)
    others = [f for f in rep.findings if f.rule.severity == "blocking" and f.rule_id != rule]
    disp = [_disp(rule, n)] + [_disp(f.rule_id, f.records) for f in others]
    rep.enforce(set(), disp, "ab" * 32, _bind(rep))   # no acknowledgement needed, none sufficient
    assert rep.releases_outputs and rep.status == "owner_disposition_recorded"
    d = rep.to_dict()
    assert d["outputs_released"] is True
    assert {"rule": rule, "records": _rec(n), "decision": "use_as_recorded"} in d[
        "owner_dispositions"]
    assert d["owner_disposition_file_sha256"] == "ab" * 32
    # the findings are still in the report: the conflict is recorded, not removed
    assert rule in [f["rule"] for f in d["findings"]]
    assert "owner disposition recorded" in rep.render()


def test_a_disposition_for_a_different_count_does_not_release():
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    n = next(f.records for f in rep.findings if f.rule_id == "V6b")
    b = _bind(rep)
    with pytest.raises(DataIntegrityError, match="different record count"):
        rep.enforce(set(), [_disp("V6b", n + 10)], "", b)   # a count the data does not have
    assert rep.releases_outputs is False


def test_a_suppressed_count_covers_only_a_finding_that_is_really_small():
    """'<5' stands for a small count; it can never cover a finding of five or more records."""
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    f = next(f for f in rep.findings if f.rule_id == "V6b")
    b = _bind(rep)
    f.records = 7                                           # as if the data held more records
    with pytest.raises(DataIntegrityError, match="different record count"):
        rep.enforce(set(), [_disp("V6b", V.SUPPRESSED)], "", b)
    rep2 = run(t)
    rep2.enforce(set(), [_disp("V6b", V.SUPPRESSED)] + [
        _disp(g.rule_id, V.SUPPRESSED) for g in rep2.findings
        if g.rule.severity == "blocking" and g.rule_id != "V6b"], "", _bind(rep2))


def test_a_disposition_for_another_rule_does_not_release():
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    with pytest.raises(DataIntegrityError, match="withheld"):
        rep.enforce(set(), [_disp("V5", 1)], "", _bind(rep))


def test_a_disposition_does_not_cover_a_second_blocking_rule():
    t = clean()
    _m_var_exact_dup(t)        # V5 and V6a
    rep = run(t)
    n5 = next(f.records for f in rep.findings if f.rule_id == "V5")
    with pytest.raises(DataIntegrityError, match="V6a"):
        rep.enforce(set(), [_disp("V5", n5)], "", _bind(rep))


def test_a_disposition_cannot_release_a_fatal_finding():
    t = clean()
    _m_zyg_unknown(t)
    rep = run(t)
    with pytest.raises(DataIntegrityError, match="input validation failed"):
        rep.enforce(set(), [_disp("V6b", 1)], "", _bind(rep))


def test_acknowledging_one_blocking_rule_does_not_waive_another():
    t = clean()
    _m_var_exact_dup(t)  # V5 and V6a
    with pytest.raises(DataIntegrityError, match="V6a"):
        run(t).enforce({"V5"})


# ------------------------------------------------- disposition file ---------
BOUND = "ef" * 32


def _write(tmp_path, entries, bound=BOUND, **extra):
    p = tmp_path / "disposition.json"
    doc = {"ruleset_sha256": V.ruleset_sha256(), "config_sha256": "cf" * 32,
           "dispositions": entries, **extra}
    if bound is not None:
        doc["input_sha256"] = bound
    p.write_text(json.dumps(doc))
    return str(p)


def test_a_valid_disposition_file_loads_with_its_digest(tmp_path):
    import hashlib
    from pathlib import Path

    path = _write(tmp_path, [_disp("V6b", 2)])
    entries, digest, bound = V.load_dispositions(path)
    assert entries[0]["rule"] == "V6b" and bound.input_sha256 == BOUND
    assert bound.ruleset_sha256 == V.ruleset_sha256() and bound.config_sha256 == "cf" * 32
    assert digest == hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.mark.parametrize(
    "bad",
    [
        _disp("X8", 1),                              # not a blocking rule
        _disp("V6b", 0), _disp("V6b", -1), _disp("V6b", True), _disp("V6b", "2"),
        _disp("V6b", 2, decision="delete_rows"),     # the software never deletes or merges
        _disp("V6b", 2, decision="merge"),
        _disp("V6b", 2, decided_by=""), _disp("V6b", 2, basis="  "),
        _disp("V6b", 2, decided_on="2 Jan 2026"),
        _disp("V6b", 2, decided_on="2026-02-30"),    # not a real calendar date
        _disp("V6b", 2, decided_on="2026-13-01"), _disp("V6b", 2, decided_on="20260102"),
        _disp("V6b", 2, decided_on="2026-W01-1"), _disp("V6b", 2, decided_on=20260102),
        _disp("V6b", 2, extra_field="x"),            # unknown field
        _disp("V6b", 2, records=3), _disp("V6b", 2, records="<6"),  # only "<5" or an integer >= 5
    ],
)
def test_malformed_dispositions_are_refused(tmp_path, bad):
    with pytest.raises(DataIntegrityError, match="owner disposition entry 1"):
        V.load_dispositions(_write(tmp_path, [bad]))


@pytest.mark.parametrize("bound", [None, "", "abc", "G" * 64, 7])
def test_a_disposition_must_name_the_input_it_was_made_for(tmp_path, bound):
    with pytest.raises(DataIntegrityError, match="input_sha256"):
        V.load_dispositions(_write(tmp_path, [_disp("V6b", 2)], bound=bound))


def test_a_disposition_made_for_other_input_files_does_not_release_even_with_equal_counts():
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    n = next(f.records for f in rep.findings if f.rule_id == "V6b")
    b = _bind(rep, "11")                                  # the files being read now
    other = V.DispositionBinding("22" * 32, b.ruleset_sha256, b.config_sha256)
    with pytest.raises(DataIntegrityError, match="different input files"):
        rep.enforce(set(), [_disp("V6b", n)], "", other)   # made for other files
    assert rep.releases_outputs is False and rep.status == "blocked"


def test_a_report_without_an_input_digest_never_matches_a_disposition():
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    n = next(f.records for f in rep.findings if f.rule_id == "V6b")
    with pytest.raises(DataIntegrityError):
        rep.enforce(set(), [_disp("V6b", n)], "", "")


@pytest.mark.parametrize("what", ["config", "ruleset"])
def test_a_disposition_made_under_other_settings_or_rules_does_not_release(what):
    """Same data, same counts, but the study settings or the validation rules differ."""
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    n = next(f.records for f in rep.findings if f.rule_id == "V6b")
    b = _bind(rep)
    other = (b._replace(config_sha256="99" * 32) if what == "config"
             else b._replace(ruleset_sha256="88" * 32))
    with pytest.raises(DataIntegrityError, match="different input files, settings or validation"):
        rep.enforce(set(), [_disp("V6b", n)], "", other)
    assert rep.releases_outputs is False and rep.status == "blocked"


def test_a_report_without_a_config_digest_never_matches_a_disposition():
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    n = next(f.records for f in rep.findings if f.rule_id == "V6b")
    b = _bind(rep)
    rep.config_sha256 = ""
    with pytest.raises(DataIntegrityError):
        rep.enforce(set(), [_disp("V6b", n)], "", b)


def test_the_rule_catalogue_digest_changes_when_a_rule_changes(monkeypatch):
    before = V.ruleset_sha256()
    rules = dict(V.RULES)
    r = rules["V6b"]
    rules["V6b"] = V.Rule(r.id, "warning", r.table, r.statement)   # a severity change
    monkeypatch.setattr(V, "RULES", rules)
    assert V.ruleset_sha256() != before


def test_the_config_digest_ignores_formatting_but_not_values():
    a = {"id_prefix": "X", "design": {"b": 1, "a": [1, 2]}}
    b = {"design": {"a": [1, 2], "b": 1}, "id_prefix": "X"}          # same values, other order
    assert V.digest_config(a) == V.digest_config(b)
    assert V.digest_config(a) != V.digest_config({**a, "id_prefix": "Y"})
    assert V.digest_config(a) != V.digest_config({**a, "design": {"b": 1, "a": [1, 3]}})


def test_duplicate_json_keys_and_non_finite_numbers_are_refused(tmp_path):
    good = json.dumps({"ruleset_sha256": V.ruleset_sha256(), "config_sha256": "cf" * 32,
                       "input_sha256": BOUND, "dispositions": [_disp("V6b", 2)]})
    p = tmp_path / "d.json"
    for text, match in (
        (good[:-1] + ', "input_sha256": "' + "aa" * 32 + '"}', "key is repeated"),
        (good.replace('"rule": "V6b"', '"rule": "V6b", "rule": "V5"'), "key is repeated"),
        (good.replace('"records": "<5"', '"records": NaN'), "not valid JSON"),
        (good.replace('"records": "<5"', '"records": Infinity'), "not valid JSON"),
    ):
        p.write_text(text)
        with pytest.raises(DataIntegrityError, match=match):
            V.load_dispositions(str(p))
    p.write_text(good)
    V.load_dispositions(str(p))                                       # the good one loads
    p.write_bytes(b"\xff\xfe not utf-8")
    with pytest.raises(DataIntegrityError, match="UTF-8"):
        V.load_dispositions(str(p))
    p.write_text(json.dumps({"input_sha256": BOUND, "bogus": 1}))
    with pytest.raises(DataIntegrityError):
        V.load_dispositions(str(p))
    p.write_text(good.replace("{", '{"surprise": 1,', 1))
    with pytest.raises(DataIntegrityError, match="unknown top-level key"):
        V.load_dispositions(str(p))


def test_input_digests_change_with_any_byte_and_name_the_file(tmp_path):
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    a.write_text("x\n"), b.write_text("y\n")
    each, comb = V.digest_inputs({"a": a, "b": b})
    b.write_text("y \n")
    each2, comb2 = V.digest_inputs({"a": a, "b": b})
    assert each["a"] == each2["a"] and each["b"] != each2["b"] and comb != comb2
    assert set(each) == {"a", "b"} and len(comb) == 64


def test_missing_field_duplicate_rule_and_bad_shape_are_refused(tmp_path):
    d = _disp("V6b", 2)
    del d["basis"]
    for entries in ([d], [_disp("V6b", 2), _disp("V6b", 3)], [], ["x"]):
        with pytest.raises(DataIntegrityError):
            V.load_dispositions(_write(tmp_path, entries))
    p = tmp_path / "bad.json"
    p.write_text("{not json")
    with pytest.raises(DataIntegrityError, match="not valid JSON"):
        V.load_dispositions(str(p))


def test_disposition_errors_and_reports_never_contain_the_owner_text(tmp_path):
    secret_name, secret_basis = "LEAKNAME", "LEAKBASIS"
    bad = _disp("V6b", 2, decided_by=secret_name, basis=secret_basis, decision="LEAKDECISION")
    with pytest.raises(DataIntegrityError) as e:
        V.load_dispositions(_write(tmp_path, [bad]))
    assert not any(x in str(e.value) for x in (secret_name, secret_basis, "LEAKDECISION"))
    t = clean()
    _m_var_conflict(t)
    rep = run(t)
    n = next(f.records for f in rep.findings if f.rule_id == "V6b")
    entries, digest, bound = V.load_dispositions(
        _write(tmp_path, [_disp("V6b", n, decided_by=secret_name, basis=secret_basis)]))
    rep.input_sha256 = bound.input_sha256
    rep.config_sha256 = bound.config_sha256
    rep.enforce(set(), entries, digest, bound)
    blob = json.dumps(rep.to_dict()) + rep.render()
    assert secret_name not in blob and secret_basis not in blob


# ------------------------------------------- later stages refuse the tables ---
def test_analysis_tables_with_an_unreleased_validation_report_are_refused(tmp_path):
    for released, ok in ((False, False), (True, True)):
        d = tmp_path / str(released)
        d.mkdir()
        (d / "patients.csv").write_text("patient_id;year_requested;n_variants\nP1;2003;0\n")
        (d / "variants_long.csv").write_text("patient_id;year_requested;variant\n")
        (d / "validation_report.json").write_text(json.dumps({"outputs_released": released}))
        if ok:
            load(d / "patients.csv", d / "variants_long.csv", DESIGN)
        else:
            with pytest.raises(DataIntegrityError, match="not released"):
                load(d / "patients.csv", d / "variants_long.csv", DESIGN)


def test_a_fatal_rule_cannot_be_acknowledged():
    t = clean()
    _m_zyg_unknown(t)
    with pytest.raises(DataIntegrityError, match="cannot acknowledge"):
        run(t).enforce({"V4"})
    with pytest.raises(DataIntegrityError, match="input validation failed"):
        run(t).enforce({"V5"})


def test_fatal_takes_precedence_over_acknowledged_blocking():
    t = clean()
    _m_var_exact_dup(t)
    _m_zyg_unknown(t)
    with pytest.raises(DataIntegrityError, match="V4"):
        run(t).enforce({"V5", "V6a"})


# ------------------------------------------------- nothing is changed -------
def test_validation_does_not_modify_its_inputs():
    t = clean()
    _m_var_conflict(t)
    _m_label_case(t)
    before = {k: v.copy() for k, v in t.items()}
    run(t)
    for k in t:
        pd.testing.assert_frame_equal(t[k], before[k])


# ---------------------------------------------------------------- privacy ----
SENTINELS = ("LEAKID", "LEAKVAR", "LEAKREG", "LEAKZYG", "LEAK2025")


def _leaky() -> dict[str, pd.DataFrame]:
    """Defects whose values are recognisable strings; none may reach any output."""
    t = clean()
    t["patients"] = row(t["patients"], patient_id="LEAKID_1", year_requested="LEAK2025")
    t["patients"].loc[2, "exons_tested"] = "R1;LEAKREG"
    t["variants"] = row(t["variants"], variant="LEAKVAR", zygosity="LEAKZYG", patient_id="LEAKID_2")
    t["variants"] = row(t["variants"], variant="LEAKVAR", exon="LEAKREG")
    t["method"].loc[1, "Patient_ID"] = "LEAKID_3"
    t["variants"] = row(t["variants"], patient_id="TST_098")  # a canonical orphan
    return t


def test_report_rendering_and_error_messages_contain_no_data_values():
    t = _leaky()
    rep = run(t)
    assert rep.status == "failed"
    texts = [rep.render(), json.dumps(rep.to_dict())]
    try:
        rep.enforce()
    except DataIntegrityError as exc:
        texts.append(str(exc))
    blob = "\n".join(texts)
    for s in (*SENTINELS, "TST_00", "SYN_", "R3"):
        assert s not in blob, f"value leaked into validation output: {s}"


def test_leaky_fixture_actually_triggers_many_rules():
    ids = {f.rule_id for f in run(_leaky()).findings}
    assert {"P2", "P5a", "E2", "V4", "V10", "M1", "X1"} <= ids


def test_to_dict_is_deterministic():
    a = json.dumps(run(_leaky()).to_dict(), sort_keys=True)
    b = json.dumps(run(_leaky()).to_dict(), sort_keys=True)
    assert a == b


# ---------------------------------------------------------------- helpers ----
@pytest.mark.parametrize(
    ("raw", "expected"),
    [("tst_1", "TST_001"), (" TST_012 ", "TST_012"), ("TST_0001", "TST_001"),
     ("TST_", None), ("BAD", None), ("", None), ("TST_1x", None), ("TST_²", None)],
)
def test_canonical_id(raw, expected):
    assert V.canonical_id(raw, PREFIX, WIDTH) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [("R1,R2", {"R1", "R2"}), ("R1", {"R1"}), ("", None), ("  ", None), ("R1, R2", None),
     ("R1;R2", None), ("R1,,R2", None), ("R1,R2,", None), ("R1,R#", None), ("R1|R2", None)],
)
def test_parse_region_set_never_reads_malformed_as_tested(value, expected):
    got = V.parse_region_set(value)
    assert (None if expected is None else set(got)) == expected


def test_parse_uint_column_reports_counts_not_values():
    s = pd.Series(["2005", " 2006 ", "20x7", "-1", "", "3.0"])
    with pytest.raises(DataIntegrityError) as exc:
        V.parse_uint_column(s, "patients.year_requested")
    assert "<5 entr" in str(exc.value)               # four entries: below the limit, not shown
    assert "20x7" not in str(exc.value)
    big = pd.Series(["x"] * 6 + ["2005"])
    with pytest.raises(DataIntegrityError, match="6 entr"):
        V.parse_uint_column(big, "patients.year_requested")
    assert list(V.parse_uint_column(pd.Series(["2005", " 2006 "]), "x")) == [2005, 2006]


def test_shape_of():
    valid = ("Het", "Hom")
    assert [V.shape_of(x, valid) for x in ("Het", " het", "", "Hemi")] == [
        "valid", "case_or_whitespace_variant_of_valid", "blank", "unrecognised"]


# ------------------------------------------- downstream genotype handling ----
def test_unrecognised_zygosity_is_rejected_not_counted_as_one_allele():
    var = pd.DataFrame({"patient_id": ["P1"], "variant": ["SYN_P01"], "zygosity": ["LEAKZYG"]})
    with pytest.raises(DataIntegrityError) as exc:
        genotype_vector(var, "SYN_P01", ["P1"])
    assert "LEAKZYG" not in str(exc.value) and "<5 entr" in str(exc.value)


def test_allele_dose_maps_exactly_het_and_hom():
    out = allele_dose(pd.Series(["Het", "Hom", "Het"]))
    assert list(out) == [1, 2, 1] and out.dtype.kind == "i"
    for bad in ("het", "", "HOM"):
        with pytest.raises(DataIntegrityError):
            allele_dose(pd.Series(["Het", bad]))


def test_load_rejects_a_non_integer_year_without_echoing_it(tmp_path):
    pat = tmp_path / "p.csv"
    var = tmp_path / "v.csv"
    pat.write_text(
        "patient_id;year_requested;sex;n_variants\nTST_001;20LEAK;male;1\n", encoding="utf-8"
    )
    var.write_text("patient_id;year_requested\nTST_001;2003\n", encoding="utf-8")
    with pytest.raises(DataIntegrityError) as exc:
        load(pat, var, DESIGN, allow_unreported=True)
    assert "20LEAK" not in str(exc.value)


def test_documented_rules_match_the_catalogue():
    """docs/validation.md lists every rule with its severity, and no other."""
    import re
    from pathlib import Path

    doc = (Path(__file__).resolve().parents[1] / "docs" / "validation.md").read_text()
    listed = dict(re.findall(r"^\| `([A-Z]\d+[a-c]?)` \| (\w+) \|", doc, flags=re.M))
    assert listed == {r.id: r.severity for r in V.RULES.values()}
