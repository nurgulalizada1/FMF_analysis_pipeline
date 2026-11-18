"""Eligibility-aware denominators, with hand-calculated expected values.

Cohort used throughout (10 individuals, invented):
    P01-P05  exons_tested "R1,R2"   (R1 and R2 tested)
    P06-P07  exons_tested "R1"      (R1 tested, R2 not)
    P08      exons_tested ""        (testing record missing -> undetermined)
    P09      exons_tested "R1,,R2"  (malformed -> undetermined)
    P10      exons_tested "R1,R2"
so for R1: tested 8 (P01-P07, P10), not tested 0, undetermined 2;
   for R2: tested 6 (P01-P05, P10), not tested 2 (P06, P07), undetermined 2.
Request years: P01-P05 2003, P06-P10 2008; primary period 2001-2006.
"""
from __future__ import annotations

import math
from types import SimpleNamespace

import pandas as pd
import pytest

from mefv import eligibility as E


@pytest.fixture
def pat():
    tested = ["R1,R2"] * 5 + ["R1", "R1", "", "R1,,R2", "R1,R2"]
    return pd.DataFrame({
        "patient_id": [f"P{i:02d}" for i in range(1, 11)],
        "exons_tested": tested,
        "year_requested": ["2003"] * 5 + ["2008"] * 5,
    })


def by_id(pat, series):
    return dict(zip(pat.patient_id, series, strict=True))


def rows(*items):
    return pd.DataFrame(items, columns=["patient_id", "variant", "exon", "zygosity"])


@pytest.fixture
def var():
    return rows(
        ("P01", "VA", "R1", "Het"), ("P02", "VA", "R1", "Hom"), ("P08", "VA", "R1", "Het"),
        ("P01", "VB", "R2", "Het"), ("P06", "VB", "R2", "Het"),
    )


# --- who is eligible ---------------------------------------------------------------------
def test_denominator_counts_by_hand(pat):
    d, mask = E.denominator(pat, "R1")
    counts = (d.n_cohort, d.n_in_scope, d.n_tested, d.n_not_tested, d.n_undetermined)
    assert counts == (10, 10, 8, 0, 2)
    assert d.n_eligible == 8 and d.allele_denominator == 16
    assert mask.sum() == 8 and not mask[7] and not mask[8]
    d2, _ = E.denominator(pat, "R2")
    assert (d2.n_tested, d2.n_not_tested, d2.n_undetermined) == (6, 2, 2)
    assert d2.allele_denominator == 12


def test_blank_and_malformed_exons_tested_are_undetermined_not_tested(pat):
    st = by_id(pat, E.record_state(pat, "R1"))
    assert st["P08"] == E.UNDETERMINED and st["P09"] == E.UNDETERMINED
    assert st["P01"] == E.TESTED


@pytest.mark.parametrize("bad", ["R1 ,R2", "R1;R2", " ", "R1,", ",R1", "R1,R2,"])
def test_malformed_values_never_make_an_individual_eligible(bad):
    p = pd.DataFrame({"patient_id": ["X"], "exons_tested": [bad], "year_requested": ["2003"]})
    d, mask = E.denominator(p, "R1")
    assert d.n_tested == 0 and d.n_undetermined == 1 and not mask.any()


def test_unknown_region_makes_everyone_undetermined(pat):
    d, mask = E.denominator(pat, None)
    assert d.n_tested == 0 and d.n_undetermined == 10 and not mask.any()


def test_panel_eligible_needs_every_panel_region(pat):
    m = E.panel_eligible(pat, ["R1", "R2"])
    assert m.tolist() == [True] * 5 + [False] * 4 + [True]
    assert E.panel_eligible(pat, ["R1"]).sum() == 8


# --- the four statuses ---------------------------------------------------------------------
def test_four_genotype_states(pat, var):
    s = by_id(pat, E.genotype_status(pat, var, "VB", "R2"))
    assert s["P01"] == E.DETECTED            # row and region tested
    assert s["P02"] == E.NOT_DETECTED        # tested, no row
    assert s["P06"] == E.NOT_TESTED          # region not in the record (row ignored as a carrier)
    assert s["P08"] == E.UNDETERMINED
    assert s["P09"] == E.UNDETERMINED
    assert set(s.values()) == set(E.STATUSES)


def test_outside_scope_is_not_tested_but_unusable_stays_undetermined(pat, var):
    scope = E.primary_period_scope(pat, SimpleNamespace(primary=SimpleNamespace(
        contains=lambda y: 2001 <= y <= 2006)))
    s = by_id(pat, E.genotype_status(pat, var, "VA", "R1", scope))
    assert s["P07"] == E.NOT_TESTED          # usable record, outside the period
    assert s["P08"] == E.UNDETERMINED
    assert s["P01"] == E.DETECTED and s["P03"] == E.NOT_DETECTED


# --- frequencies --------------------------------------------------------------------------
def test_frequency_uses_eligible_denominator_by_hand(pat, var):
    r = E.cohort_frequency(pat, var, "VA", "R1")
    # eligible 8; carriers P01 (Het), P02 (Hom): alleles 1 + 2 = 3; P08 is an ineligible carrier
    assert r["n_eligible"] == 8 and r["allele_denominator"] == 16
    assert (r["heterozygous"], r["homozygous"], r["carriers"], r["alleles"]) == (1, 1, 2, 3)
    assert r["allele_freq"] == pytest.approx(3 / 16)
    assert r["carrier_freq"] == pytest.approx(2 / 8)
    assert r["wildtype"] == 6
    assert r["n_detected_ineligible"] == 1
    assert r["eligibility_status"] == E.STATUS_EVALUABLE


def test_a_carrier_in_a_region_not_tested_is_not_counted(pat, var):
    r = E.cohort_frequency(pat, var, "VB", "R2")
    assert r["n_eligible"] == 6 and r["alleles"] == 1            # P01 only
    assert r["allele_freq"] == pytest.approx(1 / 12)
    assert r["n_detected_ineligible"] == 1                       # P06 recorded but R2 not tested


def test_no_eligible_individual_gives_nan_never_zero(pat, var):
    r = E.cohort_frequency(pat, var, "VX", "R9")
    assert r["n_eligible"] == 0 and r["eligibility_status"] == E.STATUS_NO_ELIGIBLE
    assert math.isnan(r["allele_freq"]) and math.isnan(r["carrier_freq"])
    assert math.isnan(r["wildtype"])


def test_a_real_zero_is_zero_not_nan(pat, var):
    r = E.cohort_frequency(pat, var, "VZ", "R1")                 # tested, no carriers
    assert r["allele_freq"] == 0.0 and r["carrier_freq"] == 0.0 and r["wildtype"] == 8


@pytest.mark.parametrize("source,status", [("ambiguous", E.STATUS_REGION_AMBIGUOUS),
                                           ("unknown", E.STATUS_REGION_UNKNOWN)])
def test_unknown_or_ambiguous_region_is_not_evaluable(pat, var, source, status):
    r = E.cohort_frequency(pat, var, "VA", None, region_source=source)
    assert r["eligibility_status"] == status and math.isnan(r["allele_freq"])


def test_period_scope_restricts_the_denominator_and_the_numerator(pat, var):
    design = SimpleNamespace(primary=SimpleNamespace(contains=lambda y: 2001 <= y <= 2006))
    scope = E.primary_period_scope(pat, design)
    assert scope.tolist() == [True] * 5 + [False] * 5
    r = E.cohort_frequency(pat, var, "VA", "R1", scope=scope, scope_label="primary")
    assert (r["n_in_scope"], r["n_tested"], r["n_undetermined"]) == (5, 5, 0)
    assert r["allele_denominator"] == 10 and r["alleles"] == 3
    assert r["allele_freq"] == pytest.approx(3 / 10)
    assert r["n_detected_ineligible"] == 1 and r["scope"] == "primary"


def test_more_carriers_than_eligible_individuals_is_an_integrity_error(pat):
    bad = rows(*[(f"P{i:02d}", "VQ", "R1", "Het") for i in range(1, 8)],
               *[(f"P{i:02d}", "VQ", "R1", "Het") for i in range(1, 8)])
    from mefv.validation import DataIntegrityError
    with pytest.raises(DataIntegrityError):
        E.cohort_frequency(pat.iloc[:5], bad, "VQ", "R1")


# --- region resolution --------------------------------------------------------------------
def test_region_comes_from_observed_rows_then_a_known_dictionary_region(var):
    obs = E.variant_regions(var)
    assert E.resolve_region("VA", obs, {"VA": "R9"}, {"R1", "R2"}) == ("R1", "cohort_observed")
    assert E.resolve_region("VN", obs, {"VN": "R2"}, {"R1", "R2"}) == ("R2", "dictionary")
    assert E.resolve_region("VN", obs, {"VN": "R9"}, {"R1", "R2"}) == (None, "unknown")
    assert E.resolve_region("VN", obs, {}, {"R1"}) == (None, "unknown")


def test_conflicting_regions_are_ambiguous_and_not_guessed():
    v = rows(("A", "VM", "R1", "Het"), ("B", "VM", "R2", "Het"))
    obs = E.variant_regions(v)
    assert E.resolve_region("VM", obs, {"VM": "R1"}, {"R1", "R2"}) == (None, "ambiguous")


def test_known_regions_come_only_from_usable_records(pat):
    assert E.known_regions(pat) == {"R1", "R2"}
