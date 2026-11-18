"""Statistical edge cases and conservative handling of missing reference values."""
from __future__ import annotations

import math

import numpy as np
import pytest

from mefv import acmg
from mefv import comparison as C
from mefv import statistics as S


# --- bootstrap D' ---------------------------------------------------------------------------
def test_bootstrap_with_a_monomorphic_locus_reports_every_replicate_dropped():
    r = S.bootstrap_dprime([0] * 20, [0, 1] * 10, 30, np.random.default_rng(1))
    assert (r["n_requested"], r["n_valid"], r["n_dropped"]) == (30, 0, 30)
    assert math.isnan(r["lo"]) and math.isnan(r["hi"])          # NaN, not 0


def test_bootstrap_accounting_is_consistent_and_conditional_for_rare_variants():
    gA = [1] + [0] * 29
    gB = [1] + [0] * 29
    r = S.bootstrap_dprime(gA, gB, 200, np.random.default_rng(7))
    assert r["n_valid"] + r["n_dropped"] == 200 == r["n_requested"]
    assert r["n_dropped"] > 0 and r["n_valid"] > 0              # the carrier is missed sometimes
    assert 0 <= r["lo"] <= r["hi"] <= 1


def test_bootstrap_is_deterministic_for_a_seed():
    g = ([0, 1, 2, 1, 0, 2, 1, 1, 0, 2], [0, 1, 2, 1, 0, 1, 1, 2, 0, 2])
    a = S.bootstrap_dprime(*g, 100, np.random.default_rng(3))
    b = S.bootstrap_dprime(*g, 100, np.random.default_rng(3))
    assert a == b


def test_bootstrap_with_no_individuals_is_nan():
    r = S.bootstrap_dprime([], [], 10, np.random.default_rng(0))
    assert r["n_valid"] == 0 and math.isnan(r["lo"])


# --- EM / LD --------------------------------------------------------------------------------
def test_perfect_linkage_gives_dprime_and_r2_of_one():
    g = [0, 0, 1, 1, 2, 2, 0, 1, 2, 1]
    r = S.ld_stats(g, g)
    assert r["Dprime"] == pytest.approx(1.0) and r["r2"] == pytest.approx(1.0)
    assert r["converged"] is True


def test_monomorphic_locus_gives_nan_statistics():
    r = S.ld_stats([0] * 10, [0, 1, 2, 1, 0, 2, 1, 1, 0, 2])
    assert math.isnan(r["Dprime"]) and math.isnan(r["r2"]) and math.isnan(r["D"])


def test_empty_input_to_ld_is_nan_not_an_error():
    r = S.ld_stats([], [])
    # EM was not run, so convergence is neither True nor False (it used to say False, which
    # reads as "failed to converge")
    assert math.isnan(r["Dprime"]) and r["converged"] is None


def test_em_reports_non_convergence_when_the_iteration_limit_is_hit():
    gA = [0, 1, 1, 2, 1, 0, 1, 2, 1, 1]
    gB = [0, 1, 1, 2, 0, 1, 1, 1, 1, 1]
    _, ok = S.em_haplotypes_info(gA, gB, iters=1)
    assert ok is False
    _, ok = S.em_haplotypes_info(gA, gB)
    assert ok is True


@pytest.mark.parametrize("bad", [3, -1, 1.5])
def test_invalid_dosage_is_rejected(bad):
    with pytest.raises(ValueError):
        S.em_haplotypes_info([0, bad], [0, 1])


def test_em_rejects_unequal_or_empty_vectors():
    with pytest.raises(ValueError):
        S.em_haplotypes_info([0, 1], [0])
    with pytest.raises(ValueError):
        S.em_haplotypes_info([], [])


def test_haplotype_frequencies_sum_to_one():
    h, _ = S.em_haplotypes_info([0, 1, 2, 1, 1, 0], [1, 1, 0, 2, 1, 0])
    assert sum(h.values()) == pytest.approx(1.0)


# --- multiple testing, intervals ------------------------------------------------------------
def test_bh_fdr_by_hand_with_a_missing_value():
    # finite: 0.01, 0.03, 0.04 (n=3) -> 0.03, 0.045, 0.04 -> running minimum -> 0.03, 0.04, 0.04
    adj = S.bh_fdr([0.01, 0.04, np.nan, 0.03])
    assert adj[0] == pytest.approx(0.03) and adj[1] == pytest.approx(0.04)
    assert adj[3] == pytest.approx(0.04) and math.isnan(adj[2])


def test_bh_fdr_all_missing_and_empty():
    assert np.isnan(S.bh_fdr([np.nan, np.nan])).all()
    assert len(S.bh_fdr([])) == 0


def test_wilson_interval_with_no_observations_is_nan():
    lo, hi = S.wilson_ci(0, 0)
    assert math.isnan(lo) and math.isnan(hi)


def test_wilson_interval_stays_inside_zero_one_at_the_extremes():
    for k, n in ((0, 10), (10, 10)):
        lo, hi = S.wilson_ci(k, n)
        assert 0 <= lo <= hi <= 1


# --- reference comparison -------------------------------------------------------------------
@pytest.mark.parametrize("grpmax,thr,expected", [
    (None, 0.05, "not_evaluable"), (float("nan"), 0.01, "not_evaluable"),
    (0.0, 0.05, "no"),                       # a real zero is a real value
    (0.05, 0.05, "no"), (0.0500001, 0.05, "yes"),   # strictly above the threshold
    (0.011, 0.01, "yes"),
])
def test_population_flag(grpmax, thr, expected):
    assert C.population_flag(grpmax, thr) == expected


def test_missing_grpmax_is_never_reported_as_criterion_not_met():
    assert C.population_flag(None, C.BA1_THRESHOLD) != "no"


@pytest.mark.parametrize("cohort,ref,value,status", [
    (0.02, 0.01, 2.0, "computed"),
    (0.0, 0.01, 0.0, "computed"),
    (0.02, 0.0, None, "undefined_reference_zero"),            # not infinity
    (0.02, None, None, "not_evaluable_reference_missing"),    # not zero
    (0.02, float("nan"), None, "not_evaluable_reference_missing"),
    (float("nan"), 0.01, None, "not_evaluable_cohort_frequency_undefined"),
])
def test_enrichment_statuses(cohort, ref, value, status):
    v, s = C.enrichment(cohort, ref)
    assert s == status
    if value is None:
        assert math.isnan(v)
    else:
        assert v == pytest.approx(value)
    assert not math.isinf(v)


def test_highest_group_reports_incompleteness():
    order = ["a", "b", "c"]
    assert C.highest_group({"a": 0.1, "b": 0.3, "c": 0.2}, order) == ("b", 0.3, True)
    g, f, complete = C.highest_group({"a": 0.1, "b": float("nan"), "c": 0.2}, order)
    assert (g, f, complete) == ("c", 0.2, False)
    g, f, complete = C.highest_group({"a": None, "b": None, "c": None}, order)
    assert g is None and math.isnan(f) and complete is False


# --- ACMG configuration and provenance ------------------------------------------------------
def test_thresholds_are_labelled_provisional_in_the_provenance_table():
    t = acmg.provenance_table().set_index("item")
    assert t.loc["BA1_threshold", "value"] == C.BA1_THRESHOLD == 0.05
    assert t.loc["BS1_threshold", "value"] == C.BS1_THRESHOLD == 0.01
    assert "provisional" in t.loc["BA1_threshold", "evidence_status"]
    assert "provisional" in t.loc["not_applicable_codes", "evidence_status"]
    assert t.loc["classification", "value"] == "none"


def test_not_applicable_codes_are_valid_codes():
    assert acmg.NOT_APPLICABLE <= set(acmg.CODES)
    assert len(set(acmg.CODES)) == len(acmg.CODES) == 28
