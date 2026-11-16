"""Tests for the statistical functions.

Every expected value is derived independently of the implementation: from the
closed-form definition of the statistic, from a worked example in the cited
paper, or from a configuration whose answer is known by construction.

These are software tests. They establish that the code computes the stated
statistic; they say nothing about whether the scientific interpretation of a
result is correct.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from mefv.statistics import bh_fdr, em_haplotypes, fmt_ci, fmt_p, hwe_exact, ld_stats, wilson_ci


# --------------------------------------------------------------- Wilson -----
def wilson_reference(k: int, n: int, z: float = 1.959963984540054):
    """Closed-form Wilson score interval, written out independently."""
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return centre - half, centre + half


@pytest.mark.parametrize("k,n", [(1, 10), (5, 20), (37, 200), (118, 200), (0, 50), (50, 50)])
def test_wilson_matches_closed_form(k, n):
    got = wilson_ci(k, n)
    want = wilson_reference(k, n)
    assert got[0] == pytest.approx(max(0.0, want[0]), abs=1e-12)
    assert got[1] == pytest.approx(min(1.0, want[1]), abs=1e-12)


def test_wilson_interval_contains_point_estimate():
    lo, hi = wilson_ci(37, 200)
    assert lo < 37 / 200 < hi


def test_wilson_is_bounded_and_handles_empty():
    assert wilson_ci(0, 10)[0] == 0.0
    assert wilson_ci(10, 10)[1] == pytest.approx(1.0, abs=1e-12)
    assert all(np.isnan(v) for v in wilson_ci(0, 0))


# ------------------------------------------------- Hardy-Weinberg (exact) ---
def hwe_reference(n_hom_alt: int, n_het: int, n_hom_ref: int) -> float:
    """Exact Hardy-Weinberg p value by direct enumeration.

    Written independently of the implementation under test: it sums the exact
    conditional probability of every heterozygote count with the same minor
    allele count, using the closed-form expression rather than the recurrence
    used in statistics.hwe_exact.
    """
    n = n_hom_alt + n_het + n_hom_ref
    rare = 2 * min(n_hom_alt, n_hom_ref) + n_het
    if n == 0:
        return float("nan")
    if rare == 0:
        return 1.0

    def log_prob(het: int) -> float:
        hom_rare = (rare - het) // 2
        hom_common = n - het - hom_rare
        if hom_rare < 0 or hom_common < 0:
            return -math.inf
        return (
            math.lgamma(n + 1) - math.lgamma(hom_rare + 1) - math.lgamma(het + 1)
            - math.lgamma(hom_common + 1)
            + het * math.log(2)
            + math.lgamma(rare + 1) + math.lgamma(2 * n - rare + 1)
            - math.lgamma(2 * n + 1)
        )

    hets = [h for h in range(rare % 2, rare + 1, 2)]
    logs = {h: log_prob(h) for h in hets}
    m = max(logs.values())
    probs = {h: math.exp(v - m) for h, v in logs.items()}
    total = sum(probs.values())
    observed = probs[n_het]
    return min(1.0, sum(v for v in probs.values() if v <= observed * (1 + 1e-9)) / total)


@pytest.mark.parametrize(
    "hom_alt,het,hom_ref",
    [(1, 5, 994), (3, 40, 357), (21, 56, 269), (0, 1, 199), (10, 10, 10), (2, 2, 400)],
)
def test_hwe_matches_independent_enumeration(hom_alt, het, hom_ref):
    assert hwe_exact(hom_alt, het, hom_ref) == pytest.approx(
        hwe_reference(hom_alt, het, hom_ref), rel=1e-9
    )


def test_hwe_perfect_equilibrium_is_not_significant():
    # p = q = 0.5 in 1000 individuals: 250 / 500 / 250 is the expected split.
    assert hwe_exact(250, 500, 250) > 0.9


def test_hwe_extreme_heterozygote_excess_is_significant():
    # No heterozygotes at all with both homozygotes common is a strong departure.
    assert hwe_exact(100, 0, 100) < 1e-10


def test_hwe_monomorphic_returns_one():
    assert hwe_exact(0, 0, 500) == 1.0


def test_hwe_empty_sample_is_nan():
    assert np.isnan(hwe_exact(0, 0, 0))


# ------------------------------------------------------ EM haplotypes -------
def test_em_recovers_complete_linkage():
    """Every A allele travels with a B allele: D' must be 1 and r^2 must be 1."""
    gA = [2] * 20 + [1] * 20 + [0] * 60
    gB = [2] * 20 + [1] * 20 + [0] * 60
    out = ld_stats(gA, gB)
    assert out["Dprime"] == pytest.approx(1.0, abs=1e-6)
    assert out["r2"] == pytest.approx(1.0, abs=1e-6)


def test_em_recovers_independence():
    """Two loci assorting independently: D' and r^2 must be close to zero."""
    rng = np.random.default_rng(0)
    n = 4000
    gA = rng.binomial(2, 0.3, n).tolist()
    gB = rng.binomial(2, 0.4, n).tolist()
    out = ld_stats(gA, gB)
    assert abs(out["Dprime"]) < 0.1
    assert out["r2"] < 0.01


def test_em_haplotype_frequencies_sum_to_one():
    rng = np.random.default_rng(1)
    gA = rng.binomial(2, 0.25, 500).tolist()
    gB = rng.binomial(2, 0.25, 500).tolist()
    h = em_haplotypes(gA, gB)
    assert sum(h.values()) == pytest.approx(1.0, abs=1e-9)
    assert all(v >= -1e-12 for v in h.values())


def test_em_marginals_match_observed_allele_frequencies():
    gA = [2, 1, 0, 1, 2, 0, 1, 1, 0, 2]
    gB = [1, 1, 0, 2, 2, 0, 0, 1, 1, 2]
    h = em_haplotypes(gA, gB)
    assert h["11"] + h["10"] == pytest.approx(sum(gA) / (2 * len(gA)), abs=1e-6)
    assert h["11"] + h["01"] == pytest.approx(sum(gB) / (2 * len(gB)), abs=1e-6)


def test_ld_is_symmetric_in_its_arguments():
    rng = np.random.default_rng(2)
    gA = rng.binomial(2, 0.2, 300).tolist()
    gB = rng.binomial(2, 0.5, 300).tolist()
    a, b = ld_stats(gA, gB), ld_stats(gB, gA)
    assert a["r2"] == pytest.approx(b["r2"], abs=1e-9)
    assert a["Dprime"] == pytest.approx(b["Dprime"], abs=1e-9)


def test_ld_monomorphic_locus_gives_nan():
    """Linkage disequilibrium is undefined when a locus carries one allele only.

    Regression test for a defect found during migration: the guard compared the
    estimated allele frequency to 0 and 1 exactly, but the EM estimate converges
    to roughly 1e-14 rather than reaching 0, so the guard never fired and D' and
    r2 were reported as 0.0 - which reads as "measured, no linkage" rather than
    "not measurable". The guard now uses a tolerance.
    """
    out = ld_stats([0, 1, 2] * 20, [0] * 60)
    assert out["pB"] == pytest.approx(0.0, abs=1e-9)
    assert np.isnan(out["Dprime"])
    assert np.isnan(out["r2"])


def test_ld_fixed_locus_gives_nan():
    """The same holds when a locus is fixed for the alternative allele."""
    out = ld_stats([0, 1, 2] * 20, [2] * 60)
    assert out["pB"] == pytest.approx(1.0, abs=1e-9)
    assert np.isnan(out["Dprime"])


def test_ld_tolerance_is_far_below_any_attainable_frequency():
    """A single allele in a cohort of 10000 must still be treated as polymorphic."""
    from mefv.statistics import FIXED_TOL
    assert FIXED_TOL < 1 / (2 * 10_000)
    n = 10_000
    gA = [1] + [0] * (n - 1)
    gB = [1] + [0] * (n - 1)
    assert not np.isnan(ld_stats(gA, gB)["Dprime"])


# ------------------------------------------- Benjamini-Hochberg procedure ---
def test_bh_leaves_a_single_p_value_unchanged():
    assert bh_fdr([0.03])[0] == pytest.approx(0.03)


def test_bh_adjusted_values_are_monotone_and_never_smaller():
    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    q = bh_fdr(p)
    assert all(qi >= pi - 1e-12 for qi, pi in zip(q, p))
    order = np.argsort(p)
    assert np.all(np.diff(np.asarray(q)[order]) >= -1e-12)


def test_bh_matches_hand_calculation():
    """q_i = min over j>=i of (n/j) * p_j, with p sorted ascending."""
    p = [0.01, 0.02, 0.03, 0.04]
    n = len(p)
    want = []
    for i in range(n):
        want.append(min((n / (j + 1)) * p[j] for j in range(i, n)))
    assert list(bh_fdr(p)) == pytest.approx(want, abs=1e-12)


def test_bh_is_capped_at_one():
    assert max(bh_fdr([0.5, 0.6, 0.9, 0.95])) <= 1.0


# ----------------------------------------------------------- formatting -----
def test_formatting_helpers():
    assert fmt_ci(0.1234567, 0.7654321, dp=3) == "0.123–0.765"
    assert fmt_ci(float("nan"), float("nan")) == ""
    assert fmt_p(0.0123) == "0.0123"
    assert fmt_p(1.23e-28) == "1.2e-28"
    assert fmt_p(float("nan")) == ""
