"""Comparison of cohort allele frequencies with a population reference.

Every function separates four situations that the earlier code collapsed into
one number: a value that is **zero**, one that is **missing** from the
reference, one that is **undefined** (a ratio with a zero denominator), and one
that is **not evaluable** (an input needed to compute it is itself missing). A
missing reference value is never read as zero, an undefined ratio is never
reported as infinite, and a missing ``grpmax`` never means "criterion not met".
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import pandas as pd

#: Population-frequency criteria, carried over unchanged from the original
#: analysis. Their provenance is not documented in this repository; see
#: ``mefv.acmg.RULE_PROVENANCE`` and docs/acmg_provenance.md.
BA1_THRESHOLD = 0.05  # grpmax filtering allele frequency strictly above this: BA1 flagged
BS1_THRESHOLD = 0.01  # ... strictly above this: BS1 flagged

NOT_EVALUABLE = "not_evaluable"


def _missing(x: object) -> bool:
    """None, NaN (any float type) or ``pd.NA``. A missing value is never read as zero."""
    if x is None:
        return True
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):   # a sequence or other non-scalar is not a missing scalar
        return False


def population_flag(grpmax: object, threshold: float) -> str:
    """``yes`` / ``no`` for a frequency criterion, ``not_evaluable`` if ``grpmax`` is missing."""
    if _missing(grpmax):
        return NOT_EVALUABLE
    return "yes" if float(grpmax) > threshold else "no"


def highest_group(
    values: Mapping[str, float], order: Sequence[str]
) -> tuple[str | None, float, bool]:
    """Ancestry group with the highest reference frequency, among those that have a value.

    Returns ``(group, frequency, complete)``. ``complete`` is False if any
    group's value is missing: the maximum over the remaining groups can
    understate the true maximum, so a ratio against it can overstate
    enrichment, and the caller must say so. With no value at all the group is
    None and the frequency NaN.
    """
    present = [k for k in order if not _missing(values[k])]
    if not present:
        return None, math.nan, False
    best = max(present, key=lambda k: values[k])
    return best, float(values[best]), len(present) == len(order)


def enrichment(cohort_af: float, reference_af: object) -> tuple[float, str]:
    """Cohort frequency divided by a reference frequency, with an explicit status.

    ``computed``                      a finite ratio (may be 0 if the cohort frequency is 0);
    ``undefined_reference_zero``      reference frequency is exactly 0 (ratio undefined,
                                      not "infinite");
    ``not_evaluable_reference_missing``  the reference has no value;
    ``not_evaluable_cohort_frequency_undefined``  no eligible individuals in the cohort.
    """
    if _missing(cohort_af):
        return math.nan, "not_evaluable_cohort_frequency_undefined"
    if _missing(reference_af):
        return math.nan, "not_evaluable_reference_missing"
    ref = float(reference_af)
    if ref == 0:
        return math.nan, "undefined_reference_zero"
    return float(cohort_af) / ref, "computed"
