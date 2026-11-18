"""Eligibility-aware population model: who may enter which denominator.

A frequency is only meaningful for the individuals who could have shown the
variant. Dividing by the whole cohort assumes everyone was tested for every
region, which is not true when assays differ (a panel assay covers fewer
regions than an extended one). This module is the single place that decides
who counts, so every stage uses the same rule.

Source of truth
    ``exons_tested`` (patients table): the comma-separated regions the
    laboratory recorded as tested for the individual. It is validated by
    :mod:`mefv.validation`; a blank or malformed value is never read as
    "tested".

Four testing states per individual and variant
    ``detected``       a valid variant row exists and the individual's record
                       lists the variant's region;
    ``not_detected``   no variant row, and the record lists the region. This is
                       the only situation in which an absent row is read as
                       reference, and it relies on the documented data-generation
                       rule that individuals with a usable testing record have
                       every observed variant recorded;
    ``not_tested``     a usable record that does not list the region;
    ``undetermined``   the testing record is missing, malformed, or the region
                       of the variant is itself unknown. Nothing is inferred.

Denominators (see docs/eligibility.md)
    n_cohort       every individual in the cohort;
    n_in_scope     the individuals in the analysis scope (the whole cohort, or
                   one study period for period-restricted variants);
    n_tested       in-scope individuals whose record lists the region;
    n_not_tested   in-scope individuals whose usable record omits the region;
    n_undetermined in-scope individuals with an unusable record;
    n_eligible     = n_tested (a variant is eligible when its region was tested
                   in the scope; the period restriction is part of the scope);
    allele denominator = 2 * n_eligible.

A frequency with no eligible individual is NaN with an explicit status. It is
never reported as zero.
"""

from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

from mefv.validation import parse_region_set

DETECTED = "detected"
NOT_DETECTED = "not_detected"
NOT_TESTED = "not_tested"
UNDETERMINED = "undetermined"
STATUSES = (DETECTED, NOT_DETECTED, NOT_TESTED, UNDETERMINED)

TESTED = "tested"  # testing-record state (before looking at the variant rows)

#: ``coverage_source`` written by stage 01 for individuals present only in the
#: method table, whose testing record has not been confirmed by the owner.
RESTORED_UNCONFIRMED = "restored_from_method_coverage_unconfirmed"
#: ... and the value used when the owner has confirmed the panel was tested.
RESTORED_CONFIRMED = "lab_confirmed_baseline_panel"

STATUS_EVALUABLE = "evaluable"
STATUS_NO_ELIGIBLE = "not_evaluable_no_eligible_individuals"
STATUS_REGION_UNKNOWN = "not_evaluable_region_unknown"
STATUS_REGION_AMBIGUOUS = "not_evaluable_region_ambiguous"


# ------------------------------------------------------------- regions ------
def region_sets(pat: pd.DataFrame) -> pd.Series:
    """Per individual: the frozenset of tested regions, or None if unusable."""
    return pat.exons_tested.map(parse_region_set)


def primary_period_scope(pat: pd.DataFrame, design) -> pd.Series:
    """Individuals whose request year lies in the configured primary period.

    Variants configured as period-restricted (``variant_groups.synonymous``:
    recorded completely only in the primary period) are evaluated over this
    scope, so their denominator uses the same time window as their records.
    """
    years = pat.year_requested.astype(str).str.strip().astype(int)
    return years.map(design.primary.contains).astype(bool)


def known_regions(pat: pd.DataFrame) -> set[str]:
    """Every region named in a usable ``exons_tested`` value."""
    out: set[str] = set()
    for s in region_sets(pat):
        if s is not None:
            out |= s
    return out


def record_state(pat: pd.DataFrame, region: str | None) -> pd.Series:
    """``tested`` / ``not_tested`` / ``undetermined`` for one region.

    With ``region=None`` (the variant's region is unknown) every individual is
    ``undetermined``.
    """
    if region is None:
        return pd.Series(UNDETERMINED, index=pat.index)

    def one(s):
        if s is None:
            return UNDETERMINED
        return TESTED if region in s else NOT_TESTED

    return region_sets(pat).map(one)


def panel_eligible(pat: pd.DataFrame, panel_regions: Iterable[str]) -> pd.Series:
    """Individuals whose usable record lists every configured panel region.

    These are the individuals for whom "no variant row" can be read as
    reference across the whole panel, which is what multi-variant analyses
    (linkage, haplotypes, allele burden) require.
    """
    need = frozenset(panel_regions)
    return region_sets(pat).map(lambda s: s is not None and need <= s)


def variant_regions(var: pd.DataFrame) -> pd.DataFrame:
    """Region of each observed variant, from the variant rows.

    ``region_status`` is ``single`` when every row names the same region,
    ``ambiguous`` when rows disagree (validation rule V9).
    """
    g = var.groupby("variant").exon.agg(lambda s: sorted(set(s)))
    return pd.DataFrame(
        {
            "variant": g.index,
            "region": [v[0] if len(v) == 1 else None for v in g],
            "region_status": ["single" if len(v) == 1 else "ambiguous" for v in g],
        }
    ).set_index("variant")


def resolve_region(
    variant: str,
    observed: pd.DataFrame,
    fallback: Mapping[str, str] | None,
    vocabulary: Collection[str],
) -> tuple[str | None, str]:
    """Region of a variant and where it came from.

    Observed rows win. A variant with no rows may use ``fallback`` (the external
    dictionary) only if its region is in ``vocabulary``, the set of regions the
    laboratory records actually use; otherwise the region is unknown.
    """
    if variant in observed.index:
        row = observed.loc[variant]
        if row.region_status == "single":
            return row.region, "cohort_observed"
        return None, "ambiguous"
    r = (fallback or {}).get(variant)
    if r is not None and r in vocabulary:
        return r, "dictionary"
    return None, "unknown"


# ---------------------------------------------------------- denominators ----
@dataclass(frozen=True)
class Denominator:
    n_cohort: int
    n_in_scope: int
    n_tested: int
    n_not_tested: int
    n_undetermined: int

    @property
    def n_eligible(self) -> int:
        return self.n_tested

    @property
    def allele_denominator(self) -> int:
        return 2 * self.n_tested

    @property
    def status(self) -> str:
        return STATUS_EVALUABLE if self.n_tested > 0 else STATUS_NO_ELIGIBLE

    def as_columns(self) -> dict[str, object]:
        return {
            "n_cohort": self.n_cohort,
            "n_in_scope": self.n_in_scope,
            "n_tested": self.n_tested,
            "n_not_tested": self.n_not_tested,
            "n_undetermined": self.n_undetermined,
            "n_eligible": self.n_eligible,
            "allele_denominator": self.allele_denominator,
        }


def denominator(
    pat: pd.DataFrame, region: str | None, scope: pd.Series | None = None
) -> tuple[Denominator, pd.Series]:
    """Denominator for ``region``, and the boolean mask of eligible individuals.

    ``scope`` restricts the population (for example to one study period for a
    period-restricted variant); individuals outside it are not counted in any
    denominator column except ``n_cohort``.
    """
    scope = pd.Series(True, index=pat.index) if scope is None else scope.astype(bool)
    state = record_state(pat, region)
    s = state[scope]
    eligible = (state == TESTED) & scope
    return Denominator(
        n_cohort=len(pat),
        n_in_scope=int(scope.sum()),
        n_tested=int((s == TESTED).sum()),
        n_not_tested=int((s == NOT_TESTED).sum()),
        n_undetermined=int((s == UNDETERMINED).sum()),
    ), eligible


def genotype_status(
    pat: pd.DataFrame,
    var: pd.DataFrame,
    variant: str,
    region: str | None,
    scope: pd.Series | None = None,
) -> pd.Series:
    """The four-state status of every individual for one variant.

    Individuals outside ``scope`` are ``not_tested`` (the variant is not
    recorded for them by design), unless they are unusable, which is
    ``undetermined``.
    """
    scope = pd.Series(True, index=pat.index) if scope is None else scope.astype(bool)
    state = record_state(pat, region)
    has_row = pat.patient_id.isin(set(var.loc[var.variant == variant, "patient_id"]))
    out = pd.Series(UNDETERMINED, index=pat.index, dtype=object)
    usable = (state == TESTED) & scope
    out[usable & has_row] = DETECTED
    out[usable & ~has_row] = NOT_DETECTED
    out[((state == NOT_TESTED) | ((state == TESTED) & ~scope))] = NOT_TESTED
    return out


# ------------------------------------------------------ cohort frequencies --
def cohort_frequency(
    pat: pd.DataFrame,
    var: pd.DataFrame,
    variant: str,
    region: str | None,
    region_source: str = "cohort_observed",
    scope: pd.Series | None = None,
    scope_label: str = "whole cohort",
) -> dict:
    """Counts, frequencies and denominators of one variant in its eligible population.

    Only carriers among eligible individuals enter the numerator. Carriers
    among ineligible individuals (a variant recorded where the record says the
    region was not tested, or where the record is unusable) are not counted
    toward the frequency; they are reported in ``n_detected_ineligible`` so the
    inconsistency stays visible.

    Frequencies are NaN, never zero, when no individual is eligible.
    """
    den, eligible = denominator(pat, region, scope)
    rows = var[var.variant == variant]
    elig_ids = set(pat.loc[eligible, "patient_id"])
    inside = rows[rows.patient_id.isin(elig_ids)]
    het = int((inside.zygosity == "Het").sum())
    hom = int((inside.zygosity == "Hom").sum())
    ineligible_carriers = int(rows[~rows.patient_id.isin(elig_ids)].patient_id.nunique())
    alleles, carriers = 2 * hom + het, het + hom

    if region is None:
        status = STATUS_REGION_AMBIGUOUS if region_source == "ambiguous" else STATUS_REGION_UNKNOWN
    else:
        status = den.status
    evaluable = status == STATUS_EVALUABLE
    out = {
        "variant": variant,
        "region": region if region is not None else "",
        "region_source": region_source,
        "scope": scope_label,
        **den.as_columns(),
        "heterozygous": het,
        "homozygous": hom,
        "carriers": carriers,
        "alleles": alleles,
        "wildtype": den.n_eligible - carriers if evaluable else np.nan,
        "n_detected_ineligible": ineligible_carriers,
        "allele_freq": alleles / den.allele_denominator if evaluable else np.nan,
        "carrier_freq": carriers / den.n_eligible if evaluable else np.nan,
        "eligibility_status": status,
    }
    if evaluable and out["wildtype"] < 0:
        from mefv.validation import DataIntegrityError

        raise DataIntegrityError(
            "carrier count exceeds the eligible individuals for one variant; "
            "run the stage 01 validation (rules V5, V6a, V6b)."
        )
    return out
