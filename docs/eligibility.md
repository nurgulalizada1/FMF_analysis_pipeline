# Eligibility and denominators

A frequency is only meaningful for the individuals who could have shown the
variant. Dividing by the whole cohort assumes everyone was tested for every
region. That is not true when assays differ (a panel assay covers fewer regions
than an extended one) or when a testing record is missing. This page states who
enters which denominator. The rules live in `src/mefv/eligibility.py`, the one
place every stage uses.

## Source of truth

`exons_tested` (patients table) is the laboratory's record of the regions
tested for an individual. The data owner confirmed it as a reliable laboratory
record and the primary source of eligibility. It is validated first
([`validation.md`](validation.md), rules `E1`–`E5`, `X8`, `X9`).

* A value is **usable** if it is a comma-separated list of region tokens with no
  blank token, whitespace or unsupported character.
* A **blank or malformed** value is never read as "tested" and never as "all
  regions tested". The individual is **undetermined**.

## The four statuses

For an individual and a variant:

| Status | Meaning |
| --- | --- |
| `detected` | a valid variant row exists and the individual's record lists the variant's region |
| `not_detected` | no variant row, and the record lists the region. The only case in which an absent row is read as reference. It relies on the documented data-generation rule that an individual with a usable testing record has every observed variant recorded |
| `not_tested` | a usable record that does not list the region (or the variant lies outside the individual's scope) |
| `undetermined` | the testing record is missing or malformed, or the variant's region is itself unknown or ambiguous. Nothing is inferred |

## Denominators

| Column | Definition |
| --- | --- |
| `n_cohort` | every individual in the cohort |
| `n_in_scope` | individuals in the analysis scope: the whole cohort, or the primary period for period-restricted variants |
| `n_tested` | in-scope individuals whose usable record lists the variant's region |
| `n_not_tested` | in-scope individuals whose usable record omits the region |
| `n_undetermined` | in-scope individuals whose record is unusable |
| `n_eligible` | `n_tested` |
| `allele_denominator` | `2 × n_eligible` |
| `n_detected_ineligible` | carriers among individuals who are not eligible (for example a variant row where the record says the region was not tested). They are **not** in the frequency; the count keeps the inconsistency visible |

Carrier frequency is `carriers / n_eligible`. Allele frequency is
`(het + 2 × hom) / allele_denominator`. The confidence interval uses the same
denominator.

## Missing is not zero

If no individual is eligible, or the variant's region cannot be determined, the
frequency is **NaN** and `eligibility_status` says why:

| `eligibility_status` | Meaning |
| --- | --- |
| `evaluable` | at least one eligible individual |
| `not_evaluable_no_eligible_individuals` | nobody was tested for the region in scope |
| `not_evaluable_region_unknown` | no observed rows and no usable dictionary region |
| `not_evaluable_region_ambiguous` | observed rows disagree on the region (rule `V9`) |

A frequency of exactly 0 means eligible individuals were tested and none
carried the variant.

## Where a variant's region comes from

1. The region on the observed variant rows, if they all agree.
2. If the variant has no rows, the external dictionary's region, **only** if
   that region is one the laboratory records actually use (the vocabulary of
   regions in usable `exons_tested` values).
3. Otherwise the region is unknown. The variant is `not_evaluable`, never
   guessed.

## Period-restricted variants

Variants listed in the study parameters as period-restricted
(`variant_groups.synonymous`) are recorded completely only in the primary study
period. They are evaluated within the primary-period scope, so numerator and
denominator cover the same time window. `supp_period_scope.csv` records, per
variant, the scope and how many rows fall outside it. Eligibility history
(who was tested in which year) is not invented: only existing metadata is used.
Whether other variants are period-restricted is an owner decision
([`owner_decisions.md`](owner_decisions.md)).

## Panel-wide analyses

Linkage, haplotypes, allele burden and panel performance need "no variant row"
to mean reference across the whole panel. They therefore use only individuals
whose usable record lists **every** configured panel region. The number
excluded is reported (`n_excluded_not_eligible` in tables 3 and 4). Figures 2–5
use the same population.

## Individuals present only in the method table

Individuals who appear only in the method export are restored so that the
cohort total is complete. Because no laboratory testing record exists for them
they are restored with a blank `exons_tested`, `variant_detected =
undetermined`, a missing `n_variants`, and
`coverage_source = restored_from_method_coverage_unconfirmed`. They stay
undetermined until the owner confirms that they were tested for the panel
regions and that no variant was found, by setting
`restored_individuals.coverage_confirmed` to `true` in the private study
parameters ([`configuration.md`](configuration.md)). The pipeline never makes
that decision.

## Limits

Eligibility is only as good as `exons_tested`. The pipeline cannot verify that
a recorded region was actually sequenced well enough to exclude a variant, and
it does not model coverage or quality. Passing validation does not prove the
laboratory record is right.
