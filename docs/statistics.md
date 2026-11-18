# Statistical methods and their limits

Implemented in `src/mefv/statistics.py` and `src/mefv/comparison.py`. These are
descriptive and exploratory methods for a referred clinical series. None of
them is clinically validated.

## Frequencies and intervals

* **Allele and carrier frequencies** use the eligible denominator
  ([`eligibility.md`](eligibility.md)).
* **Wilson score interval** for a binomial proportion. With zero observations the
  interval is NaN, not [0, 0].
* **Hardy–Weinberg**: exact test (Wigginton et al. 2005). Deviation is expected
  for variants that caused referral; do not pool ascertained and non-ascertained
  variants.

## Linkage disequilibrium

* **Haplotype frequencies** by EM from unphased genotypes. Genotypes must be 0,
  1 or 2 copies; anything else raises an error. The routine reports whether it
  converged (`em_converged`, reported beside every D′/r² in Table 4, the pairwise table and the block table; `not_run` when there were no individuals); an unconverged estimate is flagged, not hidden, and the stage 04 log states how many did not converge.
  EM assumes Hardy–Weinberg at the haplotype level.
* **D′ and r²** are undefined (NaN) when either locus is monomorphic
  (frequency within 1e-9 of 0 or 1).
* **Bootstrap interval for D′** (percentile). In a resample where either locus
  is monomorphic D′ is undefined, so that replicate cannot contribute. Such
  replicates are **dropped and counted**. The tables report requested, valid and
  dropped replicates (`boot_replicates_*`). For rare variants the interval is
  therefore *conditional on the replicates in which both loci segregate*. It is
  not an interval for the full sample. No minimum number of valid replicates is
  imposed (any threshold would be arbitrary); read the counts beside the
  interval. With no valid replicate the limits are NaN.
* **Benjamini–Hochberg q-values** ignore missing p-values and leave them missing.

## Comparison with a population reference

Four situations that were previously collapsed are now kept apart:

| Situation | Example | Reported as |
| --- | --- | --- |
| zero | reference frequency is exactly 0 | `undefined_reference_zero` for ratios (not "infinite") |
| missing | reference has no value | `not_evaluable_reference_missing` (never read as 0) |
| undefined | ratio with a zero denominator | NaN with a status |
| not evaluable | cohort frequency undefined (no eligible individuals) | `not_evaluable_cohort_frequency_undefined` |

* **Highest group.** The group with the highest reference frequency is chosen
  among groups that have a value. If any group's value is missing,
  `highest_group_complete` is `False`: the maximum over the remaining groups can
  understate the true maximum, so an enrichment against it can be overstated.
* **Missing `grpmax`.** A variant with no reference `grpmax` filtering allele
  frequency has BA1 and BS1 `not_evaluable`, never "no".
* **Population-affinity (Table 7).** Only variants with a defined cohort
  frequency and a reference value in every group enter, and the number used is
  reported (`n_variants_compared`, `n_variants_requested`). Fewer than two
  variants gives NaN. The 1e-6 offset added inside the logarithm (to keep log10
  defined at a frequency of exactly 0) is an inherited choice **that has not
  been validated**; correlation and mean log-ratio values depend on it for
  rare variants.

## Panel performance (Table 9)

The panel is ranked and evaluated on the same individuals (in-sample). The
output is labelled `in-sample, exploratory` and must not be read as an
estimate of how a panel would perform on new patients. The percentage is of
**eligible** individuals (`pct_of_eligible_individuals`), not the whole cohort.

## Functional concordance (Table 8)

Variants excluded from the analysis are tallied by reason in
`table8b_functional_exclusions.csv`. Comparisons with an empty group return NaN
rather than a spurious p-value. With a handful of variants the result is
exploratory.

## What is not claimed

No output states certainty that the data cannot support. Intervals describe
sampling variability in a referred series, not population frequencies.
Statistical significance is not clinical significance.

## Empty groups and missing values

A group with no individuals gives a blank or NaN value with an explicit status, never
zero and never a p value of 1: frequencies and percentages are blank (`n/a` in Table 1),
the chi-square tests of Table 1 are not defined (NaN) when a table has fewer than two rows
or columns or an empty margin, linkage statistics are NaN with `em_converged = not_run`,
the Mann-Whitney test is NaN when a group is empty, and a figure with nothing to draw
says so instead of failing. `None`, NaN and `pd.NA` are all read as missing in the
reference comparison and in `exons_tested`; a missing `exons_tested` is unusable
(undetermined), never the region called "nan". These behaviours are exercised on
degenerate synthetic inputs (no variant rows, nobody eligible, one variant row) in
`tests/test_pipeline_integration.py`.
