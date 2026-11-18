# Unresolved decisions for the data owner

For what is verified, assumed or awaiting confirmation in each analytical element, see [`scientific_audit.md`](scientific_audit.md).

The remediation applied conservative fallbacks where a scientific or data-owner
decision was needed. None of the items below has been decided by the software
or by the people who wrote it. Each lists the fallback in force and what the
owner needs to do.

| # | Decision | Fallback in force | Owner action |
| --- | --- | --- | --- |
| 1 | Were individuals present only in the method table tested for the panel regions, with no variant found? | They are `undetermined` and excluded from all eligible denominators (`coverage_source = restored_from_method_coverage_unconfirmed`) | If confirmed by the laboratory, set `restored_individuals.coverage_confirmed` to `true` in the private study parameters |
| 2 | Which is authoritative: declared `n_variants` or the observed variant rows? | Both reported (`X5a`/`X5b`); neither is preferred or corrected | State the rule, or correct the export |
| 3 | Exact duplicate rows and conflicting genotypes | Blocking: no final output. Acknowledging by rule name does not release the run. No row is ever deleted, merged or selected automatically | Correct the export, or record an explicit disposition (`use_as_recorded`) in a private file outside the repository, bound to the reviewed input files, study settings and rule catalogue ([`validation.md`](validation.md)). **Even with a disposition, the linkage stage (patient-level genotypes) refuses an individual with more than one row for a variant**: state how one genotype per individual and variant is to be derived, or correct the export |
| 4 | Variant labels that differ only by case or whitespace (`V7`) | Blocking in the same way; never normalised or merged | Correct the labels at source, or record a disposition |
| 5 | Period restriction. The configured group (`variant_groups.synonymous`) is recorded completely only in the primary period; is any other variant period-restricted? | Only existing metadata is used and no testing history is invented. The configured group is now evaluated within the primary period, so its carrier counts and denominators exclude records from other periods (previously all records counted) | Confirm the scope rule; add metadata if other variants are restricted |
| 6 | Vocabulary of extended-only regions | A dictionary region is used for a variant with no observed rows only if it appears in a usable `exons_tested` value | Provide the closed vocabulary if it matters |
| 7 | BA1 / BS1 thresholds (0.05 / 0.01) and the pre-filled `n/a` codes | Kept as provisional; written beside every evidence package ([`acmg_provenance.md`](acmg_provenance.md)) | Confirm, replace, and cite the guideline or VCEP specification |
| 8 | Log-offset 1e-6 in the population-affinity statistics (Table 7) | Kept; documented as unvalidated | Validate or replace |
| 9 | Panel-performance claims | Labelled in-sample and exploratory | Validate on independent data before any claim |
| 10 | Treating validation rules `V2`, `V3`, `V4` as fatal, and the assumed vocabularies for `P6`/`P7` | As implemented | Confirm |
| 11 | Re-review of the study-data regression baseline | Not run | Any accepted-differences record kept for the earlier outputs predates the eligibility model. Run `scripts/check_regression.py` locally on study data and review every changed cell: the denominators, extended-region frequencies and some table columns now differ by design |
| 12 | GitHub Actions pinned by commit SHA | Not done (SHAs could not be resolved here and were not guessed) | Pin `actions/checkout` and `actions/setup-python` |
| 13 | Licence | None added | Choose one |
| 14 | Private denylist for the privacy scanner | Scanner runs with built-in rules only and says so | Create a denylist **outside** the repository and set `MEFV_PRIVATE_DENYLIST` ([`privacy_policy.md`](privacy_policy.md)) |
| 15 | Real-study validation | **Not performed.** No private data were read or run | Run the full pipeline and the scanner locally on the study data |
| 16 | Disclosure control of analytical result tables | Diagnostic output (validation, scanner, errors, run log) never shows a count of 1-4 (`<5`). The study's result tables and their console printouts report exact counts, including small ones (carriers, cross-tabulations, burden groups) | Decide whether, and with what threshold, result tables need cell suppression before they leave the study team |
