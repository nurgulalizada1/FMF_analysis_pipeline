# Scientific logic audit: what is verified, what is assumed, what awaits the owner

Status of each analytical element in the remediated pipeline. "Verified" means
the code was checked, by computation on **synthetic** data (hand-calculated
unit tests and an independent recomputation), to do what it states. That is
evidence about the software only. **Synthetic tests are not evidence that a
scientific method is appropriate for the study data**, and no private data were
used, so nothing here validates a real-study result. No scientific
justification is supplied where none exists in the repository.

Categories:
**V** verified computationally (the code computes what it states);
**M** supported by source metadata;
**A** explicit methodological assumption (a choice, not a fact);
**O** requires owner confirmation ([`owner_decisions.md`](owner_decisions.md)).

| Element | Status | Basis and limits |
| --- | --- | --- |
| Denominator arithmetic (`n_tested`, `n_not_tested`, `n_undetermined`, `n_eligible`, `2 × n_eligible`, frequency, wildtype) | **V** | Hand-calculated tests (`tests/test_eligibility.py`); an independent recomputation from the processed tables (no `mefv` code) agreed with Tables 2 and 6 and the ACMG package in all 139 checks on the synthetic run, including the partition `tested + not tested + undetermined = in scope` |
| `exons_tested` as the eligibility source | **M**, **O** | Owner decision D1 (laboratory record). The pipeline cannot verify that a listed region was sequenced adequately |
| "No variant row among tested individuals means not detected" | **A** | Inherited from the data-generation rule that an individual with a usable record has every observed variant recorded. Stated in `eligibility.md`; not testable from the data |
| Blank or malformed `exons_tested` is undetermined and excluded | **A** (conservative) | Nothing is inferred. Verified by tests that no malformed value ever makes an individual eligible |
| Extended-region (non-panel) variant frequencies | **V** (arithmetic), **A** (region rule) | Denominator is the individuals whose record lists the variant's region. Region comes from the observed rows when they agree (**M**); a dictionary region is used only for a variant with no rows and only if the laboratory records use that region. The vocabulary is derived from the data, not a closed list (**O**) |
| Period-restricted variants (`variant_groups.synonymous`) | **M** (the group), **O** (scope rule) | The group and its recording period come from existing metadata. Evaluating them within the primary period, so that out-of-period records leave the numerator as well as the denominator, follows from that metadata and changes their counts materially. The owner should confirm the rule and whether any other variant is restricted |
| Panel-wide analyses use individuals tested for every panel region | **V** (population, 399 of 400 on the synthetic cohort), **A** | Needed for "no row = reference" to hold across the whole panel. The number excluded is reported |
| Restored method-only individuals | **A** (undetermined by default), **O** | Restoration comes from the existing method table. Whether they were tested is unknown, so they are excluded from eligible denominators until the owner sets `restored_individuals.coverage_confirmed` |
| EM haplotype estimation | **V** | On 40 random synthetic datasets the EM log-likelihood was never below that of a direct multi-start optimiser, and haplotype frequencies agreed to 7e-9. D′ and r² matched their formulas exactly, and a phase-known 2×2 case matched hand calculation. Assumes Hardy–Weinberg at the haplotype level (stated limitation). The tolerance (1e-12), iteration limit (1000) and the monomorphic cut-off (1e-9) are implementation choices, not validated thresholds |
| EM convergence flag | **V** | Reported as `em_converged`; non-convergence is tested (iteration limit 1). The pipeline reports it but does not exclude unconverged estimates |
| D′ bootstrap | **V** (accounting), **A** (method) | Replicate counts always sum to the number requested (tested). Monomorphic replicates are dropped and counted, so the interval is conditional on polymorphic replicates. No minimum valid-replicate threshold is imposed because any would be arbitrary; a reader must look at `boot_replicates_valid` |
| Missing and zero reference values | **V** (logic), **M** (owner decision D6) | Missing is never zero; a zero reference gives an undefined ratio, never infinity; a missing `grpmax` makes BA1/BS1 `not_evaluable`. Verified by tests |
| Highest-group enrichment with incomplete groups | **A** | Maximum over available groups, flagged `highest_group_complete = False`; the ratio can be overstated |
| ACMG BA1 / BS1 thresholds (0.05 / 0.01) and the `n/a` code list | **O** | Inherited without documented provenance; written beside every evidence package. The pipeline assigns no classification |
| Table 7 log offset (1e-6) | **O** | Inherited, unvalidated. Its effect on real data has not been assessed |
| `n_variants` (declared) versus observed variant rows | **O** | Both counts are reported (`X5a`/`X5b`, totals); the software treats neither as authoritative and corrects nothing |
| Exact duplicate rows, conflicting genotypes, labels differing by case or whitespace | **O** | Block the run. Release needs an explicit owner disposition bound to the reviewed input files, and then the row-counting stages count the rows exactly as recorded (a duplicated row counts twice); this is not a resolution. The patient-level genotype stages (linkage, cross-tabulation, pairwise tests) refuse an individual with more than one row for a variant, whatever the disposition |
| Pre-filled `n/a` ACMG codes (PS2 PM3 PM6 PP1 PP4 PP5 BS4 BP2 BP5 BP6) | **O** | Inherited without documented provenance |
| Panel performance (Table 9) | **A** | In-sample and exploratory by construction; labelled as such. Not a performance estimate |

## What was not done

* No private study data were read, run or validated.
* The sensitivity of Table 7 to the offset, and of bootstrap intervals to the
  number of valid replicates, was not assessed on study data.
* Whether the owner's `exons_tested` record matches what each assay actually
  covered cannot be established from the repository.
