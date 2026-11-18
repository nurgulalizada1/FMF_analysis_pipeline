# Outputs

All output is written under `results/` (or the directories set by
`MEFV_TABLES`, `MEFV_FIGURES` and `MEFV_ACMG`). Nothing there is tracked. Output
from synthetic data is not a study finding.

## Scientific artifacts

### Tables (`tables/`, semicolon-separated, UTF-8 with BOM)

| File | Stage | Contents |
| --- | --- | --- |
| `table1_cohort.csv` | 03 | individuals per period, sex, assay, detection, recording completeness |
| `table2_variant_frequencies.csv` | 03 | per-variant genotype counts, carrier and allele frequencies with Wilson intervals, exact HWE p, verification flag, and the denominator columns (`n_cohort`, `n_in_scope`, `n_tested`, `n_not_tested`, `n_undetermined`, `n_eligible`, `allele_denominator`, `n_detected_ineligible`, `eligibility_status`, `scope`, `region_source`); frequencies are NaN, never 0, when no individual is eligible. See [`eligibility.md`](eligibility.md) |
| `supp_period_scope.csv` | 03 | for period-restricted variants: the scope, rows in scope and rows outside it |
| `supp_hwe.csv` | 03 | the Hardy-Weinberg columns of table 2 |
| `table3_allele_burden.csv` | 03 | alleles per individual under two definitions (P/LP only; all non-synonymous except common benign), by period and assay; panel-eligible individuals only, with `n_excluded_not_eligible` |
| `table4_primary_pair.csv` | 04 | linkage of the configured primary pair, by stratum: haplotype frequencies, D′ with bootstrap interval and its replicate accounting (`boot_replicates_requested`, `_valid`, `_dropped`), EM convergence (`em_converged`: `True`, `False`, or `not_run` when the stratum has no individuals), r², Fisher p, and `n_excluded_not_eligible` |
| `table4b_crosstab_primary.csv` | 04 | genotype cross-tabulation of the primary pair, primary period |
| `table5_complex_alleles.csv` | 04 | variant pairs carried by the same individuals: candidate complex alleles |
| `supp_pairwise_associations.csv` | 04 | all pairwise associations with Benjamini-Hochberg q-values and the EM convergence flag (`em_converged`) of each D′/r² |
| `supp_haplotype_blocks.csv` | 04 | the configured haplotype blocks, with `em_converged` for each D′/r² |
| `supp_high_allele_burden.csv` | 04 | **patient-level**: one row per individual with a high allele burden, with identifier and genotype. Built from study data it must never be published |
| `table6_gnomad_comparison.csv` | 06 | cohort frequency (with the denominator columns) against each ancestry group of the reference table; `x_vs_*_status` distinguishes computed, undefined (reference zero) and not evaluable (reference or cohort value missing); `highest_group_complete` flags a maximum taken over incomplete groups; reference source, version and retrieval date come from the reference file's provenance header |
| `table7_population_affinity.csv` | 06 | similarity of the cohort frequency profile to each ancestry group, with the number of variants compared. The 1e-6 offset in the logarithms is unvalidated ([`statistics.md`](statistics.md)) |
| `table8_functional_concordance.csv` | 08 | per variant: reviewed functional class and group, cohort frequency, enrichment against the reference, filtering frequency, ClinVar and in-silico evidence |
| `table8b_functional_exclusions.csv` | 08 | variants excluded from table 8, by reason |
| `table9_panel_performance.csv` | 08 | share of eligible individuals detected by panels of increasing size (`pct_of_eligible_individuals`; replaces `pct_of_whole_cohort`); `evaluation` states that it is in-sample and exploratory |

### Figures (`figures/`, each as PDF and PNG)

| Figures | Stage | Contents |
| --- | --- | --- |
| `figure1_allele_frequencies`, `figure2_primary_pair_matrix`, `figure3_ld_heatmap`, `figure4_protein_lollipop`, `figure5_allele_burden` | 05 | frequencies, the primary pair's genotype matrix, pairwise LD, variant positions, allele burden |
| `figure6a_enrichment`, `figure6b_ancestry_heatmap` | 07 | enrichment against the reference, coloured by the dictionary's P/LP tier; ancestry-group comparison |
| `figure7_functional_concordance`, `figure8_panel_performance` | 08 | functional concordance; panel coverage |

### Classification material (`classification/`)

Stage 09 writes a per-variant evidence package (`acmg_evidence_package.*`) and
`acmg_rule_provenance.tsv` (the provisional thresholds and their status, see
[`acmg_provenance.md`](acmg_provenance.md)), and blank assessment forms for two independent reviewers
(`acmg_form_assessor_1.tsv`, `acmg_form_assessor_2.tsv`). The forms are blank by design: completed forms record
individual professional judgements and are not published.

## Intermediates

`processed/patients.csv` and `processed/variants_long.csv` (stage 01), and the
generated annotations and merged variant dictionary (stage 02). They are
deterministic and their digests are recorded, but they are inputs to later
stages rather than results. Built from study data, they are patient-level and
must never be published.

`processed/validation_report.json` (stage 01) records the validation status
(`passed`, `passed_with_findings`, `owner_disposition_recorded`,
`acknowledged_unresolved`, `blocked` or `failed`), whether the outputs were
released, any owner dispositions (rule, count, decision and the digest of the
private file) and the rule-by-rule counts (see [`validation.md`](validation.md)). It
holds counts only, never an identifier or cell value, is deterministic, and its
digest is in the manifest. It is not a result.

## The run manifest

`run_manifest.json` sits next to the tables directory, never inside an output
directory, and is never collected as an output.

| Section | Contents | Reproducible? |
| --- | --- | --- |
| `deterministic.artifacts` | path, size and SHA-256 of every table, figure and classification file | yes |
| `deterministic.artifacts_sha256` | one digest over all artifact names and digests | yes |
| `deterministic.intermediates` | the same for the intermediates | yes |
| `deterministic.settings` | seeds, bootstrap replicates, pinned build timestamp | yes |
| `deterministic.stages` | stages run and their exit codes | yes |
| `deterministic.run_complete` | true only for a full run in which every stage succeeded | yes |
| `execution` | time of the run, stage durations, absolute paths, package versions, git commit, reference-data provenance, stale files | **no**, by design |

To compare two runs: `python scripts/compare_manifests.py FIRST.json SECOND.json`.

When a stage fails, every output directory gets a file
`RUN_FAILED_outputs_are_not_current.txt`: files beside it were not produced by a
complete run of the current inputs. Only the next *full* run removes it; a partial run
(`--only`, `--start`) leaves it in place, because it does not regenerate every directory.

`--clean` empties the output directories before a full run (only directories
the runner created; hidden files such as `.gitkeep` are kept). Without it,
files in the output directories that the current full run did not write are
listed under `execution.stale_files_not_collected` and left out of the
artifacts.
