# Outputs

All output is written under `results/` (or the directories set by
`MEFV_TABLES`, `MEFV_FIGURES` and `MEFV_ACMG`). Nothing there is tracked. Output
from synthetic data is not a study finding.

## Scientific artifacts

### Tables (`tables/`, semicolon-separated, UTF-8 with BOM)

| File | Stage | Contents |
| --- | --- | --- |
| `table1_cohort.csv` | 03 | individuals per period, sex, assay, detection, recording completeness |
| `table2_variant_frequencies.csv` | 03 | per-variant genotype counts, carrier and allele frequencies with Wilson intervals, exact HWE p, verification flag; panel regions only |
| `supp_hwe.csv` | 03 | the Hardy-Weinberg columns of table 2 |
| `table3_allele_burden.csv` | 03 | alleles per individual under two definitions (P/LP only; all non-synonymous except common benign), by period and assay |
| `table4_primary_pair.csv` | 04 | linkage of the configured primary pair, by stratum: haplotype frequencies, D′ with bootstrap interval, r², Fisher p |
| `table4b_crosstab_primary.csv` | 04 | genotype cross-tabulation of the primary pair, primary period |
| `table5_complex_alleles.csv` | 04 | variant pairs carried by the same individuals: candidate complex alleles |
| `supp_pairwise_associations.csv` | 04 | all pairwise associations with Benjamini-Hochberg q-values |
| `supp_haplotype_blocks.csv` | 04 | the configured haplotype blocks |
| `supp_high_allele_burden.csv` | 04 | **patient-level**: one row per individual with a high allele burden, with identifier and genotype. Built from study data it must never be published |
| `table6_gnomad_comparison.csv` | 06 | cohort frequency against each ancestry group of the reference table |
| `table7_population_affinity.csv` | 06 | similarity of the cohort frequency profile to each ancestry group |
| `table8_functional_concordance.csv` | 08 | per variant: reviewed functional class and group, cohort frequency, enrichment against the reference, filtering frequency, ClinVar and in-silico evidence |
| `table9_panel_performance.csv` | 08 | share of positive individuals detected by panels of increasing size |

### Figures (`figures/`, each as PDF and PNG)

| Figures | Stage | Contents |
| --- | --- | --- |
| `figure1_allele_frequencies`, `figure2_primary_pair_matrix`, `figure3_ld_heatmap`, `figure4_protein_lollipop`, `figure5_allele_burden` | 05 | frequencies, the primary pair's genotype matrix, pairwise LD, variant positions, allele burden |
| `figure6a_enrichment`, `figure6b_ancestry_heatmap` | 07 | enrichment against the reference, coloured by the dictionary's P/LP tier; ancestry-group comparison |
| `figure7_functional_concordance`, `figure8_panel_performance` | 08 | functional concordance; panel coverage |

### Classification material (`classification/`)

Stage 09 writes a per-variant evidence package (`acmg_evidence_package.*`) and
blank assessment forms for two independent reviewers
(`acmg_form_assessor_1.tsv`, `acmg_form_assessor_2.tsv`). The forms are blank by design: completed forms record
individual professional judgements and are not published.

## Intermediates

`processed/patients.csv` and `processed/variants_long.csv` (stage 01), and the
generated annotations and merged variant dictionary (stage 02). They are
deterministic and their digests are recorded, but they are inputs to later
stages rather than results. Built from study data, they are patient-level and
must never be published.

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
| `execution` | time of the run, stage durations, absolute paths, package versions, git commit, reference-data provenance, stale files | **no**, by design |

To compare two runs: `python scripts/compare_manifests.py FIRST.json SECOND.json`.

`--clean` empties the output directories before a full run (only directories
the runner created; hidden files such as `.gitkeep` are kept). Without it,
files in the output directories that the current full run did not write are
listed under `execution.stale_files_not_collected` and left out of the
artifacts.
