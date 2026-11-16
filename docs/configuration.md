# Configuration

The pipeline separates **paths**, which are set with environment variables
(see the README), from **study parameters**, which describe a particular study
and live in a private JSON file. This page documents the study-parameters file.

## The study-parameters file

| | |
| --- | --- |
| Location | `study_parameters.json` in the repository root, or the path in `MEFV_STUDY_PARAMETERS` |
| Tracked in Git | **never** - excluded by `.gitignore` (`study_parameters*.json`) |
| Template | `study_parameters.example.json` (placeholders only; it validates) |
| Read by | `mefv.config` accessors; stages 01 to 08 |

The parameters are validated when they are read. A missing file, invalid JSON,
a missing key or a malformed value raises an error naming the entry. Public
code has **no default** for any of these values, so a missing setting cannot
silently change a result.

The synthetic generator writes its own, invented parameters file next to the
synthetic data; it never reads this one.

## Keys

### `study_design` (required)

Validated by `mefv.config.parse_study_design`.

| Key | Type | Meaning |
| --- | --- | --- |
| `periods.primary` | `{label, first_year, last_year}` | the period used for analyses that need complete recording (for example linkage with a benign marker) |
| `periods.sensitivity` | `{label, first_year, last_year}` | the comparison period |
| `panel_regions` | list of strings | regions (for example exons) tested in every individual, in display order. Cohort-wide frequencies are restricted to them, because only there is "no recorded variant" a valid reference genotype |
| `assays.panel` | string | label of the assay restricted to the panel regions |
| `assays.extended` | string | label of the assay that also covers further regions |
| `assay_inference` | `null` or `{until_year, assay, source_label}` | if set, a blank assay field in a record no later than `until_year` is filled with `assay`, and the record's `assay_source` is set to `source_label`. Use only to encode a documented statement from the laboratory, and report it as such |

Rules: years are integers; the two periods do not overlap and have different
labels; every analysed record must fall inside one of them, or the stage stops;
the two assay labels differ; `assay_inference.assay` is one of them.

### Record handling (stage 01)

| Key | Meaning |
| --- | --- |
| `expected_cohort_size` | number of individuals stage 01 must produce; a different count stops the stage |
| `id_prefix`, `id_width` | identifier scheme: prefix and zero-padding of record numbers |
| `assay_overrides` | `{record id: assay}` for assays the laboratory confirmed individually |
| `nomenclature_fixes` | `{label as recorded: correct label}` for confirmed transcription errors |
| `nomenclature_fix_notes` | `{label as recorded: justification}` |

### Variant groups and analysis plan

| Key | Meaning |
| --- | --- |
| `variant_groups.synonymous` | synonymous variants; excluded from allele burden |
| `variant_groups.block_synonymous` | synonymous variants analysed together as one block |
| `variant_groups.benign_common` | common benign variants excluded from the broad allele burden |
| `variant_groups.unverified` | variants whose original report has not been verified; flagged in table 2 |
| `linkage_plan.primary_pair` | `[A, B]`, the pair for the main linkage analysis |
| `linkage_plan.haplotype_blocks` | `{name: [variants]}` analysed as blocks |
| `linkage_plan.cross_locus_pairs` | `{anchor: [partners]}` for block-to-locus comparisons |
| `affinity_comparison_variants` | variants used for the population-affinity comparison in stage 06 |

### `reference_data` (optional)

| Key | Meaning |
| --- | --- |
| `population_frequencies_version` | if set, the population reference file's `# version` header line must equal this value exactly, so a different release cannot be substituted unnoticed |

## Other private files

| File | Schema |
| --- | --- |
| `data/external/source_annotations.csv` | [`annotation_schema.md`](annotation_schema.md) |
| `data/reviewed/variant_annotations_reviewed.csv` | [`annotation_schema.md`](annotation_schema.md) |
| `data/reviewed/annotation_decision_log.csv` | [`annotation_schema.md`](annotation_schema.md) |
| `data/external/population_frequencies.csv` | [`reference_data.md`](reference_data.md) |
| `data/reviewed/accepted_differences.json` | `scripts/check_regression.py`; template `data/reviewed/accepted_differences.example.json` |

## Keeping study-specific file patterns private

`.gitignore` contains only generic patterns. Patterns that would reveal a
particular laboratory's file names or identifier scheme belong in the
untracked `.git/info/exclude` of your working copy.
