# data/

**This directory is empty in the public repository, by design.**

Every subdirectory is excluded by `.gitignore` except this file. The structure
is kept so that the pipeline's expected inputs are documented and so that a
local working copy has somewhere to put them.

## Layout

| Directory | Contents | Tracked in Git |
| --- | --- | --- |
| `raw/` | Laboratory exports exactly as received, including direct identifiers | No — never |
| `interim/` | Intermediate files produced while cleaning | No |
| `processed/` | `patients.csv` and `variants_long.csv`, the analysis tables | No — never |
| `external/` | Reference data transcribed from public databases | Only `*.example.csv` |
| `reviewed/` | Reviewed annotations, the decision log, the accepted regression differences | Only `*.example.*` |
| `synthetic/` | Artificially generated demonstration data | No — regenerate with `make synthetic` |

## Expected input files

Stage 01 reads a laboratory export (`raw/method.csv`, `raw/variants.csv`,
`raw/patients.csv` by default; see `../docs/configuration.md`) and writes the
two analysis tables below. Annotation columns are defined in
`../docs/annotation_schema.md`, the reference table in
`../docs/reference_data.md`, and the private study parameters in
`../docs/configuration.md`. A data dictionary for the cohort tables is not yet
written.

**`processed/patients.csv`** — one row per tested individual. This table is the
denominator for every rate the pipeline reports, including individuals in whom
no variant was found.

**`processed/variants_long.csv`** — one row per variant observed in one
individual. Individuals with no row here are encoded as homozygous reference,
which is valid only within the configured panel regions, which every
individual was tested for.

**`external/source_annotations.csv`** — one row per distinct variant,
transcribed from the external variant database: HGVS nomenclature, exon, rsID
and classification. Stage 02 reads it. It contains no patient data, but the
*set* of variants it lists is itself a study result, so it is held back until
the investigators approve its release.

**`reviewed/variant_annotations_reviewed.csv`** — the functional-evidence
fields and any manual scientific corrections. Stage 02 reads it and never
writes it.

**`reviewed/annotation_decision_log.csv`** — one row per manual correction,
with the evidence and reason. A reviewed value that deviates from its source
without a matching log entry makes stage 02 fail.

**`external/population_frequencies.csv`** — population reference frequencies
for each observed variant, read by stage 06. The values are public, but the
list of variants is the study's spectrum, so the file is private. Its header
must record the database, release and retrieval date, which stage 06 checks;
see `../docs/reference_data.md`.

**`reviewed/accepted_differences.json`** — the differences between the current
output and the archived baseline that are explained by logged decisions. Read
by `scripts/check_regression.py`. Private, because it quotes baseline values.

**`external/variant_dictionary.csv`** — the merged 17-column dictionary that
stage 02 writes and stages 03, 08 and 09 read. See
`../docs/annotation_schema.md`.

Each input has a tracked `*.example.*` template showing the structure with
placeholder values and no cohort variants.

## Why the real data are not here

See `docs/privacy_policy.md`. In short: individual-level genetic data remain
identifying after names are removed, and no study data may be released without
the research owners' authorisation.

## Obtaining the data

Access to the underlying dataset is decided by the study's research owners.
This repository cannot provide it.

## Synthetic demonstration data

`make demo` generates an artificial cohort into `data/synthetic/` — including a
laboratory-style raw export under `data/synthetic/raw/`, the annotation and
reference files, and its own `study_parameters.json` — and runs all nine stages
on it, writing to `results/demo/`. The generated files are not tracked.

Every generated CSV file begins with a banner marking it artificial;
identifiers use a `SYNTH_` prefix; variant labels (`SYN_…`), periods, assay
labels and region labels are invented. `PROVENANCE.json` records the seed, the
invented design and allele frequencies, that the generator read no input, and a
statement that the data are not sampled from or fitted to any real cohort.
Synthetic and real data are never written to the same directory.

Results produced from synthetic data are not study findings.
