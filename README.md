# FMF_analysis_pipeline

Reproducible analysis software for *MEFV* variant data from clinical testing
records in familial Mediterranean fever (FMF).

> **Status: research software under active development.**
> The pipeline runs end to end, is covered by automated tests and ships with an
> invented demonstration dataset. It was written for a study whose data are
> private and are not part of this repository. Nothing here is a scientific
> result; see [what has and has not been checked](#what-has-and-has-not-been-checked).

---

## What it does

The pipeline takes already-reported variant calls from laboratory records - it
does not call variants from sequence data - and, in nine ordered stages:

| Stage | Script | Purpose |
| --- | --- | --- |
| 01 | `01_build_dataset.py` | assemble a two-table analysis dataset from a laboratory export: harmonise identifiers and variant labels, resolve the assay field, restore missing records |
| 02 | `02_variant_annotation.py` | merge transcribed database annotations with a reviewed overlay into a validated 17-column variant dictionary |
| 03 | `03_descriptive.py` | cohort description; variant, allele and carrier frequencies with Wilson intervals; exact Hardy-Weinberg tests; allele burden |
| 04 | `04_linkage.py` | pairwise association with Benjamini-Hochberg control; EM haplotype estimation, D′ and r² with bootstrap intervals; haplotype blocks; candidate complex alleles |
| 05 | `05_figures.py` | descriptive and linkage figures |
| 06 | `06_external_frequencies.py` | comparison with a population reference table, across ancestry groups |
| 07 | `07_external_figures.py` | figures for stage 06 |
| 08 | `08_functional_and_panel.py` | concordance with published functional evidence; diagnostic panel coverage |
| 09 | `09_classification_package.py` | evidence package and blank assessment forms for independent variant classification |

The numeric order is the execution order. Stage 02 builds the dictionary that
stages 03, 08 and 09 read.

Stage 01 first validates the raw tables against a catalogue of rules
([`docs/validation.md`](docs/validation.md)). Every later frequency uses an
eligibility-aware denominator built from the laboratory's `exons_tested`
record ([`docs/eligibility.md`](docs/eligibility.md)): individuals are
`detected`, `not_detected`, `not_tested` or `undetermined`, and a frequency
with nobody eligible is NaN, never zero.

**This is a reproducible computational research pipeline. It is not clinically
validated software and assigns no pathogenicity classification.**

## Public software, private study

The repository holds generic software. Everything that describes a particular
study is supplied locally and never committed:

| Public (in this repository) | Private (supplied locally, never committed) |
| --- | --- |
| pipeline code, schemas and validation | the laboratory export and the processed tables |
| tests and an invented demonstration dataset | `study_parameters.json`: study periods, panel regions, assay labels, variant groups, record-level corrections |
| documentation and templates (`*.example.*`) | transcribed database annotations, reviewed annotations, the decision log |
| | population reference values for the observed variants |
| | regression baselines and the record of accepted differences |

The code holds no study-specific value and no fallback for one. Every private
setting is validated when it is read, and a missing or malformed setting stops
the pipeline with an error naming it. See
[`docs/configuration.md`](docs/configuration.md) for every setting and
[`docs/privacy_policy.md`](docs/privacy_policy.md) for what may be published.
Stage 01 validates the raw tables before building anything; the rules, the
blocking policy and the diagnostic-only mode are in
[`docs/validation.md`](docs/validation.md).

**This repository contains no patient data and no study results.** Every
tracked file must be on [`release_allowlist.txt`](release_allowlist.txt), which
the tests enforce; the list of approved research outputs is currently empty.
Adding one requires the disclosure review described in
[`docs/public_release.md`](docs/public_release.md).

---

## Installation

Python 3.11 or 3.13 (both are tested; 3.12 is expected to work but is not
tested).

```bash
git clone https://github.com/nurgulalizada1/FMF_analysis_pipeline.git
cd FMF_analysis_pipeline
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements-lock-py3.13.txt     # or requirements-lock-py3.11.txt
pip install --no-deps -e .
```

The lock files pin every package, test tools included, to the versions the
test suite was run with. `pip install -e ".[dev]"` also works, with whatever
versions satisfy the ranges in `pyproject.toml`.

---

## Synthetic demonstration

```bash
make demo          # generate the synthetic data and run all nine stages
make demo-check    # run the demo twice and check the outputs are byte-identical
```

Output goes to `results/demo/`. **Results produced from synthetic data are not
study findings.** They show that the software runs.

The demonstration data are invented, not derived:

- the generator reads no file and imports no study configuration - only the
  annotation schema's column names - which `tests/test_synthetic_independence.py`
  checks by inspecting its imports and recording every file it opens;
- the study design is invented: two periods (`SYN-PERIOD-A`, `SYN-PERIOD-B`),
  two assay labels (`ASSAY_T`, `ASSAY_W`) whose use does not depend on the
  period, and region labels `R1` to `R6` instead of exon numbers;
- variant labels (`SYN_P01`, …), coding positions, protein descriptions and rsIDs
  are placeholders; allele frequencies are round numbers chosen by hand;
- identifiers have the form `SYNTH_0001`.

The synthetic laboratory export contains one instance of each problem stage 01
corrects - a record missing from the patient table, a blank assay field resolved
by a confirmed override, a mistyped variant label, and early blank assay fields
resolved by a documented inference rule - plus a linked variant pair, a pair
carried by exactly the same individuals, and variants outside the panel that
only the extended assay observes. Every correction and analysis path therefore
runs. `PROVENANCE.json` records the seed, the invented design and frequencies,
and that no input was read.

---

## Running on study data

1. Copy `study_parameters.example.json` to `study_parameters.json` and fill it
   in ([`docs/configuration.md`](docs/configuration.md)).
2. Put the private inputs where the configuration expects them, or point the
   environment variables at them (table below).
3. Run:

```bash
python scripts/run_pipeline.py --clean     # all stages, fresh output directories
python scripts/run_pipeline.py --from 03   # from stage 03 onwards
```

Stage 01 stops with exit code 2 on structural defects, and on exact duplicate
rows, conflicting genotypes or variant labels that differ by case or whitespace
until the owner has recorded an explicit disposition for each in a private file
outside the repository (`--owner-disposition`). Naming a rule with
`--accept-reviewed-findings` only marks it as seen and releases nothing.
Nothing is ever merged, normalised or deleted automatically
([`docs/validation.md`](docs/validation.md)). Before
publishing anything, run `python scripts/check_privacy.py` locally with a
private denylist kept outside the repository
([`docs/privacy_policy.md`](docs/privacy_policy.md)).

| Variable | Default | Contents |
| --- | --- | --- |
| `MEFV_RAW_AUG`, `MEFV_RAW_VARIANTS`, `MEFV_RAW_PATIENTS` | `data/raw/method.csv`, `variants.csv`, `patients.csv` | laboratory export (stage 01 input) |
| `MEFV_PROCESSED` | `data/processed` | stage 01 output: the two analysis tables |
| `MEFV_SOURCE_ANNOTATIONS` | `data/external/source_annotations.csv` | transcribed database records |
| `MEFV_ANNOTATION_OVERLAY` | `data/reviewed/variant_annotations_reviewed.csv` | reviewed annotations |
| `MEFV_ANNOTATION_DECISION_LOG` | `data/reviewed/annotation_decision_log.csv` | one row per manual correction |
| `MEFV_POPULATION_FREQUENCIES` | `data/external/population_frequencies.csv` | population reference table (stage 06) |
| `MEFV_ANNOTATIONS_GENERATED`, `MEFV_DICT` | `data/external/…` | stage 02 outputs |
| `MEFV_STUDY_PARAMETERS` | `study_parameters.json` | private study parameters |
| `MEFV_TABLES`, `MEFV_FIGURES`, `MEFV_ACMG` | `results/tables`, `results/figures`, `results/classification` | outputs |
| `MEFV_RUN_MANIFEST` | next to the tables directory | run manifest |
| `MEFV_BOOTSTRAP_REPLICATES` | `2000` | bootstrap replicates for D′ intervals |
| `MEFV_ACCEPT_REVIEWED_FINDINGS` | none | comma-separated blocking rule IDs marked as seen (does not release the run) |
| `MEFV_OWNER_DISPOSITION` | none | private JSON file recording the owner's disposition of each blocking finding; must be outside the repository |
| `MEFV_PRIVATE_DENYLIST` | none | private denylist for `check_privacy.py`; must be outside the repository |

---

## Outputs and how to read them

Tables, figures and classification material are the **scientific artifacts**;
the processed tables and the variant dictionary are **intermediates**. Each run
also writes `run_manifest.json`, which has two sections:

- `deterministic` - the SHA-256 of every artifact and intermediate, a combined
  digest, the analysis settings and each stage's exit code. Two runs on the
  same inputs give an identical section.
- `execution` - the time of the run, durations, paths, package versions, the
  git commit and the provenance of the reference data. This is intentionally
  **not** reproducible and is excluded from every reproducibility claim.

[`docs/outputs.md`](docs/outputs.md) lists each table and figure and what it
contains.

---

## Testing

```bash
make check               # lint and the full test suite (about 40 s; the two
                         # end-to-end runs take most of it)
make test-fast           # everything except the end-to-end integration tests
make privacy-check       # scan for material that should not be published
make demo-check          # the demonstration, twice, with a determinism check
make clean-clone-check   # all of the above from a fresh clone and virtualenv
```

The tests need no study data. They cover the statistical functions against
independent reference values, configuration and reference-data validation,
the annotation merge rules, input validation (one test per rule, plus checks
that no value reaches any output), hand-calculated denominators and the four
eligibility statuses, statistical edge cases (monomorphic bootstrap replicates,
missing and zero reference values, non-convergence), the privacy scanner
(sentinel strings are detected and never reproduced), the regression tool's
precision, the synthetic generator's independence, the public-release allowlist
and private-path exclusions, and all nine stages end to end on synthetic data,
run twice into the same directory. Continuous integration
([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs lint, the tests
the privacy scan (built-in rules) and `make demo-check` on Python 3.11 and 3.13 for every push to any branch and
every pull request.

---

## What has and has not been checked

These are different claims and are reported separately.

- **Software correctness.** The statistical functions are tested against
  independently derived values, and every stage runs end to end on synthetic
  data under test. This shows the code computes what it is meant to compute on
  those inputs.
- **Reproducibility.** Repeated runs on the same inputs produce byte-identical
  artifacts and intermediates, PDFs included (seeds are fixed and
  `SOURCE_DATE_EPOCH` is pinned). This is tested on synthetic data in CI and
  was checked locally on the study data.
- **Agreement with the original analysis.** On the study data, the outputs were
  compared cell by cell with the archived outputs of the original scripts
  using `scripts/check_regression.py`. Each accepted difference must match
  exactly one cell, from an exact old value to an exact new value, and must
  cite a recorded decision and its evidence status. The decisions behind the
  accepted differences are still awaiting second review.
- **Study-data status of the remediation.** The eligibility model, validation
  rules and statistical changes were tested on synthetic data only. They have
  not been run on the study data, so no real-study result has been validated by
  them, and the earlier cell-by-cell regression comparison no longer applies
  to the changed outputs until the owner re-runs and reviews it.
- **Scientific validity - not established.** None of the above shows that a
  classification, an assumption or an interpretation is correct. Variant
  classifications, biological conclusions and clinical interpretation require
  independent expert review, which has not yet taken place.

---

## Statistical methods

| Method | Reference |
| --- | --- |
| Binomial confidence interval (Wilson score) | Wilson EB. *J Am Stat Assoc.* 1927;22:209–212 |
| Exact Hardy-Weinberg test | Wigginton JE, Cutler DJ, Abecasis GR. *Am J Hum Genet.* 2005;76:887–893 |
| EM estimation of haplotype frequencies | Excoffier L, Slatkin M. *Mol Biol Evol.* 1995;12:921–927 |
| Normalised linkage disequilibrium D′ | Lewontin RC. *Genetics.* 1964;49:49–67 |
| False discovery rate control | Benjamini Y, Hochberg Y. *J R Stat Soc B.* 1995;57:289–300 |
| Variant classification framework | Richards S, et al. *Genet Med.* 2015;17:405–424 |

Implementations are in [`src/mefv/statistics.py`](src/mefv/statistics.py).

---

## Repository structure

```
.
├── src/mefv/        library: configuration interface, I/O, statistics,
│                    annotation schema, reference-data validation
├── scripts/         the nine stages, the runner, the regression check,
│                    the synthetic generator and verification helpers
├── tests/           unit, configuration, independence and integration tests
├── docs/            configuration, outputs, annotation schema, reference data,
│                    limitations, privacy policy, release scope, status
├── release_allowlist.txt   every file the public repository may track
├── data/            not tracked, except README and *.example.* templates
└── results/         not tracked; all output lands here
```

Known limitations are in [`docs/limitations.md`](docs/limitations.md); open
work is in [`docs/development_status.md`](docs/development_status.md);
decisions that belong to the data owner are in
[`docs/owner_decisions.md`](docs/owner_decisions.md), and the status of each
analytical element (verified, assumed, awaiting confirmation) in
[`docs/scientific_audit.md`](docs/scientific_audit.md). Methods:
[`docs/statistics.md`](docs/statistics.md),
[`docs/acmg_provenance.md`](docs/acmg_provenance.md). Environment and
determinism: [`docs/reproducibility.md`](docs/reproducibility.md). Privacy
safeguards, including the scanner and the private denylist:
[`docs/privacy_policy.md`](docs/privacy_policy.md).

---

## Licence

No licence has been chosen yet. Until a licence file is added, default
copyright applies: the code may be read, but not reused or redistributed
without permission.

## Author

Nurgul Alizada. The study the pipeline was written for is collaborative work;
co-investigators, and ethics-approval details if applicable, will be credited
here only once they are confirmed in writing.
