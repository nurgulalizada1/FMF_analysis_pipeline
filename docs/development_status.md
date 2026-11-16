# Development status

Last updated: 2026-11-16 (second revision)

This file separates what is implemented from what is not, so the repository is
never read as more complete than it is. It quotes no study results or study
settings; those live in private files.

## Legend

| Mark | Meaning |
| --- | --- |
| 🧪 | Implemented and covered by automated tests that run from a clean clone |
| ✅ | Also exercised locally on the private study data (not reproducible from this repository) |
| 🔲 | Planned, not implemented |

## Library — `src/mefv/`

| Module | Status | Notes |
| --- | --- | --- |
| `config.py` | 🧪 | Paths, seeds, pinned build timestamp; validated interface to the private study parameters, with no study defaults |
| `reference.py` | 🧪 | Population reference table: schema, provenance header, value checks |
| `annotation.py` | 🧪 | Annotation schema, merge rules, decision-log validation |
| `statistics.py` | 🧪 | Wilson interval, exact Hardy-Weinberg, EM haplotypes, D′/r², Benjamini-Hochberg |
| `utils.py` | 🧪 | Two-table loading with period assignment, genotype encoding |

## Analysis stages and tools — `scripts/`

| Script | Status | Notes |
| --- | --- | --- |
| `01_build_dataset.py` … `09_classification_package.py` | 🧪 ✅ | Each runs end to end on synthetic data under test, and ran on the study data |
| `run_pipeline.py` | 🧪 ✅ | Runs the stages; writes a manifest with a deterministic and an execution section; `--clean` |
| `check_regression.py` | 🧪 ✅ | Cell- and value-level comparison with an immutable baseline |
| `compare_manifests.py` | 🧪 | Compares the deterministic sections of two run manifests |
| `make_synthetic_data.py` | 🧪 | Invented demonstration data; reads no file |
| `check_clean_clone.sh` | — | Verifies a committed branch from a fresh clone and virtual environment |

🧪 on a stage means it runs end to end on synthetic data in
`tests/test_pipeline_integration.py`: it completes, produces its outputs,
honours the parameters and is deterministic. It does not mean a study result
is scientifically correct.

## Known issues

Recorded rather than silently changed.

1. **Analysis logic is not yet library code.** Each stage mixes computation,
   input/output and formatting. The stages are tested as wholes; individual
   computations inside them are not unit-tested.
2. **Lint rules are relaxed for `scripts/`.** Migrated stages are checked for
   genuine defects but not reformatted, because reformatting migrated analysis
   code risks changing a result. Dismissed findings are documented in
   `pyproject.toml`.
3. **Integration tests take about 30 seconds.** Most of it is two complete
   nine-stage runs, dominated by process start-up and figure rendering. The
   bootstrap is reduced to 200 replicates in these tests (the same code path);
   `make test-fast` skips them.
4. **Python 3.12 is not tested.** CI and the clean-clone check cover 3.11 and
   3.13.

## Changes in the second revision

- **Stage 01 assay provenance.** The `assay_source` label now follows the same
  precedence as the assay value: a record whose assay was set by a confirmed
  override is labelled `laboratory_confirmed` even when the export field was
  not blank (previously `recorded`). Table 1 counts inferred assays explicitly
  by label in both periods. On the study data this changes no table value; a
  synthetic record with a wrong exported assay now tests the case.
- **Stage 07 classification.** The enrichment figure colours variants by the
  dictionary's pathogenic and likely pathogenic tier, the classification every
  other stage uses, instead of a separate list in the study parameters that
  disagreed with the dictionary. On the study data one variant whose
  dictionary tier is VUS is no longer drawn as pathogenic in figure 6a. If its
  classification is revised, the figure follows the dictionary.
- **Release scope.** `release_allowlist.txt` lists every file that may be
  tracked; `tests/test_public_release.py` enforces it and the private-path
  exclusions. No research output is published; see `public_release.md`.

## Changes in the first revision

- Study design - periods, panel regions, assay labels and the assay inference
  rule - moved out of the code into the validated `study_design` section of the
  private parameters. A year outside both periods now stops the pipeline
  instead of being assigned silently.
- The synthetic generator has its own invented design and reads no file.
  Structural tests check its imports, record the files it opens, and check that
  the demonstration reads no private default location.
- The run manifest is written outside the output directories and is never
  collected as an output; it separates deterministic content from execution
  metadata. `--clean` gives a fresh output state; stale files are reported, not
  collected.
- Stage 06 validates the reference table's schema, values and provenance
  header, and can require a specific release.
- The regression check accepts a difference only for an exact cell and exact
  old and new values, requires a decision, a reason and an evidence status for
  each, and fails on unused entries and structural changes.
- Decisions whose evidence has not been checked against the primary source are
  no longer recorded as settled: verification flags that rested on screenshot
  text were restored, and the remaining annotation decision is marked as
  pending second review.
- CI runs on every push and pull request, on Python 3.11 and 3.13, from pinned
  lock files, and runs the demonstration twice with a determinism check.

## Earlier resolved issues

- Stage 02 produced a 12-column dictionary while later stages needed 17, so a
  re-run destroyed reviewed fields; generated and reviewed annotations are now
  separate, with validated merge rules.
- Stage 01 always reported zero nomenclature corrections (counted after
  applying them); the data were always corrected.
- `ld_stats` reported D′ = 0 for a monomorphic locus because an exact
  comparison never fired; a tolerance now detects it.
- Stage ordering, shared environment-variable names for inputs and outputs, and
  non-reproducible PDF timestamps were fixed.

## Planned work

| Item | Depends on |
| --- | --- |
| Extract analysis logic into library modules | — |
| Cohort data dictionary | Confirmation of column definitions |
| Licence and citation metadata | Choice of licence by the research owner |

## Scientific review outstanding

Not software tasks; changing code cannot resolve them.

- Second review of each logged annotation decision and of each verification
  flag, against the original publications. The publisher's supplementary
  tables could not be retrieved from here, so the decisions remain
  provisional.
- Independent variant classification by two assessors, and their concordance.
- The open questions about the provenance of the population reference values
  (the release behind some values, and the source of the ClinVar values),
  recorded in that file's header. They could not be resolved from here: the
  database cannot be queried from this environment, and the researcher's
  screenshots do not settle them.
- Which research outputs, if any, are released, after a disclosure review
  (`public_release.md`).
- How records without a result field, and variants outside the panel, are
  treated.
