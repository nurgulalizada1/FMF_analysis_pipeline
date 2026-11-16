# Privacy policy for this repository

This project analyses human genetic data from clinical testing records. This
document states what may and may not enter the repository, and why.

## Why removing names is not enough

Individual-level genetic data remain identifying after direct identifiers are
removed. Three mechanisms are relevant here:

1. **Rare variants.** A variant observed in a single individual, combined with
   a test year, a sex and an assay type, can narrow a record to one person,
   especially in a small population served by few testing laboratories.
2. **Combinations.** A list of individuals ordered by how many variants they
   carry is quasi-identifying even without identifiers, because the carried
   combination is itself distinctive.
3. **Genotype as a fingerprint.** With enough markers, a genotype is
   identifying in principle. Patient-level genotype tables are therefore never
   published, regardless of how few columns they contain.

## Five categories of material

| Category | Definition | May be public? |
| --- | --- | --- |
| Raw data | Laboratory exports as received, including direct identifiers | **Never** |
| Processed patient-level data | Cleaned tables with one row per individual | **Never** |
| Synthetic data | Artificially generated records matching the file structure only | Yes, when clearly labelled |
| Aggregate results | Group-level counts, frequencies, figures | Only after disclosure review and investigator approval |
| Documentation | Methods, provenance, limitations | Yes |

## Rules

1. `data/` and `results/` are excluded by `.gitignore` in their entirety.
   Anything that becomes public must be un-ignored explicitly.
2. Aggregate outputs are **not** released by default. A file becomes public
   only when it has been reviewed and recorded as approved, with the reviewer
   and date noted.
3. Where an aggregate table is approved for release, cells with fewer than five
   individuals are suppressed or collapsed, and rare variants are not
   cross-tabulated with year, sex or assay.
4. The repository will not disclose more than the corresponding scientific
   publication, and nothing before the research owners have authorised it.
5. Completed variant-classification forms are not published before assessor
   consensus, because they record individual professional judgements.
6. Screenshots of database queries are not published unless they have been
   checked for incidental identifying content.
7. **Code, documentation, tests and commit messages carry no study results.**
   The list of variants observed, which variants were grouped together or
   analysed as a pair, counts, frequencies, p-values, cohort size and
   record-level corrections are all study results. They live in private,
   untracked files (`study_parameters.json`, `data/`), and the tracked code
   reads them from there. Tests and examples use invented labels only.
8. **Disclosure is judged against everything ever committed, not the current
   tree.** Before a branch is published, every file version and every commit
   message on it is scanned for the study's own variant labels, aggregate
   values and identifiers — not only for generic identifier patterns, which
   miss study results that contain no names.

## Limits of `.gitignore`

`.gitignore` prevents untracked files from being added. It does **not**:

- remove a file that is already tracked;
- remove a file from Git history once it has been committed;
- stop `git add -f`.

If sensitive content is ever committed, removing it in a later commit is not
sufficient — the history must be rewritten and any exposed data treated as
disclosed. Checking `git status` before committing is part of the workflow, not
an optional extra.

## Reporting

If you believe this repository contains material that should not be public,
please open an issue without quoting the content, or contact the author
directly.
