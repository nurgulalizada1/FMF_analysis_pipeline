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

## Automated scan and the private denylist

`scripts/check_privacy.py` (library: `src/mefv/privacy.py`) is a safety net
behind review, not a substitute for it.

```bash
python scripts/check_privacy.py                 # tracked files and untracked-not-ignored files
python scripts/check_privacy.py --tracked       # tracked files only (what CI runs)
python scripts/check_privacy.py --staged        # what the next commit would add (pre-commit)
python scripts/check_privacy.py --path results/ # generated outputs, notebooks, any folder
```

**Built-in rules** (`mefv.privacy.RULES`): files under `data/`, `results/` or
`logs/` other than documented placeholders; file types that normally carry
laboratory data (spreadsheets, archives, bundles, patches, databases, images);
file names suggesting exports, patient files, credentials or private lists;
e-mail addresses, international telephone numbers, long digit sequences,
date-of-birth fields, patient-identifier tables, repeated Het/Hom genotype
columns and credential patterns, and symbolic links. A "SYNTHETIC DATA" banner inside
a file exempts nothing: real data can carry the same words. Only the documented
placeholder examples (`data/external/*.example.csv`, `data/reviewed/*.example.csv`) are exempt
from the two table-shape rules, by path. Symbolic links are reported (`PATH-SYMLINK`) and
never followed, so a link cannot carry the scan to a file outside the folder it was given.
E-mail addresses are allowed only for the listed domains and their sub-domains. Binary
files and files over 20 MB are counted as not scanned. If Git cannot be run, a Git command
fails, a path does not exist or a file cannot be read, the scan stops with exit status 2:
it never reports a clean result for files it could not see.

**Private denylist.** Only the data owner knows which strings are sensitive (the
real identifier scheme, names, laboratory file names, study variant labels).
Put them in a text file **outside the repository** and point the scanner at it
with `--denylist PATH` or the environment variable `MEFV_PRIVATE_DENYLIST`.
One entry per line; `#` starts a comment:

```
some literal text            # case-insensitive, at least 4 characters
regex:PATTERN                # a Python regular expression
sha256:HEX                   # SHA-256 of a lower-cased token, so the list need not hold plain text
```

The scanner refuses a denylist inside the repository, refuses one that Git
tracks, and rejects a malformed entry by its line number only. Without a
denylist it runs the built-in rules and says so; `--require-denylist` makes the
absence an error. A denylist that was asked for (`--denylist` or the environment
variable) but does not exist is always an error (exit 2), never a quietly weaker scan. CI never has the denylist, so it can only run the built-in
rules: run the scan with the denylist locally before publishing a branch.

**What a report contains.** File path, rule identifier and a count (a content count below
5 is shown as `<5`; path rules flag a file and carry no count). A matched
value, a line or a snippet is never printed, logged or written, and denylist
entries are named only by line number. The tests assert this with sentinel
strings.

**What it cannot do.** It finds what its patterns and the denylist describe. It
cannot recognise arbitrary personal data, it does not read binary files, and it
does not scan Git history (rule 8 above still applies). A clean result is not a
guarantee that nothing sensitive is present.

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
