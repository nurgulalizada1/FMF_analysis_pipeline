# Input validation (stage 01)

Stage 01 checks the three raw tables **before** it builds anything. The checks
live in `src/mefv/validation.py`; this page lists every rule and the policy
that goes with it.

The validator **diagnoses; it does not repair.** It never normalises a label,
merges or deletes a record, chooses between conflicting records, or decides
which of two numbers is right. Rows that are entirely empty (trailing
separators) are the only thing removed, and they are counted (`P1a`, `V1a`,
`M3`).

## What a run does

| Severity | Meaning | Effect |
|---|---|---|
| `fatal` | structural defect (blank or malformed identifier, duplicate patient rows, invalid zygosity, orphan variant rows, unparsable year, ...) | stage 01 stops with exit code 2; **cannot be waived** |
| `blocking` | ambiguous integrity problem: exact duplicate variant rows (`V5`, `V6a`), conflicting genotypes (`V6b`) and variant labels that differ only by case or whitespace (`V7`) | stage 01 stops with exit code 2 and writes no analysis tables until the owner has recorded an explicit disposition for each finding. Acknowledging a rule by identifier does **not** release the run |
| `warning` | inconsistency that needs a human decision but does not stop the run (`n_variants`, `variant_detected`, `exons_tested`, region and label inconsistencies) | reported |
| `info` | for completeness (empty rows removed, restored individuals, ...) | reported |

Every run writes `validation_report.json` next to the processed tables. It is
part of the run manifest.

### Diagnostic-only mode

```
python scripts/01_build_dataset.py --diagnostic-only
```

Runs the checks, prints the findings, writes **only** the validation report and
no analysis table, and exits 0 (no fatal or blocking finding) or 2. The report
is marked `diagnostic_only` and `outputs_released: false`, and it replaces any
earlier report: tables left from an earlier run are therefore **not** released,
and every later stage refuses them until stage 01 is run again. It is meant for reviewing a new export before any result is produced.

### Blocking findings: acknowledgement is not release

A blocking finding is a scientific conflict the software cannot resolve: it
cannot tell which of two conflicting genotypes is right, whether an exact
duplicate is a double entry or a real second observation, or whether two
spellings of a label are one variant.

| State | How it arises | Status in the report | Final outputs |
| --- | --- | --- | --- |
| unacknowledged | a blocking rule has findings | `blocked` | **withheld**, exit 2 |
| acknowledged | the rule is named with `--accept-reviewed-findings` (or `MEFV_ACCEPT_REVIEWED_FINDINGS`) | `acknowledged_unresolved` | **withheld**, exit 2 |
| owner disposition recorded | a valid disposition file covers every blocking finding | `owner_disposition_recorded` | written, flagged |

Naming a rule only marks it as seen. It changes nothing about whether the run
proceeds.

#### Stale outputs after a stopped run

* When stage 01 stops (fatal, blocked, acknowledged-only, or a bad disposition
  file) it **removes** `patients.csv` and `variants_long.csv` left by an earlier
  run (and, if the stop happened before a report could be written, the earlier
  report), so they cannot be read as current.
* Stages 03, 04, 05, 06, 08 and 09 **fail closed**: they read the processed
  tables only if the validation report beside them exists, is a readable JSON
  object, and states `"outputs_released": true` (the JSON value; a string or a
  number does not count). A missing, empty, malformed or non-authorising report
  is refused, as is `false` (which also covers a diagnostic-only run). Only
  library and test code can opt out of the "no report" case, with
  `allow_unreported=True`; no pipeline stage does.
* When any stage of `run_pipeline.py` fails, the runner writes
  `RUN_FAILED_outputs_are_not_current.txt` into the tables, figures,
  classification and processed directories. Tables and figures from an earlier
  run are not deleted (use `--clean` for that), but the note says they are not
  results of the current inputs. The note is removed at the start of the next
  *full* run only; a partial run (`--only`, `--start`) never removes it. It is never
  collected as an output. The manifest records
  `run_complete`: true only for a full run in which every stage succeeded.

#### Owner disposition

The run is released only by an explicit disposition, recorded in a private JSON
file **kept outside the repository** (a path inside it is refused) and passed
with `--owner-disposition FILE` (or `MEFV_OWNER_DISPOSITION`; `run_pipeline.py`
accepts the same option). The three digests and the counts are copied from the
`validation_report.json` of a `--diagnostic-only` run:

```json
{"input_sha256":   "<from validation_report.json>",
 "config_sha256":  "<from validation_report.json>",
 "ruleset_sha256": "<from validation_report.json>",
 "dispositions": [
  {"rule": "V7", "records": 12, "decision": "use_as_recorded",
   "decided_by": "<name>", "decided_on": "YYYY-MM-DD",
   "basis": "<why the data are to be used as recorded>"}
]}
```

* **What the disposition is bound to.** `input_sha256`: the three raw tables, byte
  for byte. `config_sha256`: the effective study settings (the parsed study
  parameters: identifier scheme, periods, panel regions, assays, nomenclature
  corrections, expected cohort size, restored-individual setting ...); it ignores
  formatting and key order but changes with any value. `ruleset_sha256`: the rule
  catalogue (identifiers, severities, statements) and its version. A change to any
  of the three, even one that leaves the finding counts unchanged, voids the
  disposition: the run stops again and a new review is needed.
* One entry per blocking rule that has findings. Every one must be covered.
* `records` is the count the report shows: an integer of 5 or more that equals the
  current count, or the text `"<5"` for a finding with fewer than 5 records (small
  counts are never shown, see below). The input digest pins the exact data, so
  `"<5"` cannot be used for a different export.
* The file is refused if it repeats a JSON key, has an unknown key or field, holds a
  non-finite number, has a malformed digest, or has a `decided_on` that is not a real
  calendar date (`2026-02-30` is rejected).
* The only supported decision is `use_as_recorded`: every row is counted
  exactly as it is in the export. Nothing is removed, merged, normalised or
  selected, and the conflict stays in the data and in the report. It means the
  owner accepts that outputs that **count rows** are computed from the data as
  they are, not that the conflict is scientifically resolved.
* **Patient-level genotype analyses stay blocked.** Linkage disequilibrium, the
  genotype cross-tabulation and the pairwise tests need exactly one genotype per
  individual and variant. `mefv.genotype_vector` therefore raises if an individual has
  more than one row for a variant (repeated or conflicting), whatever the
  disposition says: which genotype such an individual has is not defined by
  `use_as_recorded`. The linkage stage stops, the run is marked failed, and nothing
  is chosen, merged or dropped. This needs an owner decision (see
  [`owner_decisions.md`](owner_decisions.md), item 3): correct the export, or
  state a derivation rule that the software may then implement.
* The validation report records the rule, the count as shown, the decision and the
  SHA-256 of the file (so a change to the file is visible in the run manifest).
  The names and the free text stay in the private file. Table 1 shows the
  status and one `Owner disposition` row per rule.
* A malformed file is an error that names the entry and the field, never a
  value. `--diagnostic-only` never applies a disposition.

### Privacy of the diagnostics

**Small counts are suppressed.** In the report (text and JSON), the printed
summary, the stage 01 cohort summary and every error message, a count of 1 to 4
is shown as `<5`; zero is shown as zero. The software keeps the exact counts to
make its decisions. This covers diagnostic output (validation, scanner, error
messages, the run log); it does not cover the study's analytical result tables,
whose disclosure control before sharing is an owner decision
([`owner_decisions.md`](owner_decisions.md)).

The report, the printed summary and every error message contain rule
identifiers and counts only: never an identifier, variant label, region name or
other cell value. Invalid categorical values are described by shape (`blank`,
`case_or_whitespace_variant_of_valid`, `unrecognised`). The same holds for the
integer parsing in `mefv.utils.load`.

### Identifier canonicalisation

An identifier is canonical when it is the configured prefix (in any letter case,
at most once, at the start) followed by ASCII digits only, written with the
configured zero-padding. Anything else is not repaired: an embedded or repeated
prefix, a separator, a sign, a blank or a missing value (None, NaN, NA) has no
canonical form (rules `P2`, `V10`, `M1`). Different raw identifiers that would
share one canonical form are fatal (`P3`). Numbers longer than the padding width are
kept whole, never truncated.

## Rules

| Rule | Severity | Table | Statement |
|---|---|---|---|
| `S1` | fatal | any | every required column is present |
| `P1a` | info | patients | rows that are entirely empty are removed and counted |
| `P1b` | fatal | patients | a row with content has a patient_id |
| `P2` | fatal | patients | patient_id is in canonical form (configured prefix + zero-padded number) |
| `P3` | fatal | patients | distinct raw identifiers do not collapse to one canonical identifier |
| `P4a` | fatal | patients | no individual has identical duplicate rows |
| `P4b` | fatal | patients | no individual has more than one distinct row |
| `P5a` | fatal | patients | year_requested is a non-negative integer |
| `P5b` | fatal | patients | year_requested lies in a configured study period |
| `P6` | warning | patients | sex is one of the permitted labels |
| `P7` | warning | patients | variant_detected is one of the permitted labels |
| `P8` | fatal | patients | n_variants is a non-negative integer |
| `E1` | warning | patients | exons_tested is not blank |
| `E2` | warning | patients | exons_tested is comma-separated, with no whitespace or empty token |
| `E3` | warning | patients | each exons_tested token uses only the permitted characters |
| `E4a` | warning | patients | each exons_tested token is a known region |
| `E5` | warning | patients | exons_tested contains every configured panel region |
| `E6` | info | patients | number of distinct exons_tested region sets |
| `V1a` | info | variants | rows that are entirely empty are removed and counted |
| `V1b` | fatal | variants | a row with content has a patient_id |
| `V10` | fatal | variants | patient_id is in canonical form |
| `V2` | fatal | variants | variant label is not blank |
| `V3` | fatal | variants | exon/region is not blank |
| `V4` | fatal | variants | zygosity is exactly 'Het' or 'Hom' |
| `V5` | blocking | variants | no row is an exact duplicate of another |
| `V6a` | blocking | variants | an individual has one row per variant (same genotype repeated) |
| `V6b` | blocking | variants | an individual has one genotype per variant (conflicting genotypes) |
| `V7` | blocking | variants | one spelling per variant label (case/whitespace) |
| `V8` | info | variants | rows touched by documented nomenclature corrections |
| `V9` | warning | variants | a variant is reported in one region only |
| `V11` | fatal | variants | year_requested is a non-negative integer |
| `V12` | fatal | variants | year_requested lies in a configured study period |
| `M1` | fatal | method | Patient_ID can be parsed to the configured identifier scheme |
| `M2` | fatal | method | YEAR is a non-negative integer |
| `M3` | info | method | rows that are entirely empty are removed and counted |
| `X1` | fatal | cross | every variant row belongs to an individual in the patients or method table |
| `X2` | info | cross | variant rows of individuals present only in the method table |
| `X2b` | info | cross | individuals present only in the method table (restored by design) |
| `X3` | fatal | cross | every patient has a method record |
| `X4` | fatal | cross | one method row per individual |
| `X5a` | warning | cross | declared n_variants does not exceed the observed variant rows |
| `X5b` | warning | cross | declared n_variants is not below the observed variant rows |
| `X5c` | info | cross | declared n_variants equals the distinct individual+variant pairs but not the rows |
| `X6` | warning | cross | variant_detected agrees with the declared n_variants |
| `X7` | warning | cross | variant_detected agrees with the observed variant rows |
| `X8` | warning | cross | a variant's region is listed in the individual's exons_tested |
| `X9` | warning | cross | individuals with variant rows have an interpretable exons_tested |
| `X10` | warning | cross | year agrees between the patients and method tables |
| `X11` | warning | cross | sex agrees between the patients and method tables |
| `A1` | fatal | built | every individual has a resolved assay |
| `A2` | fatal | built | every assay label is one of the configured assays |
| `C1` | fatal | built | the cohort has the configured number of individuals |

`records` counts rows; `individuals` counts distinct individuals. For rules that
compare against other tables (`X*`), the counts refer to the individuals in the
patients table.

## Open questions the validator deliberately does not answer

* **Which is authoritative: declared `n_variants` or the observed rows?**
  `X5a`/`X5b` report both numbers and the totals; nothing is decided or
  corrected.
* **Which of two conflicting genotypes is correct?** `V6b` reports them and
  blocks the outputs.
* **Does a missing variant row mean "reference"?** The pipeline's data model
  treats it so (see `src/mefv/utils.py` and `docs/limitations.md`). The
  validator can only check that an individual has a usable `exons_tested`
  record (`E1`-`E5`, `X8`, `X9`); it cannot establish that the individual was
  tested for a given variant.
* **Is a case or whitespace variant of a label the same variant?** `V7`
  blocks the run until the owner reviews it; the pipeline never merges the
  spellings, so a run released by an owner disposition still treats them as
  distinct strings.
* **Eligibility for frequency denominators** is defined from `exons_tested`
  after validation, not by the validator: see [eligibility.md](eligibility.md).
  A blank or malformed `exons_tested` makes an individual *undetermined*.

## Limits

Rules check structure and internal consistency. A clean report means the tables
are well-formed and mutually consistent, **not** that the records are correct or
that the analysis is scientifically valid.
