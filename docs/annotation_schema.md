# Annotation dictionary: schema and provenance

The variant annotation dictionary is the reference table that stage 02 builds
and stages 03, 08 and 09 consume. It has 17 columns and is built from two
sources that are deliberately kept apart.

## Why the dictionary is built from two files

Stage 02 reads a transcription of an external variant database
(`data/external/source_annotations.csv`) and emits twelve columns. Re-running the stage reproduces them exactly.

Five further columns cannot be reproduced from code. Four are transcribed from
a published functional-assay study, and one records a human assessment that
may deviate from the database value. These live in a separate overlay that
stage 02 reads but never writes.

An earlier version of this pipeline kept the reviewed values only in the
generated output file. A routine re-run of the annotation stage then silently discarded
them, and the dictionary lost the columns that stages 08 and 09 require. The
split below exists to make that failure impossible.

## Files

| File | Written by | Tracked in Git | Contents |
| --- | --- | --- | --- |
| `data/external/variant_annotations_generated.csv` | stage 02 | no | 12 generated columns |
| `data/reviewed/variant_annotations_reviewed.csv` | a human | no | 5 overlay columns |
| `data/reviewed/annotation_decision_log.csv` | a human | no | one row per manual correction |
| `data/external/variant_dictionary.csv` | stage 02 | no | the merged 17-column dictionary |
| `data/external/source_annotations.csv` | a human | no | the transcribed database records |
| `data/*/*.example.*` | — | **yes** | structure only, no cohort variants |

The working files are untracked because the set of variants they list is itself
a study result. See `privacy_policy.md`.

## Merge rules

1. **Stage 02 never writes the overlay.** A routine re-run cannot overwrite a
   reviewed value.
2. **The variant sets must match exactly.** Every generated variant needs an
   overlay row — empty where no evidence exists — and the overlay may not add
   variants of its own. A mismatch raises `AnnotationError` rather than
   producing silently empty columns.
3. **A reviewed value may deviate from its source only with a logged
   decision.** `speckseq` is derived from `speckseq_assay_level_infevers`; if
   the two differ and the decision log has no entry for that variant and
   column, the merge fails.
4. **A missing overlay stops the stage.** It does not fall back to twelve
   columns.

## Schema

`origin` records how a value came to exist:

- **generated** — transcribed from the external variant database into the
  source-annotation file; reproducible by re-running stage 02.
- **derived** — computed by stage 02 from another column in the same table.
- **external** — transcribed from a published source into the overlay; not
  reproducible from code.
- **reviewed** — a human assessment that may deviate from its source.

| # | Column | Origin | Type | Definition | Empty means | Used by |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | `variant` | generated | string | Variant label as used in the laboratory records; the join key | not allowed | 03, 08, 09 |
| 2 | `infevers_name` | generated | string | Variant name in the external database | no separate name | — |
| 3 | `hgvs_full` | derived (`hgvs_c`) | string | Transcript-qualified HGVS coding description | not allowed | 09 |
| 4 | `hgvs_c` | generated | string | HGVS coding description on NM_000243.3 | not allowed | 09 |
| 5 | `hgvs_p` | generated | string | HGVS protein description | not applicable | 09 |
| 6 | `exon_infevers` | generated | string | Exon as reported by the database | not recorded | 09 |
| 7 | `rsid` | generated | string | dbSNP identifier | no dbSNP entry found | 09 |
| 8 | `infevers_class` | generated | string | Current clinical classification | not listed | 09 |
| 9 | `previous_class` | generated | string | Classification before the 2025 re-evaluation | unchanged | 09 |
| 10 | `comment` | generated | string | Note transcribed with the record | no note | 09 |
| 11 | `tier` | derived (`infevers_class`) | string | P, LP, VUS, LB, B, unclassified | not allowed | 03, 08 |
| 12 | `reclassified_2025` | derived (`previous_class`) | boolean | Whether the classification changed in 2025 | not allowed | — |
| 13 | `speckseq` | **reviewed** (`speckseq_assay_level_infevers`) | string | Functional class used by the analysis | variant not assayed | 08, 09 |
| 14 | `other_functional` | external | string | Other published functional evidence | none found | 09 |
| 15 | `speckseq_assay_level_infevers` | external | string | Functional label as recorded in the database | not assayed | — |
| 16 | `speckseq_integrated_S7` | external | string | Integrated functional class, supplementary table S7 | absent from that table | — |
| 17 | `bronnec2026_revised_S9` | external | string | Revised classification, supplementary table S9 | absent from that table | — |

### Validation enforced in code

`mefv.annotation.validate` checks: column names and order match the schema;
no duplicate variants; `variant`, `hgvs_c`, `hgvs_full` and `tier` are never
empty; `tier` uses only the six permitted values; and `hgvs_full` equals
`"NM_000243.3:" + hgvs_c` for every row.

## Decision log

One row per manual correction, with these fields:

`variant`, `column`, `previous_value`, `accepted_value`, `evidence`, `reason`,
`review_status`, `review_date`, `affects_downstream_results`.

`review_status` accepts `accepted`, `provisional` or `pending_second_review`.
`accepted` means the maintainer of the log has accepted the decision; it does
**not** mean the decision has been independently reviewed. There is
deliberately no `supervisor_approved` or `independently_reviewed` value: no
decision is recorded as reviewed unless that review exists in writing and the
vocabulary is extended knowingly. Stage 02 prints every decision with its
status and flags those that are not yet `accepted`, because results depending
on them are provisional too.

A screenshot of a source, or text recognised from one, is a lead, not a check
against the source. An entry whose evidence rests on such material stays
`provisional` or `pending_second_review` until it has been checked against the
original publication or database record.

### Recorded decisions

The decision log is a private file, because each entry names a variant observed
in the study and states its effect on study results. Its structure is shown in
`data/reviewed/annotation_decision_log.example.csv`.

Each entry separates three things, and the log is expected to keep them apart:

- **verified evidence** — what the cited source actually states, checked
  against the source itself, with a note on how it was checked;
- **inference** — any explanation of *why* the source says it that has not been
  checked, marked as such;
- **provenance of the decision** — who proposed and who accepted it, without
  claiming review that did not take place.

A decision that changes an output table must also appear in the private
accepted-differences record used by the regression check, so the change is
explained wherever it shows up.

### Functional-evidence vocabulary

The reviewed `speckseq` column uses the integrated functional class defined by
the published assay, not any single assay readout. Quoting that publication's
methods: "A single PS3/indeterminate/BS3 score retaining the highest
classification in the three assays was kept with the following exception:
variants classified as PS3 based solely on pregnanolone SpeckSeq assay were
downgraded to indeterminate due to the uncertain clinical relevance of
pregnanolone." (Bronnec et al., *J Exp Med* 2026;223(2):e20251065.) The
supplementary table that lists integrated classes calls the middle class
"Intermediate"; the article text calls it "indeterminate". They are the same
class.

Stage 08 maps the reviewed label to a functional group by lookup, so a single
reviewed change alters both the label and the derived group in its output.

## What this schema does not establish

The validation above is a check on structure and traceability. It does not
establish that a classification is scientifically correct, that the external
database was transcribed without error, or that the functional assay applies to
these variants. Those require independent review, which is tracked in
`development_status.md`.
