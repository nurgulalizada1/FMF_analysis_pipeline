# ACMG/AMP evidence: configuration and provenance

Stage 09 assembles an evidence package and blank assessor forms. **The pipeline
assigns no ACMG code and no pathogenicity classification.** Classification is a
human, multi-assessor decision.

## Provisional rules carried over from the original analysis

| Item | Value | Status |
| --- | --- | --- |
| BA1 threshold | reference `grpmax` filtering AF strictly above 0.05 | provisional |
| BS1 threshold | reference `grpmax` filtering AF strictly above 0.01 | provisional |
| Codes pre-filled `n/a` on the forms | PS2 PM3 PM6 PP1 PP4 PP5 BS4 BP2 BP5 BP6 | provisional |

These values were present in the original analysis. This repository does not
document who set them, from which guideline or specification, or on what
evidence. BS1 thresholds are normally disease-specific. They are kept so that
results do not change silently, written beside every evidence package in
`classification/acmg_rule_provenance.tsv`, and listed as owner decisions
([`owner_decisions.md`](owner_decisions.md)). They are defined in
`src/mefv/comparison.py` and `src/mefv/acmg.py`.

## Conservative rules added

* A variant without a reference `grpmax` is `not_evaluable` for BA1 and BS1, never
  "criterion not met".
* The reference source, version and retrieval date are read from the reference
  file's provenance header and written into tables 6 and the evidence package.
  They are no longer constants in the code.

## What the evidence package is

A per-variant summary of cohort frequency (with its eligibility status),
reference frequencies, ClinVar and in-silico fields as transcribed, and the
flags above. It is an aid to reviewers, not a result.
