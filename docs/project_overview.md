# Project overview

## Background

Familial Mediterranean fever (FMF) is an autoinflammatory disorder caused by
variants in *MEFV*, which encodes pyrin. It is characterised by recurrent
self-limiting episodes of fever and serositis. The variant spectrum differs
markedly between populations, and the clinical interpretation of several common
*MEFV* variants remains contested: for some frequent variants, pathogenicity is
disputed and reporting practice differs between laboratories.

Two consequences follow for diagnostic practice:

1. The **diagnostic yield** of *MEFV* testing depends as much on which variants
   a laboratory counts as reportable as on which variants are present.
2. Variants that occur together on the same chromosome (complex alleles) are
   frequently reported as independent findings, which can change how a result
   is interpreted.

## What this software does

This repository holds a computational workflow for *MEFV* variant data from
individuals referred for FMF testing. The pipeline takes
already-reported variants from laboratory records and produces:

- a cleaned, two-table analysis dataset with harmonised variant nomenclature;
- descriptive frequency statistics with confidence intervals;
- Hardy-Weinberg tests stratified so that ascertained variants (which are
  expected to deviate) are not pooled with non-ascertained ones (which are not);
- pairwise variant association with multiple-testing control;
- two-locus haplotype estimates, D′ and r², with bootstrap confidence intervals;
- comparison of cohort frequencies against population reference data;
- a variant annotation table and a structured evidence package for independent
  variant classification.

## What this software does not do

- It does not call variants from sequence data. Variant calls are the input.
- It does not infer phase. Haplotypes are statistical estimates from unphased
  genotypes.
- It does not assign clinical significance automatically. The classification
  stage produces an evidence package and blank forms for human assessors; the
  assessment itself is a separate, human task.
- It does not analyse phenotype. It has no phenotype inputs, so it performs no
  genotype-phenotype analysis.

## Study designs the software supports

The pipeline is written for a retrospective series of individuals tested for
*MEFV* variants, with these features, all configured privately per study
(`docs/configuration.md`):

- two analysis periods, a primary one and a sensitivity one, for series in
  which record-keeping changed over time;
- a set of panel regions tested in every individual, which makes a
  cohort-wide denominator valid, and an extended assay that may cover more;
- optionally, a documented rule for filling a missing assay field in older
  records from a laboratory statement.

The public documentation deliberately describes no particular study: its
periods, panel, assays, cohort size and results are private pending disclosure
review (`privacy_policy.md`).

## Relationship to the scientific publication

This repository documents the computational workflow. It is not a manuscript
and is not a substitute for peer review. Scientific claims belong in the
publication, where they can be assessed with the full methods and the complete
limitations.

## Roles and credit

The repository author wrote the pipeline. Contributions to the underlying
study, and how they are credited, will be stated here only once the research
owners have confirmed them.
