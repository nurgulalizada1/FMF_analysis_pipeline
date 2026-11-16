# Limitations

These are limitations of the methods the pipeline implements and of the kind
of data it is designed for: variant calls transcribed from clinical testing
records. They bound what any output can be used to claim. Limitations
specific to a particular study belong in that study's own documentation.

## Data of this kind

**Clinical selection.** Individuals referred for testing on clinical
suspicion are not a population sample. Frequencies computed from such a
series describe the referred series and must not be presented as population
allele frequencies.

**No phenotype.** The pipeline analyses genotypes only. It performs no
genotype-phenotype analysis, penetrance estimate or severity comparison.

**Retrospective records.** Transcribed laboratory records vary in completeness
over time. Where recording practice changed, analyses that depend on the
affected variants must be restricted to the period with complete recording.
The pipeline supports this through a configured primary period; choosing it
correctly is a scientific decision, not a software one.

**Assay information may be inferred.** If the assay is missing from older
records, the optional inference rule (`study_design.assay_inference`) fills it
from a laboratory statement. Such values are labelled as inferred in the
`assay_source` column and must be reported as a statement, not as data.

**Panel restriction.** Treating an individual with no recorded variant as
homozygous reference is valid only for regions tested in every individual.
Cohort-wide tables are therefore restricted to the configured panel regions,
and variants outside them are reported separately. The resulting spectrum is a
panel-restricted spectrum.

**Records without a result field** are counted as "no variant detected" by
the build stage. This is an assumption to be stated wherever it affects a
rate.

## Analysis

**Haplotypes are estimated, not phased.** D′, r² and haplotype frequencies are
estimated from unphased genotypes by EM. They are population-level estimates,
compatible with a *cis* configuration but not demonstrating it in any
individual. Variants carried by the same individuals are reported as
**candidate** complex alleles only; confirming one requires family segregation
or long-read sequencing.

**EM assumes Hardy-Weinberg equilibrium** at the haplotype level. In an
ascertained series this is violated for the variants on which selection acted.
The estimates remain usable, but the assumption is real.

**Hardy-Weinberg deviations are expected for ascertained variants.** For
variants that caused referral, departure from equilibrium follows from
ascertainment, not from genotyping error. Non-ascertained variants are
analysed separately as an internal check; the two groups must not be pooled.

**Classification depends on source and date.** Tiers come from an external
database as transcribed at a particular time. A later database update can
change the diagnostic yield. The transcription date must be recorded with the
source file.

**Diagnostic yield has more than one definition.** A rate over pathogenic and
likely pathogenic variants differs from one over all non-synonymous variants.
Both are reported; they answer different questions and must not be quoted
interchangeably.

**Functional evidence can be sparse and partly circular.** Published functional
assays cover only some variants, and some variants may have served as
calibration controls in the assay that classified them. Concordance analyses
using those variants are supportive at best.

**Small groups.** Comparisons involving a handful of variants are sensitive to
the reclassification of one of them and are exploratory.

## Software

Passing tests show that the code computes what it is meant to compute on the
tested inputs. They do not validate any scientific interpretation, variant
classification or clinical conclusion. The synthetic demonstration shows that
the software runs; its outputs mean nothing biologically. Open technical
issues are listed in `development_status.md`.
