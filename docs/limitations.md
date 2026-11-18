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

**Eligibility depends on the laboratory record.** Treating an individual with
no recorded variant as reference is valid only for regions the individual was
tested for. The pipeline uses `exons_tested` for this and counts only tested
individuals in each denominator ([`eligibility.md`](eligibility.md)). Individuals
with a missing or malformed record are *undetermined* and are excluded from
denominators; they are neither assumed tested nor assumed negative. Panel-wide
analyses use only individuals tested for every panel region, so they describe a
panel-restricted spectrum. The pipeline cannot verify that a recorded region
was sequenced well enough to exclude a variant.

**Individuals restored from the method table only** are undetermined until the
owner confirms their testing record ([`owner_decisions.md`](owner_decisions.md)).

**Period-restricted variants.** Variants recorded completely only in the primary
period are evaluated within it, using existing metadata only. The pipeline
does not know the testing history of individual years.

**Records without a result field** are `undetermined`, not "no variant
detected".

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

**Bootstrap intervals for D′ are conditional.** Replicates in which either locus
is monomorphic have no D′ and are dropped (and counted). For rare variants the
interval is conditional on the replicates in which both loci segregate and can
be misleadingly narrow ([`statistics.md`](statistics.md)).

**ACMG thresholds are provisional.** The BA1/BS1 frequency thresholds and the
`n/a` code list were inherited without documented provenance
([`acmg_provenance.md`](acmg_provenance.md)). The pipeline assigns no
classification.

**Panel performance is in-sample.** Panels are ranked and evaluated on the same
individuals; the result is exploratory and optimistic.

**Population-affinity statistics use an unvalidated offset** (1e-6 added before
taking logarithms) and depend on how many variants have a complete reference.

**Missing reference values are not zero.** A variant without a reference value
has no enrichment ratio; one with a reference of exactly zero has an undefined
ratio. Neither is reported as a number.

**Small groups.** Comparisons involving a handful of variants are sensitive to
the reclassification of one of them and are exploratory.

## Software

Passing tests show that the code computes what it is meant to compute on the
tested inputs. They do not validate any scientific interpretation, variant
classification or clinical conclusion. The synthetic demonstration shows that
the software runs; its outputs mean nothing biologically. Open technical
issues are listed in `development_status.md`.
