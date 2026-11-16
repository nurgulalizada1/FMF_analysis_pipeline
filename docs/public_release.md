# Public release scope

## What this release contains

Software, documentation, tests, templates with placeholder values, and an
invented demonstration dataset that is generated at run time and not tracked.
Every tracked file must match `release_allowlist.txt`;
`tests/test_public_release.py` enforces this, and also checks that no data,
results, archives, images or office files are tracked, that private paths are
ignored, and that templates hold placeholders only.

## Research outputs

**No research output is published in this release.** The allowlist of approved
research outputs is empty.

Publication of study results is permitted only subject to anonymisation,
confidentiality and the applicable ethical requirements. Until a disclosure
review against `privacy_policy.md` has been completed and recorded, and the
ethical basis for releasing aggregate results has been confirmed, no table,
figure or classification file built from study data is added.

Before any output is added, it must pass the review below and be listed here
and in `release_allowlist.txt` by name.

| Output | Level | Main disclosure risk | Eligibility |
| --- | --- | --- | --- |
| processed tables, `supp_high_allele_burden.csv`, completed classification forms | individual | direct re-identification | **never** |
| `table2_variant_frequencies.csv`, `supp_hwe.csv`, `supp_pairwise_associations.csv`, `table5_complex_alleles.csv` | per variant | counts of 1 or 2 for rare variants, which with year, sex or assay can single out a person | only after small-cell suppression or aggregation, and review |
| `table1_cohort.csv`, `table3_allele_burden.csv`, `table9_panel_performance.csv` | cohort strata | small strata (one period, one assay) | after review of every cell below the threshold in `privacy_policy.md` |
| `table4*`, `supp_haplotype_blocks.csv`, `table6`, `table7`, `table8` | per variant pair or per variant | low counts in stratified cross-tabulations; unpublished findings | after review, and not before the corresponding publication |
| figures | derived from the tables above | as for their source tables | as for their source tables |

The variant list itself - which variants were observed - is a study result.
It is published only with the article or after explicit approval.

## Recording an approval

An approval to release a specific output must name the output, the reviewer,
the date and the disclosure checks applied. Record it here and add the file
by name to the allowlist in the same commit. A general permission to publish
"subject to" conditions is not, on its own, an approval of a specific file.
