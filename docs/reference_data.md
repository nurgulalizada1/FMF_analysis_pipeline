# Reference data

Stage 06 compares cohort allele frequencies with a population reference table
transcribed from an external database. The values in such a table are public;
the list of variants it covers is a study's observed spectrum, so the table
itself is private and never committed.

## File

| | |
| --- | --- |
| Location | `data/external/population_frequencies.csv`, or `MEFV_POPULATION_FREQUENCIES` |
| Template | `data/external/population_frequencies.example.csv` |
| Validated by | `mefv.reference.load_population_frequencies` |

## Provenance header (required)

The file opens with `#` lines. Three must be present as `# key : value`:

| Key | Content |
| --- | --- |
| `source` | the database and, where relevant, the gene and transcript |
| `version` | the database release |
| `retrieved` | when the values were obtained |

Other `# key : value` lines are kept; `verified` and `open` are recommended,
for what has been checked against the source and what has not. A line starting
with `#` and two or more spaces continues the previous value.

The header is printed by stage 06 and copied into the `execution` section of
the run manifest. If the study parameters set
`reference_data.population_frequencies_version`, the header's version must
equal it exactly.

## Columns (required, in this order)

| Column | Meaning |
| --- | --- |
| `variant` | label, matching the analysis tables; unique |
| `mid`, `nfe`, `fin`, `asj`, `sas`, `eas`, `afr`, `amr`, `rmi`, `ami` | allele frequency per genetic ancestry group, from the population table |
| `glob` | global allele frequency |
| `grpmax` | the group-maximum **filtering** allele frequency, the value used for the ACMG BA1 and BS1 criteria |
| `hom` | homozygote count |
| `xx`, `xy` | allele frequency by chromosomal sex |
| `revel`, `cadd` | in-silico scores (may be empty) |
| `cv`, `stars` | ClinVar assertion and review-status stars |

The per-group columns are population-table frequencies, not filtering
frequencies; only `grpmax` is a filtering value. Treating them as the same
would misapply BA1/BS1.

## What stops stage 06

A missing file; a missing or empty `source`, `version` or `retrieved` line; a
version different from the one the study parameters expect; columns that
differ from the list above; an empty or duplicated variant label; a
non-numeric value in a numeric column; a frequency outside [0, 1]. There is no
built-in fallback table.

## What validation does not establish

The checks above concern structure and documentation. They cannot show that
the values were transcribed correctly, or that every value came from the
release the header names. Any such open question belongs in the file's `open`
header line, where it travels with the data; the study's own open questions
are recorded there.
