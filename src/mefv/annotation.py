"""Variant annotation schema, provenance and merge rules.

The annotation dictionary has two distinct origins, and keeping them apart is
the point of this module.

**Generated annotations** are produced by stage 04 from a transcribed record of
an external database. Re-running the stage reproduces them exactly.

**Reviewed annotations** are not reproducible from code. They carry external
functional-evidence fields transcribed from the literature, and manual
scientific corrections made by a human assessor. They live in a separate
overlay file that stage 04 reads but never writes.

The two are joined on ``variant`` to form the dictionary that later stages
consume. Three rules make the join safe:

1. Stage 04 writes the generated file and the merged dictionary. It never
   writes the overlay, so a routine re-run cannot overwrite a reviewed value.
2. Every variant in the generated table must appear in the overlay, and the
   overlay may not introduce variants of its own. A mismatch is an error, not
   a silently empty column.
3. A reviewed field may deviate from the external field it was derived from
   only when the deviation is recorded in the annotation decision log. An
   undocumented deviation is an error.

Rule 3 is what makes a manual correction traceable: the value, its source, the
evidence and the reason are all recorded, and a value that quietly drifts away
from its source cannot pass validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mefv.utils import _leading_comment_lines

KEY = "variant"


@dataclass(frozen=True)
class Column:
    """One column of the annotation dictionary."""

    name: str
    definition: str
    origin: str  # "generated" | "external" | "reviewed" | "derived"
    dtype: str
    missing: str  # what an empty value means
    derived_from: str | None = None
    consumers: tuple[str, ...] = ()


#: The complete schema. ``origin`` records provenance:
#:
#: generated  transcribed from an external database into stage 04 and emitted
#:            by it; reproducible by re-running the stage
#: external   transcribed from a published source into the reviewed overlay;
#:            not reproducible from code
#: reviewed   a human assessment that may deviate from its source, in which
#:            case the deviation must appear in the decision log
#: derived    computed by stage 04 from another column in this table
SCHEMA: tuple[Column, ...] = (
    Column("variant", "Variant label as used in the laboratory records; the join key.",
           "generated", "string", "not allowed",
           consumers=("02", "08", "09")),
    Column("infevers_name", "Variant name as given by the external variant database.",
           "generated", "string", "no separate database name"),
    Column("hgvs_full", "Transcript-qualified HGVS coding description.",
           "derived", "string", "not allowed", derived_from="hgvs_c",
           consumers=("09",)),
    Column("hgvs_c", "HGVS coding description on NM_000243.3.",
           "generated", "string", "not allowed", consumers=("09",)),
    Column("hgvs_p", "HGVS protein description.",
           "generated", "string", "not applicable to this variant",
           consumers=("09",)),
    Column("exon_infevers", "Exon as reported by the external database.",
           "generated", "string", "not recorded", consumers=("09",)),
    Column("rsid", "dbSNP identifier.",
           "generated", "string", "no dbSNP entry found", consumers=("09",)),
    Column("infevers_class", "Current clinical classification in the external database.",
           "generated", "string", "not listed in the database", consumers=("09",)),
    Column("previous_class", "Classification before the database's 2025 re-evaluation.",
           "generated", "string", "classification unchanged", consumers=("09",)),
    Column("comment", "Free-text note transcribed with the database record.",
           "generated", "string", "no note", consumers=("09",)),
    Column("tier", "Compact classification tier: P, LP, VUS, LB, B or unclassified.",
           "derived", "string", "not allowed", derived_from="infevers_class",
           consumers=("02", "08")),
    Column("reclassified_2025", "Whether the classification changed at the 2025 re-evaluation.",
           "derived", "boolean", "not allowed", derived_from="previous_class"),
    Column("speckseq", "Functional class used by the analysis.",
           "reviewed", "string", "variant not assayed",
           derived_from="speckseq_assay_level_infevers", consumers=("08", "09")),
    Column("other_functional", "Other published functional evidence, as free text.",
           "external", "string", "no other functional evidence found",
           consumers=("09",)),
    Column("speckseq_assay_level_infevers",
           "Functional label as recorded in the external variant database.",
           "external", "string", "variant not assayed"),
    Column("speckseq_integrated_S7",
           "Integrated functional class from the assay publication's supplementary table S7.",
           "external", "string", "variant absent from that table"),
    Column("bronnec2026_revised_S9",
           "Revised clinical classification from the assay publication's supplementary table S9.",
           "external", "string", "variant absent from that table"),
)

COLUMNS: tuple[str, ...] = tuple(c.name for c in SCHEMA)
GENERATED_COLUMNS: tuple[str, ...] = tuple(
    c.name for c in SCHEMA if c.origin in ("generated", "derived")
)
OVERLAY_COLUMNS: tuple[str, ...] = tuple(
    c.name for c in SCHEMA if c.origin in ("external", "reviewed")
)

DECISION_LOG_COLUMNS: tuple[str, ...] = (
    "variant",
    "column",
    "previous_value",
    "accepted_value",
    "evidence",
    "reason",
    "review_status",
    "review_date",
    "affects_downstream_results",
)

#: Allowed values for ``review_status``. "supervisor_approved" is deliberately
#: absent: no decision may be recorded as approved by a supervisor unless that
#: approval is documented elsewhere and the vocabulary is extended knowingly.
REVIEW_STATUS_VALUES = ("accepted", "provisional", "pending_second_review")


class AnnotationError(ValueError):
    """Raised when the annotation tables violate the schema or the merge rules."""


def _read(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix == ".tsv" else ";"
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False,
                       encoding="utf-8-sig", skiprows=_leading_comment_lines(path))


def load_overlay(path: Path) -> pd.DataFrame:
    """Read the reviewed-annotation overlay and check its shape."""
    if not Path(path).exists():
        raise AnnotationError(
            f"reviewed annotation overlay not found at {path}. Stage 04 cannot emit a "
            "complete dictionary without it, and emitting an incomplete one would "
            "silently drop the functional-evidence columns that stages 08 and 09 "
            "require. Copy data/reviewed/variant_annotations_reviewed.example.csv "
            "and fill it in, or set MEFV_ANNOTATION_OVERLAY."
        )
    df = _read(Path(path))
    expected = (KEY, *OVERLAY_COLUMNS)
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise AnnotationError(f"overlay {path} is missing column(s): {missing}")
    extra = [c for c in df.columns if c not in expected]
    if extra:
        raise AnnotationError(f"overlay {path} has unexpected column(s): {extra}")
    if df[KEY].duplicated().any():
        dupes = sorted(df.loc[df[KEY].duplicated(), KEY])
        raise AnnotationError(f"overlay {path} has duplicate variants: {dupes}")
    return df


def load_decision_log(path: Path) -> pd.DataFrame:
    """Read the annotation decision log, or return an empty log if absent."""
    if not Path(path).exists():
        return pd.DataFrame(columns=list(DECISION_LOG_COLUMNS))
    df = _read(Path(path))
    missing = [c for c in DECISION_LOG_COLUMNS if c not in df.columns]
    if missing:
        raise AnnotationError(f"decision log {path} is missing column(s): {missing}")
    bad = sorted(set(df.review_status) - set(REVIEW_STATUS_VALUES))
    if bad:
        raise AnnotationError(
            f"decision log {path} uses unrecognised review_status value(s): {bad}. "
            f"Allowed: {list(REVIEW_STATUS_VALUES)}"
        )
    return df


def check_reviewed_deviations(merged: pd.DataFrame, log: pd.DataFrame) -> list[str]:
    """Return a message for every undocumented reviewed deviation.

    A reviewed column may differ from the external column it was derived from
    only when the decision log records that variant and column.
    """
    documented = {(r.variant, r.column) for r in log.itertuples()}
    problems: list[str] = []
    for col in SCHEMA:
        if col.origin != "reviewed" or not col.derived_from:
            continue
        for row in merged.itertuples():
            value = getattr(row, col.name)
            source = getattr(row, col.derived_from)
            if value != source and (getattr(row, KEY), col.name) not in documented:
                problems.append(
                    f"{getattr(row, KEY)}: reviewed '{col.name}' is {value!r} but its source "
                    f"'{col.derived_from}' is {source!r}, and the decision log has no entry "
                    "for this variant and column"
                )
    return problems


def merge(generated: pd.DataFrame, overlay: pd.DataFrame, log: pd.DataFrame) -> pd.DataFrame:
    """Join generated annotations with the reviewed overlay and validate the result."""
    gen_missing = [c for c in GENERATED_COLUMNS if c not in generated.columns]
    if gen_missing:
        raise AnnotationError(f"generated annotations are missing column(s): {gen_missing}")

    only_generated = sorted(set(generated[KEY]) - set(overlay[KEY]))
    if only_generated:
        raise AnnotationError(
            f"{len(only_generated)} variant(s) have no row in the reviewed overlay: "
            f"{only_generated}. Add them to the overlay, leaving reviewed fields empty "
            "where no evidence exists, rather than letting the merge fill them silently."
        )
    only_overlay = sorted(set(overlay[KEY]) - set(generated[KEY]))
    if only_overlay:
        raise AnnotationError(
            f"the reviewed overlay contains variant(s) absent from the generated "
            f"annotations: {only_overlay}"
        )

    merged = generated.merge(overlay, on=KEY, how="left", validate="one_to_one")
    merged = merged[list(COLUMNS)]

    problems = check_reviewed_deviations(merged, log)
    if problems:
        raise AnnotationError(
            "undocumented deviation(s) between reviewed and external annotation:\n  "
            + "\n  ".join(problems)
        )
    return merged


def validate(dictionary: pd.DataFrame) -> None:
    """Check a finished dictionary against the schema. Raises on any violation."""
    if list(dictionary.columns) != list(COLUMNS):
        raise AnnotationError(
            f"dictionary columns do not match the schema.\n  expected: {list(COLUMNS)}\n"
            f"  found:    {list(dictionary.columns)}"
        )
    if dictionary[KEY].duplicated().any():
        raise AnnotationError("dictionary has duplicate variants")
    for col in ("variant", "hgvs_c", "hgvs_full", "tier"):
        empty = dictionary.loc[dictionary[col].astype(str).str.len() == 0, KEY].tolist()
        if empty:
            raise AnnotationError(f"column '{col}' must not be empty; empty for {empty}")
    bad_tier = sorted(set(dictionary.tier) - {"P", "LP", "VUS", "LB", "B", "unclassified"})
    if bad_tier:
        raise AnnotationError(f"unrecognised tier value(s): {bad_tier}")
    wrong = dictionary.loc[
        dictionary.hgvs_full != "NM_000243.3:" + dictionary.hgvs_c, KEY
    ].tolist()
    if wrong:
        raise AnnotationError(f"hgvs_full is inconsistent with hgvs_c for {wrong}")
