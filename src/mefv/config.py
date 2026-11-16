"""File paths, seeds and the interface to the private study parameters.

This module holds no study-specific value. Everything that describes a
particular study - its periods, the genomic regions every individual was
tested for, the assays used, the variants observed and how they are grouped -
is read from a local parameters file that is never committed
(``study_parameters.json``; structure in ``study_parameters.example.json``).
Accessors validate that file and raise a ``StudyParameterError`` naming the
missing or malformed entry. There are no fallback values: a missing setting
stops the pipeline rather than silently changing a result.

Paths can be overridden with environment variables, which is how the test
suite and the synthetic-data demonstration redirect inputs and outputs away
from any real dataset.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

NAVY, TEAL, AMBER, GREY = "#1F3864", "#1D9E75", "#BA7517", "#888780"


# ---------------------------------------------------------------- paths -----
ROOT = Path(__file__).resolve().parents[2]

# Stage 01 inputs: the laboratory export, as received. The default file names
# are generic; point the variables at the real export.
RAW_DIR = Path(os.environ.get("MEFV_RAW", ROOT / "data/raw"))
RAW_AUG = Path(os.environ.get("MEFV_RAW_AUG", RAW_DIR / "method.csv"))
RAW_VARIANTS = Path(os.environ.get("MEFV_RAW_VARIANTS", RAW_DIR / "variants.csv"))
RAW_PATIENTS = Path(os.environ.get("MEFV_RAW_PATIENTS", RAW_DIR / "patients.csv"))

# Stage 01 outputs, which are the inputs of every later stage.
PROCESSED_DIR = Path(os.environ.get("MEFV_PROCESSED", ROOT / "data/processed"))
PATIENTS = Path(os.environ.get("MEFV_PATIENTS", PROCESSED_DIR / "patients.csv"))
VARIANTS = Path(os.environ.get("MEFV_VARIANTS", PROCESSED_DIR / "variants_long.csv"))
DICTIONARY = Path(os.environ.get("MEFV_DICT", ROOT / "data/external/variant_dictionary.csv"))

# Stage 02: generated annotations, and the reviewed overlay it must never write.
ANNOTATIONS_GENERATED = Path(
    os.environ.get("MEFV_ANNOTATIONS_GENERATED",
                   ROOT / "data/external/variant_annotations_generated.csv")
)
POPULATION_FREQUENCIES = Path(
    os.environ.get("MEFV_POPULATION_FREQUENCIES",
                   ROOT / "data/external/population_frequencies.csv")
)
SOURCE_ANNOTATIONS = Path(
    os.environ.get("MEFV_SOURCE_ANNOTATIONS",
                   ROOT / "data/external/source_annotations.csv")
)
ANNOTATION_OVERLAY = Path(
    os.environ.get("MEFV_ANNOTATION_OVERLAY",
                   ROOT / "data/reviewed/variant_annotations_reviewed.csv")
)
ANNOTATION_DECISION_LOG = Path(
    os.environ.get("MEFV_ANNOTATION_DECISION_LOG",
                   ROOT / "data/reviewed/annotation_decision_log.csv")
)

TABLES = Path(os.environ.get("MEFV_TABLES", ROOT / "results/tables"))
FIGURES = Path(os.environ.get("MEFV_FIGURES", ROOT / "results/figures"))
CLASSIFICATION = Path(os.environ.get("MEFV_ACMG", ROOT / "results/classification"))

#: Execution record written by scripts/run_pipeline.py. It sits outside the
#: output directories so that one run's manifest is never collected as an
#: output of the next.
RUN_MANIFEST = Path(os.environ.get("MEFV_RUN_MANIFEST", TABLES.parent / "run_manifest.json"))

# ------------------------------------------------------- random number seeds -
# Every stochastic step uses an explicit seed so that repeated runs of the
# same input produce byte-identical output.
SEED_BOOTSTRAP = 20260803   # bootstrap confidence intervals for D'
SEED_FIGURES = 7            # jitter in scatter plots

BOOTSTRAP_REPLICATES = int(os.environ.get("MEFV_BOOTSTRAP_REPLICATES", "2000"))

#: Fixed timestamp for file formats that embed a creation date.
#: Matplotlib writes /CreationDate into every PDF, which makes two otherwise
#: identical runs produce different bytes. Setting SOURCE_DATE_EPOCH - the
#: reproducible-builds convention, which matplotlib honours - pins it, so a
#: digest comparison between runs is meaningful for PDF output too.
SOURCE_DATE_EPOCH = os.environ.get("SOURCE_DATE_EPOCH", "1767225600")  # 2026-01-01T00:00:00Z


def pin_build_timestamp() -> None:
    """Pin the embedded creation date used by figure output formats."""
    os.environ.setdefault("SOURCE_DATE_EPOCH", SOURCE_DATE_EPOCH)


# ------------------------------------------------------ study parameters ----
STUDY_PARAMETERS = Path(
    os.environ.get("MEFV_STUDY_PARAMETERS", ROOT / "study_parameters.json")
)


#: Groupings that variant_groups() must supply, and what each one means.
VARIANT_GROUP_KEYS = {
    "synonymous": "synonymous variants; excluded from allele burden",
    "block_synonymous": "synonymous variants analysed together as one block",
    "benign_common": "common benign variants excluded from allele burden",
    "unverified": "variants whose original report could not be traced",
}

#: Keys that linkage_plan() must supply.
LINKAGE_PLAN_KEYS = {
    "primary_pair": "the two variants of the main linkage analysis, as [A, B]",
    "haplotype_blocks": "named blocks of variants analysed together",
    "cross_locus_pairs": "{anchor variant: [partners]} for the block-to-locus comparison",
}


class StudyParameterError(ValueError):
    """A required study parameter is missing or malformed."""


#: Entries that study_design() requires, and what each one means.
STUDY_DESIGN_KEYS = {
    "periods": "{'primary': PERIOD, 'sensitivity': PERIOD}; PERIOD is "
               "{'label': str, 'first_year': int, 'last_year': int}",
    "panel_regions": "ordered list of the regions (e.g. exons) tested in every individual",
    "assays": "{'panel': label of the assay restricted to the panel regions, "
              "'extended': label of the assay covering further regions}",
    "assay_inference": "null, or {'until_year': int, 'assay': one of the assay labels, "
                       "'source_label': str}: blank assay fields up to that year are "
                       "filled from a documented laboratory statement",
}


@dataclass(frozen=True)
class Period:
    label: str
    first_year: int
    last_year: int

    def contains(self, year: int) -> bool:
        return self.first_year <= year <= self.last_year


@dataclass(frozen=True)
class StudyDesign:
    """Validated study design, read from the private parameters file."""

    primary: Period
    sensitivity: Period
    panel_regions: tuple[str, ...]
    panel_assay: str
    extended_assay: str
    inference_until_year: int | None
    inferred_assay: str | None
    inferred_source_label: str | None

    @property
    def assays(self) -> tuple[str, str]:
        return (self.panel_assay, self.extended_assay)

    def period_of(self, year: int) -> str:
        """Label of the period containing ``year``; an uncovered year is an error."""
        for p in (self.primary, self.sensitivity):
            if p.contains(int(year)):
                return p.label
        raise StudyParameterError(
            f"year {year} is outside both study periods "
            f"({self.primary.first_year}-{self.primary.last_year}, "
            f"{self.sensitivity.first_year}-{self.sensitivity.last_year})"
        )

    def region_rank(self, region: str) -> int:
        """Position of a region in the configured panel order (unknown regions last)."""
        try:
            return self.panel_regions.index(str(region))
        except ValueError:
            return len(self.panel_regions)


def _require(mapping: dict, key: str, where: str):
    if not isinstance(mapping, dict) or key not in mapping:
        raise StudyParameterError(
            f"'{where}{key}' is missing from {STUDY_PARAMETERS}. "
            "See study_parameters.example.json and docs/configuration.md."
        )
    return mapping[key]


def _period(raw: object, name: str) -> Period:
    where = f"study_design.periods.{name}."
    label = _require(raw, "label", where)
    first, last = _require(raw, "first_year", where), _require(raw, "last_year", where)
    if not isinstance(label, str) or not label:
        raise StudyParameterError(f"'{where}label' must be a non-empty string")
    if not all(isinstance(y, int) and not isinstance(y, bool) for y in (first, last)):
        raise StudyParameterError(f"'{where}first_year' and 'last_year' must be integers")
    if first > last:
        raise StudyParameterError(f"'{where}first_year' is after 'last_year'")
    return Period(label, first, last)


def parse_study_design(params: dict) -> StudyDesign:
    """Validate the ``study_design`` section of a parameters mapping."""
    raw = _require(params, "study_design", "")
    periods = _require(raw, "periods", "study_design.")
    primary = _period(_require(periods, "primary", "study_design.periods."), "primary")
    sensitivity = _period(_require(periods, "sensitivity", "study_design.periods."),
                          "sensitivity")
    if primary.label == sensitivity.label:
        raise StudyParameterError("the two study periods must have different labels")
    if not (primary.last_year < sensitivity.first_year
            or sensitivity.last_year < primary.first_year):
        raise StudyParameterError("the two study periods overlap")

    regions = _require(raw, "panel_regions", "study_design.")
    if (not isinstance(regions, list) or not regions
            or not all(isinstance(r, str) and r for r in regions)
            or len(set(regions)) != len(regions)):
        raise StudyParameterError(
            "'study_design.panel_regions' must be a non-empty list of distinct strings")

    assays = _require(raw, "assays", "study_design.")
    panel, extended = (_require(assays, k, "study_design.assays.") for k in ("panel", "extended"))
    if not all(isinstance(a, str) and a for a in (panel, extended)) or panel == extended:
        raise StudyParameterError(
            "'study_design.assays' needs two different, non-empty labels")

    inference = _require(raw, "assay_inference", "study_design.")
    if inference is None:
        until = assay = label = None
    else:
        until = _require(inference, "until_year", "study_design.assay_inference.")
        assay = _require(inference, "assay", "study_design.assay_inference.")
        label = _require(inference, "source_label", "study_design.assay_inference.")
        if not isinstance(until, int) or isinstance(until, bool):
            raise StudyParameterError(
                "'study_design.assay_inference.until_year' must be an integer")
        if assay not in (panel, extended):
            raise StudyParameterError(
                "'study_design.assay_inference.assay' must be one of the two assay labels")
        if not isinstance(label, str) or not label:
            raise StudyParameterError(
                "'study_design.assay_inference.source_label' must be a non-empty string")

    return StudyDesign(primary, sensitivity, tuple(regions), panel, extended,
                       until, assay, label)


def study_design() -> StudyDesign:
    """The study design, from the study-parameters file (validated, no defaults)."""
    return parse_study_design(load_study_parameters())


def variant_groups() -> dict[str, object]:
    """Cohort-specific variant groupings, from the study-parameters file."""
    params = load_study_parameters()
    groups = params.get("variant_groups")
    if groups is None:
        raise StudyParameterError(
            f"'variant_groups' is missing from {STUDY_PARAMETERS}. It must provide "
            f"{sorted(VARIANT_GROUP_KEYS)}. See study_parameters.example.json."
        )
    missing = sorted(set(VARIANT_GROUP_KEYS) - set(groups))
    if missing:
        raise StudyParameterError(f"'variant_groups' in {STUDY_PARAMETERS} is missing: {missing}")
    return {k: set(v) if isinstance(v, list) else v for k, v in groups.items()}


def linkage_plan() -> dict[str, object]:
    """Which variants the linkage stage analyses, from the study-parameters file."""
    params = load_study_parameters()
    plan = params.get("linkage_plan")
    if plan is None:
        raise StudyParameterError(
            f"'linkage_plan' is missing from {STUDY_PARAMETERS}. It must provide "
            f"{sorted(LINKAGE_PLAN_KEYS)}. See study_parameters.example.json."
        )
    missing = sorted(set(LINKAGE_PLAN_KEYS) - set(plan))
    if missing:
        raise StudyParameterError(f"'linkage_plan' in {STUDY_PARAMETERS} is missing: {missing}")
    return plan


def load_study_parameters() -> dict:
    """Return the local, non-public study parameters.

    Everything study-specific - the design, the expected number of records,
    the record-level corrections confirmed by the laboratory and the variants
    observed - is kept out of the repository in a local JSON file excluded by
    .gitignore; ``study_parameters.example.json`` documents its shape and
    ``docs/configuration.md`` each key.

    The function raises rather than returning defaults. A missing file must
    stop the pipeline, because silently substituting an empty override table
    would change the result without any visible sign.
    """
    return _load_parameters(str(STUDY_PARAMETERS))


def _load_parameters(path: str) -> dict:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"study parameters not found at {p}. "
            "Copy study_parameters.example.json to study_parameters.json and "
            "fill it in, or set MEFV_STUDY_PARAMETERS. See docs/configuration.md."
        )
    try:
        params = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise StudyParameterError(f"{p} is not valid JSON: {exc}") from exc
    if not isinstance(params, dict):
        raise StudyParameterError(f"{p} must contain a JSON object")
    return params
