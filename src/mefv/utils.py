"""Input/output helpers and genotype encoding.

The analysis uses a two-table design:

* ``patients``       one row per tested individual; this is the denominator
                     for every rate reported by the pipeline.
* ``variants_long``  one row per observed variant in one patient.

Individuals without an observed variant have no row in ``variants_long`` and
are therefore encoded as homozygous reference. This is only valid for regions
that were tested in every individual - the configured panel regions - which is
why analyses that rely on it restrict themselves to those regions; see
docs/limitations.md.
"""
from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from mefv.config import StudyDesign
from mefv.disclosure import text as _n
from mefv.validation import DataIntegrityError, parse_uint_column


def _leading_comment_lines(path: str | Path) -> int:
    """Count the ``#`` banner lines at the top of a file, if any.

    Synthetic data files carry a banner marking them as artificial. Only
    leading lines are considered, so a ``#`` inside the data is untouched.
    """
    n = 0
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                n += 1
            else:
                break
    return n


def read_csv(path: str | Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep=";", encoding="utf-8-sig", dtype=str,
                    keep_default_na=False, engine="python",
                    skiprows=_leading_comment_lines(path))
    return d.loc[:, [c for c in d.columns if not c.startswith("Unnamed")]]


def require_released(pat_path: str | Path, *, allow_unreported: bool = False) -> None:
    """Refuse analysis tables unless their own validation report released them.

    Fails closed. The tables are accepted only if a validation report sits next
    to them, parses as a JSON object, and states ``"outputs_released": true``
    (the JSON value, not a string or a number). A missing report, an unreadable
    or malformed one, one without that key, or one with any other value is a
    refusal, as is ``false``. A blocked, failed or diagnostic-only stage-01 run,
    or one with an acknowledged but undisposed finding, therefore cannot feed
    later stages, even if tables from an earlier run are still in the folder.

    ``allow_unreported=True`` is for library and test use on tables that were
    not made by stage 01 and so have no report. No pipeline stage sets it. It
    skips only the "no report" case: a report that exists is always judged.
    """
    report_path = Path(pat_path).parent / "validation_report.json"
    if not report_path.is_file():
        if allow_unreported:
            return
        raise DataIntegrityError(
            "no validation report sits next to the analysis tables, so they are not released; "
            "refusing to analyse them. Run stage 01. See docs/validation.md."
        )
    import json

    try:
        report = json.loads(report_path.read_text())
    except (ValueError, OSError):
        raise DataIntegrityError(
            "the validation report next to the analysis tables is unreadable, so they are not "
            "released; refusing to analyse them"
        ) from None
    if not isinstance(report, dict) or report.get("outputs_released") is not True:
        raise DataIntegrityError(
            "the validation report does not state that the final outputs were released (blocking "
            "findings without an owner disposition, a failed validation, or no release recorded): "
            "the tables are not released; refusing to analyse them. See docs/validation.md."
        )


def load(pat_path: str | Path, var_path: str | Path,
         design: StudyDesign | None = None, *,
         allow_unreported: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the two analysis tables and label each row with its study period.

    ``design`` defaults to the one in the study-parameters file. A year outside
    both configured periods raises rather than being assigned silently.
    """
    if design is None:
        from mefv.config import study_design
        design = study_design()
    require_released(pat_path, allow_unreported=allow_unreported)
    pat = read_csv(pat_path)
    var = read_csv(var_path)
    pat = pat[pat.patient_id.str.strip() != ""].copy()
    # Integer parsing reports counts, never the offending value.
    pat["year"] = parse_uint_column(pat.year_requested, "patients.year_requested")
    var["year"] = parse_uint_column(var.year_requested, "variants.year_requested")
    # blank = declared count unknown (individuals restored from the method table)
    pat["n_variants"] = parse_uint_column(pat.n_variants, "patients.n_variants",
                                          allow_blank=True)
    for d in (pat, var):
        d["era"] = [design.period_of(y) for y in d.year]
    return pat, var


DOSE = {"Het": 1, "Hom": 2}


def allele_dose(zygosity: pd.Series, what: str = "variants.zygosity") -> pd.Series:
    """Copies of the variant allele (1 or 2) for each zygosity label.

    A label other than exactly ``Het`` or ``Hom`` raises. It is not counted as
    one copy and not dropped: which allele count it stands for is a decision
    for the data owner, and the error reports how many entries are affected
    without echoing them.
    """
    dose = zygosity.map(DOSE)
    if dose.isna().any():
        raise DataIntegrityError(
            f"{what}: {_n(int(dose.isna().sum()))} entr(ies) are not 'Het' or 'Hom'; "
            "run the stage 01 validation (rule V4).")
    return dose.astype(int)


def genotype_vector(
    var: pd.DataFrame, variant: str, patient_ids: Sequence[str]
) -> list[int]:
    """Copies of the variant allele per patient, in the order of patient_ids.

    Returns 0, 1 or 2 for every identifier in ``patient_ids``. A patient with
    no row for this variant is encoded as 0, which is the documented
    data-generation rule (see the module docstring and docs/limitations.md);
    it is valid only for individuals tested for the variant's region.

    Exactly one genotype per individual and variant is required. If an
    individual has more than one row for this variant, whether the rows repeat
    the same genotype or contradict each other, the function raises
    :class:`DataIntegrityError`: choosing, merging or summing them is a decision
    for the data owner (validation rules V5, V6a, V6b), and a vector built from
    an arbitrary one of them would silently misstate the genotype. A repeated
    identifier in ``patient_ids`` raises for the same reason. An unrecognised
    zygosity label raises too.
    """
    ids = list(patient_ids)
    if len(set(ids)) != len(ids):
        raise DataIntegrityError(
            "genotype_vector: an individual is listed more than once in the requested "
            "identifiers; a genotype vector needs one position per individual."
        )
    sub = var[var.variant == variant]
    repeated = sub.patient_id.duplicated(keep=False)
    if repeated.any():
        raise DataIntegrityError(
            f"genotype_vector: {_n(int(sub.patient_id[repeated].nunique()))} individual(s) have "
            "more than one row for the same variant (validation rules V5/V6a/V6b). The "
            "genotype of such an individual is not defined until the data owner decides how "
            "it is to be derived; nothing was chosen, merged or dropped. See "
            "docs/validation.md."
        )
    d = dict(zip(sub.patient_id, allele_dose(sub.zygosity)))
    return [int(d.get(p, 0)) for p in ids]
