"""Tests for the annotation schema and the generated/reviewed merge rules.

These tests use small fixtures rather than the study dictionary, so they run
from a clean clone with no private data present.
"""
from __future__ import annotations

import pandas as pd
import pytest

from mefv import annotation as A


# ------------------------------------------------------------- fixtures ----
def generated_frame(variants=("V1", "V2")) -> pd.DataFrame:
    rows = []
    for i, v in enumerate(variants, start=1):
        rows.append({
            "variant": v,
            "infevers_name": f"p.Name{i}",
            "hgvs_full": f"NM_000243.3:c.{i}A>G",
            "hgvs_c": f"c.{i}A>G",
            "hgvs_p": f"p.(Xaa{i}Yaa)",
            "exon_infevers": "10",
            "rsid": f"rs{i}",
            "infevers_class": "Pathogenic",
            "previous_class": "",
            "comment": "",
            "tier": "P",
            "reclassified_2025": "False",
        })
    return pd.DataFrame(rows)[list(A.GENERATED_COLUMNS)]


def overlay_frame(variants=("V1", "V2"), speckseq=None) -> pd.DataFrame:
    rows = []
    for v in variants:
        label = "Not hypermorphic"
        rows.append({
            "variant": v,
            "speckseq": (speckseq or {}).get(v, label),
            "other_functional": "",
            "speckseq_assay_level_infevers": label,
            "speckseq_integrated_S7": "BS3",
            "bronnec2026_revised_S9": "LB",
        })
    return pd.DataFrame(rows)[[A.KEY, *A.OVERLAY_COLUMNS]]


def empty_log() -> pd.DataFrame:
    return pd.DataFrame(columns=list(A.DECISION_LOG_COLUMNS))


def log_with(variant: str, column: str) -> pd.DataFrame:
    return pd.DataFrame([{
        "variant": variant, "column": column,
        "previous_value": "a", "accepted_value": "b",
        "evidence": "citation", "reason": "reason",
        "review_status": "accepted", "review_date": "2026-01-01",
        "affects_downstream_results": "no",
    }], columns=list(A.DECISION_LOG_COLUMNS))


# ---------------------------------------------------------------- schema ---
def test_schema_has_seventeen_columns_split_twelve_and_five():
    assert len(A.COLUMNS) == 17
    assert len(A.GENERATED_COLUMNS) == 12
    assert len(A.OVERLAY_COLUMNS) == 5
    assert set(A.GENERATED_COLUMNS) | set(A.OVERLAY_COLUMNS) == set(A.COLUMNS)


def test_every_column_declares_a_known_origin():
    assert {c.origin for c in A.SCHEMA} <= {"generated", "external", "reviewed", "derived"}


def test_derived_columns_name_their_source():
    for c in A.SCHEMA:
        if c.origin in ("derived", "reviewed"):
            assert c.derived_from, f"{c.name} must record what it is derived from"
            assert c.derived_from in A.COLUMNS


def test_columns_required_by_downstream_stages_are_present():
    # 08 reads tier and speckseq; 09 reads the HGVS, database and evidence fields.
    for name in ("tier", "speckseq", "other_functional", "hgvs_c", "hgvs_p",
                 "exon_infevers", "rsid", "infevers_class", "previous_class", "comment"):
        assert name in A.COLUMNS


def test_review_status_vocabulary_excludes_supervisor_approval():
    assert "supervisor_approved" not in A.REVIEW_STATUS_VALUES


# ----------------------------------------------------------------- merge ---
def test_merge_produces_the_full_schema_in_order():
    merged = A.merge(generated_frame(), overlay_frame(), empty_log())
    assert list(merged.columns) == list(A.COLUMNS)
    assert len(merged) == 2
    A.validate(merged)


def test_merge_is_deterministic():
    a = A.merge(generated_frame(), overlay_frame(), empty_log())
    b = A.merge(generated_frame(), overlay_frame(), empty_log())
    assert a.equals(b)


def test_variant_missing_from_the_overlay_is_an_error():
    with pytest.raises(A.AnnotationError, match="no row in the reviewed overlay"):
        A.merge(generated_frame(("V1", "V2")), overlay_frame(("V1",)), empty_log())


def test_overlay_may_not_introduce_variants():
    with pytest.raises(A.AnnotationError, match="absent from the generated"):
        A.merge(generated_frame(("V1",)), overlay_frame(("V1", "V2")), empty_log())


def test_generated_table_missing_a_column_is_an_error():
    gen = generated_frame().drop(columns=["tier"])
    with pytest.raises(A.AnnotationError, match="missing column"):
        A.merge(gen, overlay_frame(), empty_log())


# --------------------------------------- reviewed deviations and the log ---
def test_undocumented_reviewed_deviation_is_rejected():
    overlay = overlay_frame(speckseq={"V1": "Uncertain"})
    with pytest.raises(A.AnnotationError, match="undocumented deviation"):
        A.merge(generated_frame(), overlay, empty_log())


def test_documented_reviewed_deviation_is_accepted():
    overlay = overlay_frame(speckseq={"V1": "Uncertain"})
    merged = A.merge(generated_frame(), overlay, log_with("V1", "speckseq"))
    assert merged.loc[merged.variant == "V1", "speckseq"].item() == "Uncertain"


def test_a_log_entry_for_another_column_does_not_excuse_the_deviation():
    overlay = overlay_frame(speckseq={"V1": "Uncertain"})
    with pytest.raises(A.AnnotationError, match="undocumented deviation"):
        A.merge(generated_frame(), overlay, log_with("V1", "other_functional"))


def test_regeneration_cannot_overwrite_a_reviewed_value():
    """Re-running the merge leaves the reviewed value untouched."""
    overlay = overlay_frame(speckseq={"V1": "Uncertain"})
    log = log_with("V1", "speckseq")
    first = A.merge(generated_frame(), overlay, log)
    second = A.merge(generated_frame(), overlay, log)
    assert first.loc[first.variant == "V1", "speckseq"].item() == "Uncertain"
    assert second.equals(first)


def test_decision_log_rejects_an_unknown_review_status(tmp_path):
    p = tmp_path / "log.csv"
    bad = log_with("V1", "speckseq")
    bad.loc[0, "review_status"] = "supervisor_approved"
    bad.to_csv(p, sep=";", index=False, encoding="utf-8-sig")
    with pytest.raises(A.AnnotationError, match="unrecognised review_status"):
        A.load_decision_log(p)


def test_absent_decision_log_reads_as_empty(tmp_path):
    log = A.load_decision_log(tmp_path / "missing.csv")
    assert log.empty
    assert list(log.columns) == list(A.DECISION_LOG_COLUMNS)


# ------------------------------------------------------------- overlay io --
def test_missing_overlay_raises_a_actionable_error(tmp_path):
    with pytest.raises(A.AnnotationError, match="Copy data/reviewed"):
        A.load_overlay(tmp_path / "missing.csv")


def test_overlay_with_wrong_columns_is_rejected(tmp_path):
    p = tmp_path / "overlay.csv"
    overlay_frame().drop(columns=["speckseq"]).to_csv(p, sep=";", index=False, encoding="utf-8-sig")
    with pytest.raises(A.AnnotationError, match="missing column"):
        A.load_overlay(p)


def test_overlay_with_duplicate_variants_is_rejected(tmp_path):
    p = tmp_path / "overlay.csv"
    pd.concat([overlay_frame(("V1",)), overlay_frame(("V1",))]).to_csv(
        p, sep=";", index=False, encoding="utf-8-sig")
    with pytest.raises(A.AnnotationError, match="duplicate variants"):
        A.load_overlay(p)


# -------------------------------------------------------------- validate ---
def test_validate_rejects_a_wrong_column_order():
    merged = A.merge(generated_frame(), overlay_frame(), empty_log())
    shuffled = merged[list(reversed(A.COLUMNS))]
    with pytest.raises(A.AnnotationError, match="columns do not match"):
        A.validate(shuffled)


def test_validate_rejects_an_empty_required_field():
    merged = A.merge(generated_frame(), overlay_frame(), empty_log())
    merged.loc[0, "tier"] = ""
    with pytest.raises(A.AnnotationError, match="must not be empty"):
        A.validate(merged)


def test_validate_rejects_an_unknown_tier():
    merged = A.merge(generated_frame(), overlay_frame(), empty_log())
    merged.loc[0, "tier"] = "probably pathogenic"
    with pytest.raises(A.AnnotationError, match="unrecognised tier"):
        A.validate(merged)


def test_validate_rejects_hgvs_full_inconsistent_with_hgvs_c():
    merged = A.merge(generated_frame(), overlay_frame(), empty_log())
    merged.loc[0, "hgvs_full"] = "NM_000243.3:c.999A>G"
    with pytest.raises(A.AnnotationError, match="hgvs_full is inconsistent"):
        A.validate(merged)
