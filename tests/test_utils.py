"""Tests for genotype encoding and the two-table loader."""
from __future__ import annotations

import pandas as pd
import pytest

from mefv.utils import DOSE, genotype_vector


def make_variants(rows: list[tuple[str, str, str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["patient_id", "variant", "zygosity"])


def test_genotype_vector_encodes_dosage():
    var = make_variants(
        [("P1", "SYN_P01", "Hom"), ("P2", "SYN_P01", "Het"), ("P3", "SYN_B01", "Het")]
    )
    assert genotype_vector(var, "SYN_P01", ["P1", "P2", "P3", "P4"]) == [2, 1, 0, 0]


def test_genotype_vector_preserves_the_requested_order():
    var = make_variants([("P1", "SYN_P01", "Het"), ("P2", "SYN_P01", "Hom")])
    assert genotype_vector(var, "SYN_P01", ["P2", "P1"]) == [2, 1]


def test_patients_without_a_row_are_homozygous_reference():
    """The denominator includes individuals in whom nothing was found."""
    var = make_variants([("P1", "SYN_P01", "Het")])
    assert genotype_vector(var, "SYN_P01", [f"P{i}" for i in range(1, 11)]) == [1] + [0] * 9


def test_unknown_variant_returns_all_zeros():
    var = make_variants([("P1", "SYN_P01", "Het")])
    assert genotype_vector(var, "NOT_OBSERVED", ["P1", "P2"]) == [0, 0]


def test_unrecognised_zygosity_warns_and_counts_one_copy():
    """Preserves the original behaviour, but no longer silently."""
    var = make_variants([("P1", "SYN_P01", "???")])
    with pytest.warns(UserWarning, match="unrecognised zygosity"):
        assert genotype_vector(var, "SYN_P01", ["P1"]) == [1]


def test_dose_mapping_is_the_documented_one():
    assert DOSE == {"Het": 1, "Hom": 2}
