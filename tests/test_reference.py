"""Tests for population reference-data validation (mefv.reference)."""
from __future__ import annotations

import pytest

from mefv.reference import COLUMNS, ReferenceDataError, load_population_frequencies, read_provenance

HEADER = (
    "# source      : Example database\n"
    "#               (continued on a second line)\n"
    "# version     : 9.9\n"
    "# retrieved   : 2000-01-01\n"
    "# open        : nothing outstanding\n"
)


def row(variant: str = "EX_1", **override) -> str:
    vals = {c: "0.01" for c in COLUMNS}
    vals.update(variant=variant, hom="3", cv="Benign", stars="1", revel="", cadd="")
    vals.update(override)
    return ";".join(vals[c] for c in COLUMNS)


def write(tmp_path, header=HEADER, rows=None, columns=COLUMNS):
    p = tmp_path / "ref.csv"
    body = ";".join(columns) + "\n" + "\n".join(rows or [row()]) + "\n"
    p.write_text(header + body, encoding="utf-8-sig")
    return p


def test_a_valid_file_loads_with_its_provenance(tmp_path):
    df, prov = load_population_frequencies(write(tmp_path, rows=[row("EX_1"), row("EX_2")]))
    assert list(df.variant) == ["EX_1", "EX_2"]
    assert df.mid.dtype.kind == "f" and df.revel.isna().all()
    assert prov["version"] == "9.9"
    assert prov["source"] == "Example database (continued on a second line)"


def test_missing_file_is_a_clear_error(tmp_path):
    with pytest.raises(ReferenceDataError, match="not found"):
        load_population_frequencies(tmp_path / "absent.csv")


@pytest.mark.parametrize("drop", ["source", "version", "retrieved"])
def test_each_provenance_line_is_required(tmp_path, drop):
    header = "".join(line + "\n" for line in HEADER.splitlines()
                     if not line.startswith(f"# {drop}"))
    with pytest.raises(ReferenceDataError, match=drop):
        load_population_frequencies(write(tmp_path, header=header))


def test_a_different_version_is_not_substituted(tmp_path):
    with pytest.raises(ReferenceDataError, match="expect"):
        load_population_frequencies(write(tmp_path), expected_version="1.0")
    load_population_frequencies(write(tmp_path), expected_version="9.9")


def test_wrong_columns_are_rejected(tmp_path):
    cols = tuple(c for c in COLUMNS if c != "grpmax")
    with pytest.raises(ReferenceDataError, match="columns"):
        load_population_frequencies(write(tmp_path, columns=cols, rows=["x"]))


def test_duplicated_variants_are_rejected(tmp_path):
    with pytest.raises(ReferenceDataError, match="duplicated"):
        load_population_frequencies(write(tmp_path, rows=[row("EX_1"), row("EX_1")]))


def test_non_numeric_values_are_rejected(tmp_path):
    with pytest.raises(ReferenceDataError, match="non-numeric"):
        load_population_frequencies(write(tmp_path, rows=[row(nfe="n/a")]))


def test_frequencies_outside_the_unit_interval_are_rejected(tmp_path):
    with pytest.raises(ReferenceDataError, match="outside"):
        load_population_frequencies(write(tmp_path, rows=[row(grpmax="1.5")]))


def test_provenance_reader_ignores_free_text_blocks(tmp_path):
    header = HEADER + "#\n# columns\n#   mid  frequency : not a key\n"
    prov = read_provenance(write(tmp_path, header=header))
    assert set(prov) == {"source", "version", "retrieved", "open"}
