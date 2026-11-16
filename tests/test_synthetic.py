"""Tests for the synthetic data generator.

These run from a clean clone: nothing here touches study data.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from mefv import annotation as A

SPEC = importlib.util.spec_from_file_location(
    "make_synthetic_data",
    Path(__file__).resolve().parents[1] / "scripts" / "make_synthetic_data.py",
)
syn = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(syn)


def test_generation_is_reproducible_for_a_given_seed():
    a = syn.generate(n=120, seed=7)
    b = syn.generate(n=120, seed=7)
    for key in ("patients", "variants_long", "source", "overlay"):
        assert a[key].equals(b[key]), key


def test_a_different_seed_gives_different_data():
    a = syn.generate(n=120, seed=7)
    b = syn.generate(n=120, seed=8)
    assert not a["variants_long"].equals(b["variants_long"])


def test_identifiers_are_marked_synthetic_and_cannot_collide():
    d = syn.generate(n=50, seed=1)
    assert d["patients"].patient_id.str.fullmatch(r"SYNTH_\d{4}").all()


def test_provenance_states_that_the_data_are_artificial():
    d = syn.generate(n=50, seed=1)
    p = d["provenance"]
    assert p["synthetic"] is True
    assert "not sampled from" in p["statement"].lower()
    assert p["seed"] == 1


def test_written_files_carry_a_synthetic_banner(tmp_path):
    syn.write(tmp_path, syn.generate(n=40, seed=3))
    for rel in ("raw/synthetic_method.csv", "raw/synthetic_variants.csv",
                "raw/synthetic_patients.csv",
                "external/source_annotations.csv",
                "reviewed/variant_annotations_reviewed.csv"):
        head = (tmp_path / rel).read_text(encoding="utf-8-sig").splitlines()[0]
        assert head.startswith("#") and "SYNTHETIC" in head, rel
    assert json.loads((tmp_path / "PROVENANCE.json").read_text())["synthetic"] is True


def test_written_files_are_readable_by_the_pipeline_reader(tmp_path):
    """The banner must not break the loader."""
    from mefv.utils import read_csv

    syn.write(tmp_path, syn.generate(n=40, seed=3))
    pat = read_csv(tmp_path / "raw/synthetic_patients.csv")
    var = read_csv(tmp_path / "raw/synthetic_variants.csv")
    assert len(pat) == 39  # one record is deliberately left out, for stage 01 to restore
    assert "patient_id" in pat.columns and "patient_id" in var.columns
    assert set(var.zygosity) <= {"Het", "Hom"}


def test_overlay_matches_the_annotation_schema(tmp_path):
    syn.write(tmp_path, syn.generate(n=40, seed=3))
    overlay = A.load_overlay(tmp_path / "reviewed/variant_annotations_reviewed.csv")
    assert list(overlay.columns) == [A.KEY, *A.OVERLAY_COLUMNS]


def test_source_and_overlay_cover_the_same_variants():
    d = syn.generate(n=40, seed=3)
    assert set(d["source"].variant) == set(d["overlay"].variant)


def test_the_linked_pair_shows_an_association():
    """The generator imposes linkage so the linkage stage has a signal."""
    from mefv import genotype_vector, ld_stats

    d = syn.generate(n=800, seed=11)
    ids = list(d["patients"].patient_id)
    a, b = syn.LD_PAIR
    out = ld_stats(genotype_vector(d["variants_long"], a, ids),
                   genotype_vector(d["variants_long"], b, ids))
    assert out["Dprime"] > 0.5


def test_synthetic_frequencies_are_the_declared_ones_not_a_real_cohort():
    """Observed frequencies track the invented parameters, within sampling error."""
    d = syn.generate(n=4000, seed=5)
    var = d["variants_long"]
    n = len(d["patients"])
    for name, _, declared, _ in syn.CATALOGUE:
        if name == syn.LD_PAIR[1] or name in syn.EXTENDED_ONLY:
            continue  # thinned on purpose, or observable only with one assay
        sub = var[var.variant == name]
        observed = (sub.zygosity.map({"Het": 1, "Hom": 2}).sum()) / (2 * n)
        assert observed == pytest.approx(declared, abs=0.02), name


def test_the_complex_allele_pair_co_occurs_perfectly():
    """Every carrier of one member carries the other, so stage 04 has a candidate."""
    d = syn.generate(n=400, seed=13)
    var = d["variants_long"]
    a, b = syn.COMPLEX_PAIR
    carriers_a = set(var[var.variant == a].patient_id)
    carriers_b = set(var[var.variant == b].patient_id)
    assert len(carriers_a) >= 2
    assert carriers_a == carriers_b


def test_raw_export_contains_one_instance_of_each_known_defect():
    """Stage 01 has five correction paths; the synthetic export triggers each."""
    d = syn.generate(n=200, seed=17)
    # a record present in the method table but missing from the patient table
    assert d["restored_id"] not in set(d["raw_patients"].patient_id)
    assert d["restored_id"].lower() in {x.lower() for x in d["raw_method"].Patient_ID}
    # a late-era record with a blank method, resolved by an override
    m = d["raw_method"].set_index(d["raw_method"].Patient_ID.str.upper())
    assert m.loc[d["override_id"], "Method"] == ""
    assert d["override_id"] in d["parameters"]["assay_overrides"]
    # a record whose exported assay is wrong, also corrected by an override
    confirmed = d["parameters"]["assay_overrides"][d["corrected_id"]]
    assert m.loc[d["corrected_id"], "Method"] not in ("", confirmed)
    # a mistyped label, resolved by a nomenclature correction
    assert (d["raw_variants"].variant == syn.TYPO[0]).sum() == 1
    assert d["parameters"]["nomenclature_fixes"] == {syn.TYPO[0]: syn.TYPO[1]}
    # blank assay fields in early records, for the documented inference rule
    rule = syn.SYNTHETIC_DESIGN["assay_inference"]
    early_blank = (d["raw_method"].Method == "") & (d["raw_method"].YEAR <= rule["until_year"])
    assert early_blank.sum() > 0


def test_extended_only_variants_occur_only_with_the_extended_assay():
    d = syn.generate(n=600, seed=19)
    var = d["variants_long"]
    sub = var[var.variant.isin(syn.EXTENDED_ONLY)]
    assert len(sub) > 0
    assert set(sub.assay) == {syn.SYNTHETIC_DESIGN["assays"]["extended"]}
