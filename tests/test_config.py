"""Tests for the configuration boundary between public code and private study settings."""
from __future__ import annotations

import ast
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from mefv import config
from mefv.config import StudyParameterError, parse_study_design
from mefv.utils import load

ROOT = Path(__file__).resolve().parents[1]

DESIGN = {
    "periods": {
        "primary": {"label": "P-A", "first_year": 100, "last_year": 109},
        "sensitivity": {"label": "P-B", "first_year": 110, "last_year": 119},
    },
    "panel_regions": ["X1", "X2"],
    "assays": {"panel": "A1", "extended": "A2"},
    "assay_inference": {"until_year": 104, "assay": "A1", "source_label": "inferred"},
}


def params(**changes) -> dict:
    d = {"study_design": copy.deepcopy(DESIGN)}
    for path, value in changes.items():
        node = d["study_design"]
        *parents, leaf = path.split("__")
        for k in parents:
            node = node[k]
        if value is ...:
            del node[leaf]
        else:
            node[leaf] = value
    return d


def test_a_valid_design_is_parsed():
    d = parse_study_design(params())
    assert d.primary.label == "P-A" and d.sensitivity.last_year == 119
    assert d.panel_regions == ("X1", "X2")
    assert d.assays == ("A1", "A2")
    assert d.period_of(105) == "P-A" and d.period_of(119) == "P-B"
    assert d.region_rank("X2") == 1 and d.region_rank("other") == 2


def test_the_section_is_required():
    with pytest.raises(StudyParameterError, match="study_design"):
        parse_study_design({})


@pytest.mark.parametrize("change, message", [
    ({"periods__primary": ...}, "periods.primary"),
    ({"periods__primary__label": ""}, "label"),
    ({"periods__primary__first_year": "100"}, "integers"),
    ({"periods__primary__first_year": 120}, "after"),
    ({"periods__sensitivity__first_year": 105}, "overlap"),
    ({"periods__sensitivity__label": "P-A"}, "different labels"),
    ({"panel_regions": []}, "panel_regions"),
    ({"panel_regions": ["X1", "X1"]}, "panel_regions"),
    ({"assays__extended": "A1"}, "two different"),
    ({"assays__panel": ...}, "assays.panel"),
    ({"assay_inference__assay": "A9"}, "one of the two assay labels"),
    ({"assay_inference__until_year": ...}, "until_year"),
    ({"assay_inference": ...}, "assay_inference"),
])
def test_malformed_designs_are_rejected(change, message):
    with pytest.raises(StudyParameterError, match=message):
        parse_study_design(params(**change))


def test_assay_inference_may_be_switched_off_explicitly():
    d = parse_study_design(params(assay_inference=None))
    assert d.inference_until_year is None and d.inferred_assay is None


def test_a_year_outside_both_periods_is_an_error():
    d = parse_study_design(params())
    with pytest.raises(StudyParameterError, match="outside both study periods"):
        d.period_of(99)


def test_loader_assigns_periods_from_the_design_and_rejects_uncovered_years(tmp_path):
    d = parse_study_design(params())
    pat = pd.DataFrame({"patient_id": ["a", "b"], "year_requested": ["101", "115"],
                        "n_variants": ["0", "1"]})
    var = pd.DataFrame({"patient_id": ["b"], "variant": ["V"], "year_requested": ["115"]})
    pat.to_csv(tmp_path / "p.csv", sep=";", index=False)
    var.to_csv(tmp_path / "v.csv", sep=";", index=False)
    p, v = load(tmp_path / "p.csv", tmp_path / "v.csv", d)
    assert list(p.era) == ["P-A", "P-B"] and list(v.era) == ["P-B"]
    pat.loc[0, "year_requested"] = "150"
    pat.to_csv(tmp_path / "p.csv", sep=";", index=False)
    with pytest.raises(StudyParameterError):
        load(tmp_path / "p.csv", tmp_path / "v.csv", d)


def test_missing_parameters_file_has_a_clear_error(tmp_path):
    code = "from mefv import config; config.study_design()"
    env = dict(os.environ, MEFV_STUDY_PARAMETERS=str(tmp_path / "absent.json"),
               PYTHONPATH=str(ROOT / "src"))
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.returncode != 0
    assert "study parameters not found" in r.stderr


def test_invalid_json_has_a_clear_error(tmp_path):
    bad = tmp_path / "p.json"
    bad.write_text("{not json")
    code = "from mefv import config; config.load_study_parameters()"
    env = dict(os.environ, MEFV_STUDY_PARAMETERS=str(bad), PYTHONPATH=str(ROOT / "src"))
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.returncode != 0 and "not valid JSON" in r.stderr


def test_package_imports_without_any_private_file(tmp_path):
    code = ("import mefv, mefv.config, mefv.reference, mefv.annotation, mefv.statistics, "
            "mefv.utils; print('ok')")
    env = dict(os.environ, MEFV_STUDY_PARAMETERS=str(tmp_path / "absent.json"),
               PYTHONPATH=str(ROOT / "src"))
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                       env=env, cwd=tmp_path)
    assert r.returncode == 0 and r.stdout.strip() == "ok", r.stderr


def test_config_module_holds_no_study_design_values():
    """No module-level constant may carry study periods, regions or assays."""
    tree = ast.parse((ROOT / "src/mefv/config.py").read_text())
    names = {t.id for n in tree.body if isinstance(n, ast.Assign)
             for t in n.targets if isinstance(t, ast.Name)}
    for forbidden in ("PANEL_EXONS", "ERA_EARLY", "ERA_LATE", "COMPLETE_ERA_MAX",
                      "SYNONYMOUS", "UNVERIFIED", "BENIGN_COMMON"):
        assert forbidden not in names
    assert not hasattr(config, "PANEL_EXONS")


def test_example_parameters_file_is_complete_and_valid():
    example = json.loads((ROOT / "study_parameters.example.json").read_text())
    parse_study_design(example)
    for key in ("variant_groups", "linkage_plan"):
        assert key in example
    assert set(config.VARIANT_GROUP_KEYS) <= set(example["variant_groups"])
    assert set(config.LINKAGE_PLAN_KEYS) <= set(example["linkage_plan"])
    labels = json.dumps(example)
    assert "EXAMPLE" in labels
