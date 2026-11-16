"""Structural checks that the synthetic demonstration is independent.

A banner saying "synthetic" proves nothing about where the numbers came from.
These tests check the generator itself: what it imports, what it reads, and
whether its design is defined by itself rather than taken from a study
configuration. They cannot compare against the private study design, which is
not in the repository; they check that the generator has no path by which to
obtain it.
"""
from __future__ import annotations

import ast
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from mefv import config

ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "scripts" / "make_synthetic_data.py"

def load_generator():
    spec = importlib.util.spec_from_file_location("make_synthetic_data", GENERATOR)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ALLOWED_IMPORTS = {"__future__", "argparse", "json", "pathlib", "sys",
                   "numpy", "pandas", "mefv", "mefv.annotation"}


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text())
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            names.add(mod)
            names |= {f"{mod}.{a.name}" for a in node.names if mod == "mefv"}
    return names


def test_generator_imports_no_study_configuration():
    imported = _imports(GENERATOR)
    assert imported <= ALLOWED_IMPORTS, sorted(imported - ALLOWED_IMPORTS)
    assert "mefv.config" not in imported


def test_generator_source_contains_no_file_reading_call():
    tree = ast.parse(GENERATOR.read_text())
    readers = {"read_csv", "read_text", "read_bytes", "read_excel", "load", "loads"}
    called = {n.func.attr for n in ast.walk(tree)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)}
    called |= {n.func.id for n in ast.walk(tree)
               if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert not called & readers, sorted(called & readers)
    assert "open" not in called


def test_generator_opens_no_data_file_at_run_time(tmp_path):
    """Run the generator with an audit hook and list every file it reads."""
    out = tmp_path / "out"
    probe = f"""
import sys, os
reads = []
def hook(event, args):
    if event == "open" and args and not isinstance(args[0], int):
        p = os.path.realpath(os.fsdecode(args[0]))
        if not p.startswith({str(out.resolve())!r}):
            reads.append(p)
import runpy
mod = runpy.run_path({str(GENERATOR)!r})
sys.addaudithook(hook)
data = mod["generate"](n=60, seed=3)
mod["write"](__import__("pathlib").Path({str(out)!r}), data)
data_like = [p for p in reads if p.endswith((".csv", ".tsv", ".json", ".xlsx", ".txt"))]
print("READS:" + repr(data_like))
"""
    env = dict(os.environ, MEFV_STUDY_PARAMETERS=str(tmp_path / "absent.json"))
    r = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True,
                       env=env, cwd=tmp_path)
    assert r.returncode == 0, r.stderr[-2000:]
    line = next(x for x in r.stdout.splitlines() if x.startswith("READS:"))
    assert line == "READS:[]"


def test_synthetic_variant_labels_are_invented():
    syn = load_generator()
    for name, region, _, _ in syn.CATALOGUE:
        assert re.fullmatch(r"SYN_[A-Z]\d{2}", name), name
        # nothing in the shape of a real protein-level variant name (e.g. A123B)
        assert not re.fullmatch(r"[A-Z]\d+[A-Z]", name.removeprefix("SYN_"))
        assert re.fullmatch(r"R\d+", region), region


def test_synthetic_design_is_self_defined_and_valid():
    syn = load_generator()
    d = syn.generate(n=150, seed=5)
    params = d["parameters"]
    assert params["study_design"] is syn.SYNTHETIC_DESIGN
    design = config.parse_study_design(params)
    # every generated year falls in the invented periods
    years = d["patients"].year_requested.astype(int)
    assert all(design.primary.contains(y) or design.sensitivity.contains(y) for y in years)
    # assays and regions are the invented labels
    assert set(d["patients"].assay) == set(design.assays)
    assert set(d["variants_long"].exon) <= {*design.panel_regions, *syn.EXTENDED_ONLY_REGIONS}
    # assay use is not tied to the period: both assays occur in both periods
    pat = d["patients"].assign(period=[design.period_of(y) for y in years])
    assert pat.groupby("period").assay.nunique().min() == 2


def test_provenance_records_independence(tmp_path):
    syn = load_generator()
    syn.write(tmp_path, syn.generate(n=40, seed=2))
    prov = json.loads((tmp_path / "PROVENANCE.json").read_text())
    assert prov["synthetic"] is True
    assert prov["inputs_read"] == []
    assert prov["design"] == syn.SYNTHETIC_DESIGN
    assert "invented" in prov["statement"]
