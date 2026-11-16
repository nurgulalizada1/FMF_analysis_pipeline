#!/usr/bin/env python3
"""Generate an artificial dataset so the pipeline can be run without study data.

The output is **not** study data and must never be presented as a study result.

Independence
------------
The generator reads no file and imports no study configuration: its only
inputs are the constants in this module, the command-line arguments and a
random seed. It imports ``mefv.annotation`` solely for the column names of the
annotation schema. ``tests/test_synthetic_independence.py`` checks this
structurally, by inspecting the imports and by recording every file the
generator opens.

Everything below was invented for the demonstration and designed without
reference to any real study:

* **Study design** (``SYNTHETIC_DESIGN``): two invented periods, two invented
  assay labels, and invented region labels (``R1`` ... ``R6``) rather than real
  exon numbers. Assay use does not depend on the period.
* **Variant labels** are invented (``SYN_P01`` ...); coding positions, protein
  descriptions and rsIDs are placeholders.
* **Allele frequencies** are round numbers chosen by hand.
* **Genotypes** are drawn independently per individual, except for one pair
  placed in deliberate linkage, one pair given identical genotypes (a candidate
  complex allele), and two variants located outside the panel regions that are
  observable only with the extended assay.
* **Identifiers** have the form ``SYNTH_0001``.
* **Export defects**: one record missing from the patient table, one blank assay
  field and one wrong assay field, each resolved by a confirmed override, one
  mistyped label, and blank assay fields in early records resolved by the
  documented inference rule. These exist only
  to exercise each correction path of stage 01.

The generator is seeded, so the same seed reproduces the same files.

Usage
-----
    python scripts/make_synthetic_data.py --out data/synthetic --n 400 --seed 20260101
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

# Make the package importable when the script is run without installing it.
import pathlib as _pathlib
import sys as _sys

_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from mefv import annotation  # noqa: E402  (schema column names only)

BANNER = (
    "# SYNTHETIC DATA - artificially generated for software demonstration only. "
    "Contains no real patient records and reproduces no study result."
)

#: The invented study design, in the format of the study-parameters file.
SYNTHETIC_DESIGN = {
    "periods": {
        "primary": {"label": "SYN-PERIOD-A", "first_year": 2001, "last_year": 2006},
        "sensitivity": {"label": "SYN-PERIOD-B", "first_year": 2007, "last_year": 2010},
    },
    "panel_regions": ["R1", "R2", "R3", "R4"],
    "assays": {"panel": "ASSAY_T", "extended": "ASSAY_W"},
    "assay_inference": {
        "until_year": 2003,
        "assay": "ASSAY_T",
        "source_label": "inferred_from_synthetic_statement",
    },
}

#: Regions covered only by the extended assay.
EXTENDED_ONLY_REGIONS = ("R5", "R6")

#: Probability that an individual was tested with the extended assay, the same
#: in every year.
P_EXTENDED_ASSAY = 0.35

#: Probability that an early, panel-assay record has a blank assay field.
P_BLANK_EARLY = 0.6

#: Label, region, allele frequency, and synthetic classification. All invented.
CATALOGUE: tuple[tuple[str, str, float, str], ...] = (
    ("SYN_P01", "R3", 0.07, "Pathogenic"),
    ("SYN_P02", "R1", 0.05, "Pathogenic"),
    ("SYN_P03", "R4", 0.04, "Pathogenic"),
    ("SYN_P04", "R2", 0.03, "Likely pathogenic"),
    ("SYN_V01", "R2", 0.09, "Uncertain significance"),
    ("SYN_V02", "R2", 0.06, "Uncertain significance"),
    ("SYN_V03", "R4", 0.03, "Uncertain significance"),
    ("SYN_B01", "R1", 0.18, "Likely benign"),
    ("SYN_B02", "R3", 0.07, "Likely benign"),
    ("SYN_B03", "R4", 0.03, "Benign"),
    ("SYN_S01", "R3", 0.12, "Benign"),
    ("SYN_S02", "R3", 0.12, "Benign"),
    ("SYN_S03", "R1", 0.09, "Benign"),
    ("SYN_S04", "R4", 0.11, "Benign"),
    ("SYN_S05", "R2", 0.08, "Benign"),
    ("SYN_S06", "R4", 0.06, "Benign"),
    # carried by exactly the same individuals: a candidate complex allele
    ("SYN_C01", "R4", 0.02, "Uncertain significance"),
    ("SYN_C02", "R4", 0.02, "Uncertain significance"),
    # outside the panel regions: observable only with the extended assay
    ("SYN_X01", "R5", 0.06, "Uncertain significance"),
    ("SYN_X02", "R6", 0.04, "Likely benign"),
)

SYNONYMOUS = [name for name, _, _, _ in CATALOGUE if name.startswith("SYN_S")]
EXTENDED_ONLY = [name for name, region, _, _ in CATALOGUE if region in EXTENDED_ONLY_REGIONS]

#: One pair is given deliberate linkage disequilibrium so that the linkage stage
#: has a signal to detect. The strength carries no biological meaning.
LD_PAIR = ("SYN_V01", "SYN_V02")
LD_STRENGTH = 0.85

#: Two invented variants given identical genotypes.
COMPLEX_PAIR = ("SYN_C01", "SYN_C02")

#: Invented functional labels, using the published assay vocabulary.
SPECKSEQ_DEMO = {
    "SYN_P01": "Strongly hypermorphic (FMF-like)",
    "SYN_P02": "Moderately hypermorphic (FMF-like)",
    "SYN_P03": "Strongly hypermorphic (FMF-like)",
    "SYN_P04": "Strongly hypermorphic (PAAND-like)",
    "SYN_V01": "Uncertain",
    "SYN_V02": "Uncertain",
    "SYN_V03": "Not hypermorphic",
    "SYN_B01": "Uncertain",
    "SYN_B02": "Not hypermorphic",
    "SYN_B03": "Not hypermorphic",
}

ID_PREFIX = "SYNTH_"
ID_WIDTH = 4

#: A deliberately mistyped label, so stage 01's nomenclature correction runs.
TYPO = ("SYN_P0l", "SYN_P01")   # lower-case L standing in for the digit one

POPULATION_COLUMNS = ("mid", "nfe", "fin", "asj", "sas", "eas", "afr", "amr", "rmi", "ami")

#: Provenance header in the form stage 06 requires of a reference file.
POPULATION_PROVENANCE = (
    "# source      : none - values drawn at random by the synthetic generator\n"
    "# version     : synthetic\n"
    "# retrieved   : not applicable\n"
)


def generate(n: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    design = SYNTHETIC_DESIGN
    periods = design["periods"]
    panel, extended = design["assays"]["panel"], design["assays"]["extended"]
    regions = design["panel_regions"]
    inference = design["assay_inference"]

    ids = [f"{ID_PREFIX}{i:0{ID_WIDTH}d}" for i in range(1, n + 1)]
    years = rng.integers(periods["primary"]["first_year"],
                         periods["sensitivity"]["last_year"] + 1, size=n)
    sex = rng.choice(["female", "male"], size=n)
    assay = np.where(rng.random(n) < P_EXTENDED_ASSAY, extended, panel)

    dosage: dict[str, np.ndarray] = {}
    for name, _, f, _ in CATALOGUE:
        dosage[name] = rng.binomial(2, f, size=n)

    a, b = LD_PAIR
    drop = (dosage[a] == 0) & (rng.random(n) < LD_STRENGTH)
    dosage[b] = np.where(drop, 0, dosage[b])

    c1, c2 = COMPLEX_PAIR
    dosage[c2] = dosage[c1].copy()

    # Regions outside the panel are only seen by the extended assay.
    for name in EXTENDED_ONLY:
        dosage[name] = np.where(assay == extended, dosage[name], 0)

    position = {pid: i for i, pid in enumerate(ids)}
    rows = []
    for name, region, _, tier in CATALOGUE:
        for pid, d in zip(ids, dosage[name], strict=True):
            if d:
                i = position[pid]
                rows.append({
                    "patient_id": pid, "variant": name, "exon": region,
                    "zygosity": "Hom" if d == 2 else "Het",
                    "legacy_class": "Pathogenic" if "athogenic" in tier else "Benign",
                    "year_requested": int(years[i]), "sex": sex[i],
                    "assay": assay[i], "assay_source": "synthetic",
                })
    variants = pd.DataFrame(rows).sort_values(["patient_id", "variant"]).reset_index(drop=True)

    tested = {panel: ",".join(regions),
              extended: ",".join([*regions, *EXTENDED_ONLY_REGIONS])}
    counts = variants.groupby("patient_id").size()
    patients = pd.DataFrame({
        "patient_id": ids,
        "year_requested": years,
        "sex": sex,
        "n_variants": [int(counts.get(p, 0)) for p in ids],
        "variant_detected": ["yes" if counts.get(p, 0) else "no" for p in ids],
        "exons_tested": [tested[x] for x in assay],
        "coverage_source": "synthetic",
        "assay": assay,
        "assay_source": "synthetic",
    })

    # ---- raw laboratory-style export, for stage 01 ----------------------------
    early = years <= inference["until_year"]
    blank = early & (assay == inference["assay"]) & (rng.random(n) < P_BLANK_EARLY)
    no_variant = [p for p in ids if counts.get(p, 0) == 0]
    restored_id = no_variant[0]
    override_id = next(p for p in ids
                       if years[position[p]] > inference["until_year"] and p != restored_id)
    # A record exported with the wrong assay, corrected by a confirmed override.
    corrected_id = next(p for p in ids
                        if years[position[p]] > inference["until_year"]
                        and p not in (restored_id, override_id)
                        and assay[position[p]] == extended)

    raw_variants = variants.drop(columns=["assay", "assay_source"]).copy()
    typo_rows = raw_variants.index[raw_variants.variant == TYPO[1]]
    raw_variants.loc[typo_rows[:1], "variant"] = TYPO[0]

    raw_patients = patients.drop(columns=["assay", "assay_source"])
    raw_patients = raw_patients[raw_patients.patient_id != restored_id].copy()

    method = ["" if blank[i] else assay[i] for i in range(n)]
    method[position[override_id]] = ""
    method[position[corrected_id]] = panel     # wrong in the export; the override fixes it
    raw_method = pd.DataFrame({
        # identifiers are written in mixed case, as exports often are
        "Patient_ID": [(p.lower() if i % 7 == 0 else p) for i, p in enumerate(ids)],
        "YEAR": years,
        "sex": sex,
        "Method": method,
    })

    source = pd.DataFrame([{
        "variant": name,
        "infevers_name": name,
        "hgvs_c": f"c.{100 + i}A>G",
        "hgvs_p": f"p.(Xaa{100 + i}Yaa)",
        "exon_infevers": region,
        "rsid": f"rs{900000 + i}",
        "infevers_class": tier,
        "previous_class": "",
        "comment": "synthetic entry; not a database record",
    } for i, (name, region, _, tier) in enumerate(CATALOGUE, start=1)])

    overlay = pd.DataFrame([{
        "variant": name,
        "speckseq": SPECKSEQ_DEMO.get(name, ""),
        "other_functional": "",
        "speckseq_assay_level_infevers": SPECKSEQ_DEMO.get(name, ""),
        "speckseq_integrated_S7": "not tested",
        "bronnec2026_revised_S9": "not listed",
    } for name, _, _, _ in CATALOGUE])[[annotation.KEY, *annotation.OVERLAY_COLUMNS]]

    non_syn = [name for name, _, _, _ in CATALOGUE if name not in SYNONYMOUS]
    pop_rows = []
    for name in non_syn:
        vals = {c: round(float(rng.uniform(1e-4, 0.05)), 6) for c in POPULATION_COLUMNS}
        pop_rows.append({
            "variant": name, **vals,
            "glob": round(float(np.mean(list(vals.values()))), 6),
            "grpmax": max(vals.values()),
            "hom": int(rng.integers(0, 50)),
            "xx": round(float(rng.uniform(1e-4, 0.05)), 6),
            "xy": round(float(rng.uniform(1e-4, 0.05)), 6),
            "revel": round(float(rng.uniform(0, 1)), 3),
            "cadd": round(float(rng.uniform(0, 30)), 1),
            "cv": "synthetic",
            "stars": int(rng.integers(0, 3)),
        })
    population = pd.DataFrame(pop_rows)

    pathogenic = [name for name, _, _, tier in CATALOGUE
                  if tier in ("Pathogenic", "Likely pathogenic")]
    parameters = {
        "_comment": "SYNTHETIC parameters written by scripts/make_synthetic_data.py.",
        "study_design": SYNTHETIC_DESIGN,
        "expected_cohort_size": n,
        "id_prefix": ID_PREFIX,
        "id_width": ID_WIDTH,
        "cohort_label": f"SYNTHETIC (n={n})",
        "assay_overrides": {override_id: str(assay[position[override_id]]),
                            corrected_id: str(assay[position[corrected_id]])},
        "nomenclature_fixes": {TYPO[0]: TYPO[1]},
        "nomenclature_fix_notes": {TYPO[0]: "deliberate synthetic typo, exercising the correction path"},
        "variant_groups": {
            "synonymous": SYNONYMOUS,
            "block_synonymous": ["SYN_S01", "SYN_S02"],
            "benign_common": ["SYN_B01"],
            "unverified": [],
        },
        "linkage_plan": {
            "primary_pair": list(LD_PAIR),
            "haplotype_blocks": {"block A": ["SYN_S01", "SYN_S02"],
                                 "block B": ["SYN_S04", "SYN_S06"]},
            "cross_locus_pairs": {"SYN_S01": list(LD_PAIR), "SYN_S04": list(LD_PAIR)},
        },
        "affinity_comparison_variants": [v for v in non_syn if v not in pathogenic],
    }

    provenance = {
        "synthetic": True,
        "generator": "scripts/make_synthetic_data.py",
        "seed": seed,
        "n_individuals": n,
        "inputs_read": [],
        "design": SYNTHETIC_DESIGN,
        "allele_frequencies": {name: f for name, _, f, _ in CATALOGUE},
        "linkage_pair": list(LD_PAIR),
        "linkage_strength": LD_STRENGTH,
        "complex_allele_pair": list(COMPLEX_PAIR),
        "extended_assay_only_variants": EXTENDED_ONLY,
        "statement": (
            "Artificially generated data for software demonstration. Not sampled from, "
            "fitted to, or conditioned on any real cohort; the study design, labels and "
            "frequencies are invented. The generator reads no file. Results produced "
            "from these files are not study findings."
        ),
    }
    return {
        "raw_method": raw_method,
        "raw_variants": raw_variants,
        "raw_patients": raw_patients,
        "restored_id": restored_id,
        "override_id": override_id,
        "corrected_id": corrected_id,
        "population": population,
        "source": source,
        "patients": patients,
        "variants_long": variants,
        "overlay": overlay,
        "parameters": parameters,
        "provenance": provenance,
    }


def write(out: Path, data: dict) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for sub in ("raw", "reviewed", "external"):
        (out / sub).mkdir(exist_ok=True)
    for name, target in (("raw_method", out / "raw/synthetic_method.csv"),
                         ("raw_variants", out / "raw/synthetic_variants.csv"),
                         ("raw_patients", out / "raw/synthetic_patients.csv"),
                         ("source", out / "external/source_annotations.csv"),
                         ("population", out / "external/population_frequencies.csv"),
                         ("overlay", out / "reviewed/variant_annotations_reviewed.csv")):
        body = data[name].to_csv(sep=";", index=False)
        header = BANNER + "\n"
        if name == "population":
            header += POPULATION_PROVENANCE
        target.write_text(header + body, encoding="utf-8-sig")

    (out / "reviewed/annotation_decision_log.csv").write_text(
        BANNER + "\n" + ";".join(annotation.DECISION_LOG_COLUMNS) + "\n", encoding="utf-8-sig"
    )
    (out / "study_parameters.json").write_text(
        json.dumps(data["parameters"], indent=2) + "\n")
    (out / "PROVENANCE.json").write_text(
        json.dumps(data["provenance"], indent=2, sort_keys=True) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("data/synthetic"))
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20260101)
    args = ap.parse_args()

    data = generate(args.n, args.seed)
    write(args.out, data)
    print(f"synthetic cohort: {args.n} individuals, seed {args.seed}")
    print(f"  variant observations : {len(data['variants_long'])}")
    print(f"  with >=1 variant     : {int((data['patients'].n_variants > 0).sum())}")
    print(f"  written to           : {args.out}")
    print("\nThese files are artificial. Results produced from them are not study findings.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
