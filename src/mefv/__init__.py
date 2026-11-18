"""Analysis library for MEFV variant data in familial Mediterranean fever.

The package holds the reusable parts of the analysis: study constants, data
loading and genotype encoding, and the statistical functions. The ordered
analysis stages live in ``scripts/`` and import from here.
"""
from mefv import annotation
from mefv.statistics import (
    bh_fdr,
    bootstrap_dprime,
    em_haplotypes,
    fmt_ci,
    fmt_p,
    hwe_exact,
    ld_stats,
    wilson_ci,
)
from mefv.utils import DOSE, allele_dose, genotype_vector, load, read_csv, require_released
from mefv.validation import DataIntegrityError

__version__ = "0.1.0"

__all__ = [
    "DOSE",
    "DataIntegrityError",
    "allele_dose",
    "annotation",
    "bh_fdr",
    "bootstrap_dprime",
    "em_haplotypes",
    "fmt_ci",
    "fmt_p",
    "genotype_vector",
    "hwe_exact",
    "ld_stats",
    "load",
    "read_csv",
    "require_released",
    "wilson_ci",
]
