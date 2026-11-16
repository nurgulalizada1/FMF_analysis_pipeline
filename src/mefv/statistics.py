"""Statistical functions used by the analysis scripts.

Each function is a direct implementation of a published method and is covered
by tests in ``tests/test_statistics.py`` against values that can be derived
independently.

References
----------
Wilson EB. J Am Stat Assoc. 1927;22:209-212.
    Score interval for a binomial proportion.
Wigginton JE, Cutler DJ, Abecasis GR. Am J Hum Genet. 2005;76:887-893.
    Exact test of Hardy-Weinberg equilibrium.
Excoffier L, Slatkin M. Mol Biol Evol. 1995;12:921-927.
    Maximum-likelihood estimation of haplotype frequencies by the EM algorithm.
Lewontin RC. Genetics. 1964;49:49-67.
    Normalised linkage disequilibrium coefficient D'.
Benjamini Y, Hochberg Y. J R Stat Soc B. 1995;57:289-300.
    Control of the false discovery rate.
"""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import stats


def wilson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (np.nan, np.nan)
    z = stats.norm.ppf(1 - alpha / 2)
    p, d = k / n, 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def hwe_exact(n_hom_alt: int, n_het: int, n_hom_ref: int) -> float:
    """Exact test of Hardy-Weinberg equilibrium (Wigginton et al. 2005)."""
    a, h, b = int(n_hom_alt), int(n_het), int(n_hom_ref)
    n = a + h + b
    if n == 0:
        return np.nan
    rare = 2 * min(a, b) + h
    if rare == 0:
        return 1.0
    probs = np.zeros(rare + 1)
    mid = int(rare * (2 * n - rare) / (2 * n))
    if mid % 2 != rare % 2:
        mid += 1
    probs[mid] = 1.0
    total = 1.0
    het, hr = mid, int((rare - mid) / 2)
    ha = n - het - hr
    while het >= 2:
        probs[het - 2] = probs[het] * het * (het - 1.0) / (4.0 * (hr + 1.0) * (ha + 1.0))
        total += probs[het - 2]
        het, hr, ha = het - 2, hr + 1, ha + 1
    het, hr = mid, int((rare - mid) / 2)
    ha = n - het - hr
    while het <= rare - 2:
        probs[het + 2] = probs[het] * 4.0 * hr * ha / ((het + 2.0) * (het + 1.0))
        total += probs[het + 2]
        het, hr, ha = het + 2, hr - 1, ha - 1
    probs /= total
    return float(min(1.0, probs[probs <= probs[h] * (1 + 1e-9)].sum()))


def em_haplotypes(
    gA: Sequence[int], gB: Sequence[int], iters: int = 1000, tol: float = 1e-12
) -> dict[str, float]:
    """EM estimation of two-locus haplotype frequencies from unphased data."""
    n = len(gA)
    pA, pB = sum(gA) / (2 * n), sum(gB) / (2 * n)
    h = {"11": pA * pB, "10": pA * (1 - pB),
         "01": (1 - pA) * pB, "00": (1 - pA) * (1 - pB)}
    for _ in range(iters):
        c = {k: 1e-12 for k in h}
        for a, b in zip(gA, gB):
            if a == 1 and b == 1:                    # the only ambiguous case
                w_cis, w_trans = h["11"] * h["00"], h["10"] * h["01"]
                s = w_cis + w_trans
                if s <= 0:
                    w_cis = w_trans = 0.5
                    s = 1.0
                c["11"] += w_cis / s
                c["00"] += w_cis / s
                c["10"] += w_trans / s
                c["01"] += w_trans / s
            elif a == 1:
                hb = "1" if b == 2 else "0"
                c["1" + hb] += 1
                c["0" + hb] += 1
            elif b == 1:
                ha = "1" if a == 2 else "0"
                c[ha + "1"] += 1
                c[ha + "0"] += 1
            else:
                ha = "1" if a == 2 else "0"
                hb = "1" if b == 2 else "0"
                c[ha + hb] += 2
        tot = sum(c.values())
        new = {k: v / tot for k, v in c.items()}
        if max(abs(new[k] - h[k]) for k in h) < tol:
            return new
        h = new
    return h


#: Allele frequencies this close to 0 or 1 are treated as monomorphic.
#: The EM estimate converges to within about 1e-14 of an exact 0 or 1 rather
#: than reaching it, so an exact comparison never fires. The tolerance sits far
#: above that residue and far below any attainable frequency: the smallest
#: frequency a real locus can show is 1/(2n), which is 5e-5 even for ten
#: thousand individuals - more than four orders of magnitude above the tolerance.
FIXED_TOL = 1e-9


def _is_fixed(p: float) -> bool:
    """Whether an estimated allele frequency is monomorphic within tolerance."""
    return p <= FIXED_TOL or p >= 1.0 - FIXED_TOL


def ld_stats(gA: Sequence[int], gB: Sequence[int]) -> dict[str, float]:
    """Linkage disequilibrium between two biallelic loci from unphased genotypes.

    Returns D, D', r^2, the estimated haplotype frequencies and the two allele
    frequencies. Linkage disequilibrium is undefined when either locus is
    monomorphic, and the statistics are NaN in that case.
    """
    h = em_haplotypes(gA, gB)
    pA, pB = h["11"] + h["10"], h["11"] + h["01"]
    if _is_fixed(pA) or _is_fixed(pB):
        return dict(D=np.nan, Dprime=np.nan, r2=np.nan, h=h, pA=pA, pB=pB)
    D = h["11"] - pA * pB
    dmax = min(pA * (1 - pB), (1 - pA) * pB) if D > 0 else min(pA * pB, (1 - pA) * (1 - pB))
    return dict(D=D, Dprime=D / dmax if dmax > 0 else np.nan,
                r2=D**2 / (pA * (1 - pA) * pB * (1 - pB)), h=h, pA=pA, pB=pB)


def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p values."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


# ------------------------------------------------------------ formatting ----
def fmt_ci(lo: float, hi: float, dp: int = 4) -> str:
    return "" if np.isnan(lo) else f"{lo:.{dp}f}–{hi:.{dp}f}"


def fmt_p(p: float) -> str:
    if np.isnan(p):
        return ""
    return f"{p:.1e}" if p < 1e-4 else f"{p:.4f}"
