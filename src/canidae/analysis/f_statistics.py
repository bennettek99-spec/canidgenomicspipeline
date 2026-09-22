"""Frequency-based f4, D and f4-ratio statistics with a weighted block jackknife.

Genotypes are alt-allele counts (0/1/2) with ``-1`` for a missing call, shaped
``(n_sites, n_samples)``. Population allele frequencies ignore missing calls; a
site enters a statistic only when every population in it has at least one call.

Uncertainty uses the delete-one weighted block jackknife of Busing, Meijer and
van der Leeden (1999), as in ADMIXTOOLS: blocks are contiguous genomic windows
or chromosomes, so linked sites are never treated as independent.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class JackknifeEstimate:
    """A ratio-of-sums statistic with its block-jackknife standard error."""

    estimate: float
    se: float
    n_sites: int
    n_blocks: int

    @property
    def z(self) -> float:
        return self.estimate / self.se if self.se > 0 else float("nan")


def population_frequencies(
    genotypes: np.ndarray, samples: Sequence[str], groups: Mapping[str, Sequence[str]]
) -> dict[str, np.ndarray]:
    """Alt-allele frequency per site for each named group (NaN where uncalled)."""
    index = {sample: i for i, sample in enumerate(samples)}
    freqs: dict[str, np.ndarray] = {}
    for name, members in groups.items():
        cols = [index[m] for m in members]
        block: np.ndarray = genotypes[:, cols].astype(float)
        called = block >= 0
        alt = np.where(called, block, 0.0).sum(axis=1)
        n_alleles = 2.0 * called.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            freqs[name] = np.where(n_alleles > 0, alt / n_alleles, np.nan)
    return freqs


def weighted_block_jackknife(
    numerator: np.ndarray, denominator: np.ndarray, blocks: np.ndarray
) -> JackknifeEstimate:
    """Estimate ``sum(numerator) / sum(denominator)`` with a weighted block jackknife.

    Sites with a non-finite numerator or denominator are dropped. At least two
    non-empty blocks are required for a standard error.
    """
    ok = np.isfinite(numerator) & np.isfinite(denominator)
    num, den, blk = numerator[ok], denominator[ok], blocks[ok]
    n = int(ok.sum())
    total_num, total_den = float(num.sum()), float(den.sum())
    if n == 0 or total_den == 0:
        return JackknifeEstimate(float("nan"), float("nan"), n, 0)
    theta = total_num / total_den
    labels, inverse = np.unique(blk, return_inverse=True)
    g = len(labels)
    if g < 2:
        return JackknifeEstimate(theta, float("nan"), n, g)
    block_num = np.bincount(inverse, weights=num, minlength=g)
    block_den = np.bincount(inverse, weights=den, minlength=g)
    sizes: np.ndarray = np.bincount(inverse, minlength=g).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        theta_minus = (total_num - block_num) / (total_den - block_den)
    h = n / sizes
    theta_j = g * theta - np.sum((1.0 - sizes / n) * theta_minus)
    pseudo = h * theta - (h - 1.0) * theta_minus
    # A block holding every site (h == 1) carries no leave-one-out information.
    informative = h > 1.0
    variance = float(np.sum((pseudo[informative] - theta_j) ** 2 / (h[informative] - 1.0)) / g)
    return JackknifeEstimate(theta, float(np.sqrt(variance)), n, g)


def f4_terms(pa: np.ndarray, pb: np.ndarray, pc: np.ndarray, pd: np.ndarray) -> np.ndarray:
    """Per-site f4(A,B;C,D) = (a-b)(c-d); NaN where any population is uncalled."""
    return (pa - pb) * (pc - pd)


def f4(
    freqs: Mapping[str, np.ndarray], a: str, b: str, c: str, d: str, blocks: np.ndarray
) -> JackknifeEstimate:
    """Mean f4(A,B;C,D) over sites, with block-jackknife SE."""
    terms = f4_terms(freqs[a], freqs[b], freqs[c], freqs[d])
    return weighted_block_jackknife(terms, np.where(np.isfinite(terms), 1.0, np.nan), blocks)


def d_statistic(
    freqs: Mapping[str, np.ndarray], p1: str, p2: str, p3: str, p4: str, blocks: np.ndarray
) -> JackknifeEstimate:
    """Patterson's D(P1,P2;P3,P4) in allele-frequency form.

    Positive values mean P1 shares more derived drift with P3 than P2 does
    (equivalently P2 with P4), given P4 is the outgroup.
    """
    w, x, y, z = freqs[p1], freqs[p2], freqs[p3], freqs[p4]
    num = (w - x) * (y - z)
    den = (w + x - 2 * w * x) * (y + z - 2 * y * z)
    return weighted_block_jackknife(num, den, blocks)


def f4_ratio(
    freqs: Mapping[str, np.ndarray],
    a: str,
    o: str,
    x: str,
    b: str,
    c: str,
    blocks: np.ndarray,
) -> JackknifeEstimate:
    """Patterson's f4-ratio: alpha = f4(A,O;X,C) / f4(A,O;B,C).

    Under the tree (O, ((C, X_C), (A, (B, X_B)))) with X = alpha*X_B + (1-alpha)*X_C,
    alpha is the ancestry X derives from the B side. Only sites where all five
    populations are called contribute to either sum, so both share a site set.
    """
    num = f4_terms(freqs[a], freqs[o], freqs[x], freqs[c])
    den = f4_terms(freqs[a], freqs[o], freqs[b], freqs[c])
    keep = np.isfinite(num) & np.isfinite(den)
    return weighted_block_jackknife(
        np.where(keep, num, np.nan), np.where(keep, den, np.nan), blocks
    )
