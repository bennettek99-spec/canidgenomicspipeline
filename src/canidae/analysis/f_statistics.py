"""Frequency-based f4, D, f4-ratio and qpAdm statistics with a weighted block jackknife.

Genotypes are alt-allele counts (0/1/2) with ``-1`` for a missing call, shaped
``(n_sites, n_samples)``. Population allele frequencies ignore missing calls; a
site enters a statistic only when every population in it has at least one call.
For low-coverage data, :func:`gl_population_frequencies` estimates the same
frequencies by maximum likelihood from phred-scaled genotype likelihoods instead.

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


@dataclass(frozen=True)
class QpAdmResult:
    """Mixture weights for a target, one per source, with a model-fit test."""

    sources: tuple[str, ...]
    weights: np.ndarray
    se: np.ndarray
    chisq: float
    dof: int
    p_value: float
    n_sites: int
    n_blocks: int


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


def gl_population_frequencies(
    phred_likelihoods: np.ndarray,
    samples: Sequence[str],
    groups: Mapping[str, Sequence[str]],
    *,
    n_iter: int = 60,
    tol: float = 1e-7,
) -> dict[str, np.ndarray]:
    """Maximum-likelihood alt-allele frequency per site from genotype likelihoods.

    ``phred_likelihoods`` is ``(n_sites, n_samples, 3)`` phred-scaled PL values for
    genotypes 0/0, 0/1, 1/1. An all-zero PL is uninformative and treated as missing.
    Frequencies are fitted by EM under Hardy-Weinberg proportions (Kim et al. 2011),
    so low-depth heterozygotes are neither forced to a hard call nor discarded.
    """
    index = {sample: i for i, sample in enumerate(samples)}
    freqs: dict[str, np.ndarray] = {}
    for name, members in groups.items():
        pl: np.ndarray = phred_likelihoods[:, [index[m] for m in members], :].astype(float)
        informative = pl.sum(axis=2) > 0
        likelihood = np.power(10.0, -pl / 10.0) * informative[:, :, None]
        n_informative = informative.sum(axis=1)
        p = np.full(pl.shape[0], 0.5)
        for _ in range(n_iter):
            prior = np.stack([(1 - p) ** 2, 2 * p * (1 - p), p**2], axis=1)[:, None, :]
            joint = likelihood * prior
            total = joint.sum(axis=2)
            with np.errstate(invalid="ignore", divide="ignore"):
                dosage = np.where(total > 0, (joint[:, :, 1] + 2 * joint[:, :, 2]) / total, 0.0)
                updated = np.where(
                    n_informative > 0, dosage.sum(axis=1) / (2.0 * n_informative), 0.5
                )
            converged = np.max(np.abs(updated - p)) < tol
            p = updated
            if converged:
                break
        freqs[name] = np.where(n_informative > 0, p, np.nan)
    return freqs


def _jackknife_pseudovalues(
    theta: np.ndarray, theta_minus: np.ndarray, sizes: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pseudovalues, their jackknife mean and h for (possibly vector) estimates.

    ``theta_minus`` has one leading row per block; ``sizes`` are block site counts.
    """
    g = len(sizes)
    n = float(sizes.sum())
    h = (n / sizes).reshape((g,) + (1,) * (theta_minus.ndim - 1))
    theta_j = g * theta - np.sum((1.0 - 1.0 / h) * theta_minus, axis=0)
    pseudo = h * theta - (h - 1.0) * theta_minus
    return pseudo, theta_j, h


def _jackknife_se(theta: np.ndarray, theta_minus: np.ndarray, sizes: np.ndarray) -> np.ndarray:
    pseudo, theta_j, h = _jackknife_pseudovalues(theta, theta_minus, sizes)
    # A block holding every site (h == 1) carries no leave-one-out information.
    with np.errstate(invalid="ignore", divide="ignore"):
        terms = np.where(h > 1.0, (pseudo - theta_j) ** 2 / (h - 1.0), 0.0)
    return np.sqrt(terms.sum(axis=0) / len(sizes))


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
    se = _jackknife_se(np.array(theta), theta_minus, sizes)
    return JackknifeEstimate(theta, float(se), n, g)


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


def _qpadm_fit(y: np.ndarray, m: np.ndarray) -> np.ndarray:
    """Weights for sources 1..n-1 (source 0 takes the remainder) by least squares."""
    return np.linalg.lstsq(m, y, rcond=None)[0]


def qpadm(
    freqs: Mapping[str, np.ndarray],
    target: str,
    sources: Sequence[str],
    rights: Sequence[str],
    blocks: np.ndarray,
) -> QpAdmResult:
    """Fit ``target`` as a mixture of ``sources`` using ``rights`` as outgroups.

    With statistics y_i = f4(target, S0; R0, Ri) and M_is = f4(S_s, S0; R0, Ri)
    for i = 1..m and s = 1..n-1, a correct model satisfies y = M w, where w are the
    weights of sources 1..n-1 and source 0 takes 1 - sum(w) (Haak et al. 2015). The
    weights come from least squares; their SEs and the residual covariance from the
    block jackknife. ``chisq`` tests the m - (n - 1) remaining constraints: a small
    p-value rejects the model. ``rights[0]`` is the base outgroup R0.
    """
    if len(sources) < 2 or len(rights) < len(sources):
        raise ValueError("qpadm needs >= 2 sources and at least as many rights as sources")
    names = [target, *sources, *rights]
    keep = np.all([np.isfinite(freqs[name]) for name in names], axis=0)
    f = {name: freqs[name][keep] for name in names}
    blk = blocks[keep]
    s0, r0 = sources[0], rights[0]
    y_terms = np.stack([(f[target] - f[s0]) * (f[r0] - f[r]) for r in rights[1:]], axis=1)
    m_terms = np.stack(
        [
            np.stack([(f[s] - f[s0]) * (f[r0] - f[r]) for s in sources[1:]], axis=1)
            for r in rights[1:]
        ],
        axis=1,
    )  # (sites, rights - 1, sources - 1)
    labels, inverse = np.unique(blk, return_inverse=True)
    g, n = len(labels), int(keep.sum())
    if g < 2:
        raise ValueError("qpadm needs at least two jackknife blocks")
    sizes: np.ndarray = np.bincount(inverse, minlength=g).astype(float)
    y_blocks = np.stack([np.bincount(inverse, weights=col, minlength=g) for col in y_terms.T], 1)
    m_blocks = np.zeros((g, *m_terms.shape[1:]))
    for i in range(m_terms.shape[1]):
        for j in range(m_terms.shape[2]):
            m_blocks[:, i, j] = np.bincount(inverse, weights=m_terms[:, i, j], minlength=g)
    y = y_blocks.sum(0) / n
    m = m_blocks.sum(0) / n
    if np.linalg.matrix_rank(m) < len(sources) - 1:
        raise ValueError("rights do not distinguish the sources (rank-deficient f4 matrix)")

    def full_weights(partial: np.ndarray) -> np.ndarray:
        return np.concatenate([[1.0 - partial.sum()], partial])

    w = _qpadm_fit(y, m)
    residual = y - m @ w
    w_minus = np.empty((g, len(sources)))
    r_minus = np.empty((g, len(y)))
    for b in range(g):
        y_b = (y_blocks.sum(0) - y_blocks[b]) / (n - sizes[b])
        m_b = (m_blocks.sum(0) - m_blocks[b]) / (n - sizes[b])
        w_b = _qpadm_fit(y_b, m_b)
        w_minus[b] = full_weights(w_b)
        r_minus[b] = y_b - m_b @ w_b
    weights = full_weights(w)
    se = _jackknife_se(weights, w_minus, sizes)

    dof = len(rights) - len(sources)
    if dof > 0:
        pseudo, mean, h = _jackknife_pseudovalues(residual, r_minus, sizes)
        centered = (pseudo - mean) / np.sqrt(h - 1.0)
        cov = centered.T @ centered / g
        # The fitted residual lives in a dof-dimensional subspace; invert only that.
        eigval, eigvec = np.linalg.eigh(cov)
        top = np.argsort(eigval)[::-1][:dof]
        projected = eigvec[:, top].T @ residual
        chisq = float(np.sum(projected**2 / eigval[top]))
        from scipy.stats import chi2

        p_value = float(chi2.sf(chisq, dof))
    else:
        chisq, p_value = 0.0, float("nan")
    return QpAdmResult(tuple(sources), weights, se, chisq, dof, p_value, n, g)
