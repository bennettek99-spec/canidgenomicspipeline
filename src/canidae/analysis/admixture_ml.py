"""Three-way coyote/wolf/dog maximum-likelihood admixture on bridge loci."""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

from canidae.analysis.genotypes import allele_count
from canidae.analysis.reference_mixture import MIN_EFFECTIVE_P

GRID_STEP = 0.01
BOOTSTRAP_N = 200
# Empirical-Bayes shrinkage: every panel's allele frequency is pulled toward
# the pooled all-canid frequency at that locus with the strength of
# SHRINKAGE_PSEUDOCOUNT virtual observations.
SHRINKAGE_PSEUDOCOUNT = 6.0
# Loci must separate coyote from at least one admixture source (on raw
# frequencies) to inform the fit; uninformative loci only add platform noise.
MIN_PANEL_SEPARATION = 0.15
MIN_DOG_FRACTION_FOR_BREED = 0.05


def wgs_allele_count_matrix(
    keys: list[tuple[str, int]],
    records: dict[tuple[str, int], list[str | None]],
) -> np.ndarray:
    """ALT allele-count matrix (all WGS samples x loci) built once and reused."""
    matrix = np.full((len(records[keys[0]]), len(keys)), np.nan)
    for row in range(matrix.shape[0]):
        for position, key in enumerate(keys):
            genotype = records[key][row]
            if genotype is not None:
                matrix[row, position] = allele_count(genotype)
    return matrix


def panel_frequencies(
    matrix: np.ndarray, indices: list[int], pooled_prior: np.ndarray
) -> np.ndarray:
    """Shrunk allele frequencies; panels of different sizes stay comparable."""
    panel = matrix[indices]
    counts = np.nansum(panel, axis=0)
    called = np.sum(~np.isnan(panel), axis=0)
    freqs = (counts + SHRINKAGE_PSEUDOCOUNT * pooled_prior) / (
        2.0 * called + SHRINKAGE_PSEUDOCOUNT
    )
    freqs[called == 0] = np.nan
    return freqs


def loglik_grid(
    alt_counts: np.ndarray, p_coyote: np.ndarray, p_wolf: np.ndarray, p_dog: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Log-likelihood over a (f_wolf, f_dog) grid for one sample.

    Vectorized: p_eff = (1 - fw - fd) * p_c + fw * p_w + fd * p_d per locus,
    with the pairwise-independence binomial likelihood on observed ALT allele
    counts in {0, 1, 2} (NaN marks an uncalled locus).
    """
    steps = round(1.0 / GRID_STEP) + 1
    grid = np.linspace(0.0, 1.0, steps)
    fw, fd = np.meshgrid(grid, grid, indexing="ij")
    valid = fw + fd <= 1.0 + 1e-9
    fw_v, fd_v = fw[valid], fd[valid]
    called = ~np.isnan(alt_counts)
    observed = alt_counts[called]
    p_obs_c = p_coyote[called]
    p_obs_w = p_wolf[called]
    p_obs_d = p_dog[called]
    p_eff = (
        (1.0 - fw_v[:, None] - fd_v[:, None]) * p_obs_c[None, :]
        + fw_v[:, None] * p_obs_w[None, :]
        + fd_v[:, None] * p_obs_d[None, :]
    )
    p_eff = np.clip(p_eff, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
    alt2 = observed[None, :] > 1.0
    het = (observed[None, :] > 0.0) & ~alt2
    ref2 = ~(alt2 | het)
    ll = np.zeros_like(p_eff)
    ll += alt2 * (2.0 * np.log(p_eff))
    ll += het * (np.log(2.0) + np.log(p_eff) + np.log1p(-p_eff))
    ll += ref2 * (2.0 * np.log1p(-p_eff))
    total = ll.sum(axis=1)
    return fw_v, fd_v, total


def best_mixture(
    alt_counts: np.ndarray, p_coyote: np.ndarray, p_wolf: np.ndarray, p_dog: np.ndarray
) -> tuple[float, float, float]:
    fw_v, fd_v, total = loglik_grid(alt_counts, p_coyote, p_wolf, p_dog)
    best = int(np.argmax(total))
    return float(fw_v[best]), float(fd_v[best]), float(total[best])


def cohort_grid(
    alt_counts: np.ndarray,
    chromosome_counts: np.ndarray,
    p_coyote: np.ndarray,
    p_wolf: np.ndarray,
    p_dog: np.ndarray,
) -> tuple[float, float, float]:
    """ML mixture fit for a pooled cohort (binomial on per-locus alt counts)."""
    steps = round(1.0 / GRID_STEP) + 1
    grid = np.linspace(0.0, 1.0, steps)
    fw, fd = np.meshgrid(grid, grid, indexing="ij")
    valid = fw + fd <= 1.0 + 1e-9
    fw_v, fd_v = fw[valid], fd[valid]
    use = chromosome_counts > 0
    n = chromosome_counts[use]
    k = np.rint(alt_counts[use]).astype(np.int64)
    pc, pw, pd = p_coyote[use], p_wolf[use], p_dog[use]
    p_eff = (
        (1.0 - fw_v[:, None] - fd_v[:, None]) * pc[None, :]
        + fw_v[:, None] * pw[None, :]
        + fd_v[:, None] * pd[None, :]
    )
    p_eff = np.clip(p_eff, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
    ll = (
        gammaln(n + 1)
        - gammaln(k + 1)
        - gammaln(n - k + 1)
        + k * np.log(p_eff)
        + (n - k) * np.log1p(-p_eff)
    ).sum(axis=1)
    best = int(np.argmax(ll))
    return float(fw_v[best]), float(fd_v[best]), float(ll[best])


def bootstrap_ci(
    alt_counts: np.ndarray,
    p_coyote: np.ndarray,
    p_wolf: np.ndarray,
    p_dog: np.ndarray,
    rng: np.random.Generator,
    *,
    bootstrap_n: int = BOOTSTRAP_N,
) -> dict[str, list[float]]:
    """Percentile CIs for (f_wolf, f_dog) from a locus bootstrap."""
    keep = ~np.isnan(alt_counts)
    n = int(np.sum(keep))
    wolf_draws: list[float] = []
    dog_draws: list[float] = []
    for _ in range(bootstrap_n):
        pick = rng.integers(0, n, n)
        d = alt_counts[keep][pick]
        pc, pw, pd = p_coyote[keep][pick], p_wolf[keep][pick], p_dog[keep][pick]
        fw, fd, _ = best_mixture(d, pc, pw, pd)
        wolf_draws.append(fw)
        dog_draws.append(fd)
    return {
        "wolf": [round(v, 3) for v in np.percentile(wolf_draws, [2.5, 97.5])],
        "dog": [round(v, 3) for v in np.percentile(dog_draws, [2.5, 97.5])],
    }


def calls_mixture_loglik(
    alt_counts: np.ndarray,
    p_coyote: np.ndarray,
    p_wolf: np.ndarray,
    p_dog_breed: np.ndarray,
    f_wolf: float,
    f_dog: float,
) -> float:
    """Log-likelihood of one sample under a fixed mixture and a given dog panel."""
    keep = ~np.isnan(alt_counts)
    observed = alt_counts[keep]
    p_eff = (
        (1.0 - f_wolf - f_dog) * p_coyote[keep]
        + f_wolf * p_wolf[keep]
        + f_dog * p_dog_breed[keep]
    )
    p_eff = np.clip(p_eff, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
    alt2 = observed > 1.0
    het = (observed > 0.0) & ~alt2
    ref2 = ~(alt2 | het)
    ll = float(
        2.0 * np.sum(np.log(p_eff[alt2]))
        + np.sum(np.log(2.0) + np.log(p_eff[het]) + np.log1p(-p_eff[het]))
        + 2.0 * np.sum(np.log1p(-p_eff[ref2]))
    )
    return ll


def two_way_fit(
    alt_counts: np.ndarray,
    chromosome_counts: np.ndarray,
    p_from: np.ndarray,
    p_to: np.ndarray,
) -> tuple[float, float]:
    """Constrained 2-way fit (f of p_from -> p_to); returns (f, loglik)."""
    use = chromosome_counts > 0
    n = chromosome_counts[use]
    k = np.rint(alt_counts[use]).astype(np.int64)
    a, b = p_from[use], p_to[use]
    best = (0.0, -np.inf)
    for f in np.linspace(0.0, 1.0, 101):
        p_eff = np.clip((1.0 - f) * a + f * b, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
        ll = float(
            (
                gammaln(n + 1)
                - gammaln(k + 1)
                - gammaln(n - k + 1)
                + k * np.log(p_eff)
                + (n - k) * np.log1p(-p_eff)
            ).sum()
        )
        if ll > best[1]:
            best = (float(f), ll)
    return best
