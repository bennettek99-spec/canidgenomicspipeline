"""A lightweight HMM for windowed local-ancestry inference.

Given allele frequencies for a set of *source* populations, each genomic window of a target
individual is assigned the source ancestry that best explains its genotypes, with an HMM
transition penalty smoothing assignments along each chromosome (few ancestry switches). The
HMM is restarted at every chromosome boundary; chromosomes are processed sequentially. This is
a genotype-based, single-state-per-window simplification of full diploid local-ancestry
inference (RFMix/Loter), implemented in pure NumPy so it runs anywhere — no GPU, no external
tools. The RFMix/Loter binaries can slot in behind the same stage as alternative backends.
"""

from __future__ import annotations

import numpy as np

_EPS = 1e-6


def genotype_loglik(n_alt: np.ndarray, freq: np.ndarray) -> np.ndarray:
    """Per-site log P(genotype | source allele frequency) under HWE.

    ``n_alt`` in {0,1,2} (ALT-allele count), -1 for missing (contributes 0). The
    binomial coefficient is dropped as it is constant across sources.
    """
    p = np.clip(freq, _EPS, 1 - _EPS)
    ll = np.select(
        [n_alt == 0, n_alt == 1, n_alt == 2],
        [2 * np.log(1 - p), np.log(2) + np.log(p) + np.log(1 - p), 2 * np.log(p)],
        default=0.0,
    )
    return np.where(n_alt < 0, 0.0, ll)


def window_emissions(
    n_alt: np.ndarray,
    chrom: np.ndarray,
    pos: np.ndarray,
    source_freqs: dict[str, np.ndarray],
    *,
    window_bp: int,
    min_sites: int = 5,
) -> tuple[list[str], list[tuple[str, int, int]], np.ndarray]:
    """Aggregate per-site log-likelihoods into per-window emissions.

    Returns (source order, window intervals, emission matrix of shape (n_windows,
    n_sources)).
    """
    sources = list(source_freqs)
    per_site = {s: genotype_loglik(n_alt, source_freqs[s]) for s in sources}
    windows: list[tuple[str, int, int]] = []
    rows: list[list[float]] = []
    for contig in dict.fromkeys(chrom):
        m = chrom == contig
        cpos = pos[m]
        if cpos.size == 0:
            continue
        for w0 in range(int(cpos.min()), int(cpos.max()) + 1, window_bp):
            wmask = (cpos >= w0) & (cpos < w0 + window_bp)
            if int(wmask.sum()) < min_sites:
                continue
            windows.append((str(contig), int(w0), int(w0 + window_bp)))
            rows.append([float(per_site[s][m][wmask].sum()) for s in sources])
    return sources, windows, np.asarray(rows, dtype=float)


def viterbi(emissions: np.ndarray, switch_prob: float) -> np.ndarray:
    """MAP state path over windows given a (n_windows, n_states) emission log-likelihood
    matrix and a per-window ancestry-switch probability."""
    n, k = emissions.shape
    if n == 0:
        return np.empty(0, dtype=int)
    if k == 1:
        return np.zeros(n, dtype=int)
    trans = np.full((k, k), np.log(switch_prob / (k - 1)))
    np.fill_diagonal(trans, np.log(1 - switch_prob))

    v = emissions[0].copy()
    back = np.zeros((n, k), dtype=int)
    for t in range(1, n):
        scores = v[:, None] + trans  # (prev, cur)
        best_prev = np.argmax(scores, axis=0)
        v = emissions[t] + scores[best_prev, np.arange(k)]
        back[t] = best_prev
    path = np.zeros(n, dtype=int)
    path[-1] = int(np.argmax(v))
    for t in range(n - 1, 0, -1):
        path[t - 1] = back[t, path[t]]
    return path


def infer_local_ancestry(
    n_alt: np.ndarray,
    chrom: np.ndarray,
    pos: np.ndarray,
    source_freqs: dict[str, np.ndarray],
    *,
    window_bp: int,
    switch_prob: float = 0.02,
    min_sites: int = 5,
) -> list[dict[str, object]]:
    """Return one MAP ancestry call per genomic window for one individual.

    The Viterbi state is intentionally reset at every chromosome boundary: an ancestry state
    at the end of one chromosome has no biological transition relationship to the next one.
    Processing is also bounded to a single chromosome's emission matrix at a time.
    """
    n_alt_arr = np.asarray(n_alt)
    chrom_arr = np.asarray(chrom, dtype=str)
    pos_arr = np.asarray(pos, dtype=np.int64)
    if n_alt_arr.ndim != 1:
        raise ValueError("n_alt must be one-dimensional")
    if chrom_arr.ndim != 1 or pos_arr.ndim != 1:
        raise ValueError("chrom and pos must be one-dimensional")
    if chrom_arr.size != n_alt_arr.size or pos_arr.size != n_alt_arr.size:
        raise ValueError("n_alt, chrom, and pos must have matching lengths")
    if not source_freqs:
        raise ValueError("at least one source population is required")

    frequencies = {name: np.asarray(freq, dtype=float) for name, freq in source_freqs.items()}
    if any(freq.ndim != 1 or freq.size != n_alt_arr.size for freq in frequencies.values()):
        raise ValueError("source frequency vectors must match n_alt length")

    calls: list[dict[str, object]] = []
    for contig in dict.fromkeys(chrom_arr.tolist()):
        on_contig = chrom_arr == contig
        contig_freqs = {name: freq[on_contig] for name, freq in frequencies.items()}
        sources, windows, emissions = window_emissions(
            n_alt_arr[on_contig],
            chrom_arr[on_contig],
            pos_arr[on_contig],
            contig_freqs,
            window_bp=window_bp,
            min_sites=min_sites,
        )
        if not windows:
            continue
        path = viterbi(emissions, switch_prob)
        calls.extend(
            {"chrom": w[0], "start": w[1], "end": w[2], "ancestry": sources[state]}
            for w, state in zip(windows, path, strict=True)
        )
    return calls
