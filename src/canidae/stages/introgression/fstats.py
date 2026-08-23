"""Core allele-frequency statistics for detecting introgression.

Pure-NumPy implementations of the standard four-population statistics — Patterson's D
(ABBA-BABA), f3, f4, the f4-ratio, and the windowed f_d admixture scan — together with a
block jackknife for standard errors and Z-scores. The jackknife can use legacy
equal-site-count blocks, whole chromosomes, or fixed genomic blocks. Coordinate-aware
calculations visit one chromosome at a time so temporary statistic arrays remain bounded by
the largest chromosome rather than the full callset.

References: Green et al. 2010; Durand et al. 2011 (D with frequencies); Patterson et al.
2012 (f-statistics); Martin, Davey & Jiggins 2015 (f_d).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np

from canidae.stages.popgen.store import Genotypes

BlockMode = Literal["site_count", "chromosome", "fixed_bp"]
_TermFunction = Callable[..., tuple[np.ndarray, np.ndarray]]


@dataclass(slots=True)
class JackknifeResult:
    estimate: float
    se: float
    z: float
    n_sites: int
    n_blocks: int = 0


def allele_frequencies(geno: Genotypes, groups: dict[str, list[int]]) -> dict[str, np.ndarray]:
    """Per-population alt-allele frequency at each site (NaN where a population has no
    called genotypes)."""
    freqs: dict[str, np.ndarray] = {}
    for pop, idx in groups.items():
        ac = geno.calls.count_alleles(subpop=idx)
        totals = ac.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.where(totals > 0, ac[:, 1] / totals, np.nan)
        freqs[pop] = f.astype(float)
    return freqs


def _finite_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(arrays[0].shape[0], dtype=bool)
    for a in arrays:
        mask &= np.isfinite(a)
    return mask


def _abba_baba(
    p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, po: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    abba = (1 - p1) * p2 * p3 * (1 - po) + p1 * (1 - p2) * (1 - p3) * po
    baba = p1 * (1 - p2) * p3 * (1 - po) + (1 - p1) * p2 * (1 - p3) * po
    return abba, baba


def d_statistic(
    p1: np.ndarray,
    p2: np.ndarray,
    p3: np.ndarray,
    po: np.ndarray,
    *,
    n_blocks: int = 20,
    chrom: np.ndarray | None = None,
    pos: np.ndarray | None = None,
    block_mode: BlockMode = "site_count",
    block_size_bp: int = 1_000_000,
) -> JackknifeResult:
    """Patterson's D for the topology ((P1,P2),P3),O. D>0 => excess allele sharing between
    P2 and P3 (gene flow); D<0 => between P1 and P3.

    Coordinate-aware modes require one chromosome label per site. Fixed blocks also require
    1-based positions.
    """
    return _ratio_with_blocks(
        (p1, p2, p3, po),
        _d_terms,
        n_blocks=n_blocks,
        chrom=chrom,
        pos=pos,
        block_mode=block_mode,
        block_size_bp=block_size_bp,
    )


def f4(
    pa: np.ndarray,
    pb: np.ndarray,
    pc: np.ndarray,
    pd: np.ndarray,
    *,
    n_blocks: int = 20,
    chrom: np.ndarray | None = None,
    pos: np.ndarray | None = None,
    block_mode: BlockMode = "site_count",
    block_size_bp: int = 1_000_000,
) -> JackknifeResult:
    """f4(A,B;C,D) = <(a-b)(c-d)>. Zero under a tree with no gene flow."""
    return _ratio_with_blocks(
        (pa, pb, pc, pd),
        lambda a, b, c, d: ((a - b) * (c - d), np.ones_like(a, dtype=float)),
        n_blocks=n_blocks,
        chrom=chrom,
        pos=pos,
        block_mode=block_mode,
        block_size_bp=block_size_bp,
    )


def f3(
    pa: np.ndarray,
    pb: np.ndarray,
    pc: np.ndarray,
    *,
    n_blocks: int = 20,
    chrom: np.ndarray | None = None,
    pos: np.ndarray | None = None,
    block_mode: BlockMode = "site_count",
    block_size_bp: int = 1_000_000,
) -> JackknifeResult:
    """f3(A;B,C) = <(a-b)(a-c)>. Significantly negative => A is admixed between B and C."""
    return _ratio_with_blocks(
        (pa, pb, pc),
        lambda a, b, c: ((a - b) * (a - c), np.ones_like(a, dtype=float)),
        n_blocks=n_blocks,
        chrom=chrom,
        pos=pos,
        block_mode=block_mode,
        block_size_bp=block_size_bp,
    )


def outgroup_f3(
    po: np.ndarray,
    pa: np.ndarray,
    pb: np.ndarray,
    *,
    n_blocks: int = 20,
    chrom: np.ndarray | None = None,
    pos: np.ndarray | None = None,
    block_mode: BlockMode = "site_count",
    block_size_bp: int = 1_000_000,
) -> JackknifeResult:
    """Outgroup f3(O;A,B) = <(o-a)(o-b)>: shared genetic drift of A and B since O."""
    return _ratio_with_blocks(
        (po, pa, pb),
        lambda o, a, b: ((o - a) * (o - b), np.ones_like(o, dtype=float)),
        n_blocks=n_blocks,
        chrom=chrom,
        pos=pos,
        block_mode=block_mode,
        block_size_bp=block_size_bp,
    )


def f4_ratio(
    pa: np.ndarray,
    pb: np.ndarray,
    pc: np.ndarray,
    po: np.ndarray,
    px: np.ndarray,
    *,
    n_blocks: int = 20,
    chrom: np.ndarray | None = None,
    pos: np.ndarray | None = None,
    block_mode: BlockMode = "site_count",
    block_size_bp: int = 1_000_000,
) -> JackknifeResult:
    """Admixture proportion alpha via the f4-ratio f4(O,A;X,C)/f4(O,A;B,C)."""
    return _ratio_with_blocks(
        (pa, pb, pc, po, px),
        lambda a, b, c, o, x: (
            (o - a) * (x - c),
            (o - a) * (b - c),
        ),
        n_blocks=n_blocks,
        chrom=chrom,
        pos=pos,
        block_mode=block_mode,
        block_size_bp=block_size_bp,
    )


def f_d_windows(
    freqs: dict[str, np.ndarray],
    quartet: tuple[str, str, str, str],
    chrom: np.ndarray,
    pos: np.ndarray,
    *,
    window_bp: int = 100_000,
    min_sites: int = 10,
) -> list[dict[str, object]]:
    """Windowed f_d (Martin et al. 2015): a localized, sign-consistent introgression
    estimate that (unlike D) is bounded and comparable across windows."""
    p1, p2, p3, po = (np.asarray(freqs[q], dtype=float) for q in quartet)
    vectors = _validate_vectors((p1, p2, p3, po))
    chrom_arr = np.asarray(chrom, dtype=str)
    pos_arr = np.asarray(pos, dtype=np.int64)
    _validate_coordinates(chrom_arr, pos_arr, vectors[0].size)

    out: list[dict[str, object]] = []
    for contig in dict.fromkeys(chrom_arr.tolist()):
        cm = chrom_arr == contig
        cp = pos_arr[cm]
        c1, c2, c3, co = (v[cm] for v in vectors)
        finite = _finite_mask(c1, c2, c3, co)
        cp, i1, i2, i3, io = cp[finite], c1[finite], c2[finite], c3[finite], co[finite]
        if cp.size == 0:
            continue
        start = int(cp.min())
        for w0 in range(start, int(cp.max()) + 1, window_bp):
            w = (cp >= w0) & (cp < w0 + window_bp)
            if int(w.sum()) < min_sites:
                continue
            out.append(
                {
                    "chrom": str(contig),
                    "start": w0,
                    "end": w0 + window_bp,
                    "n_sites": int(w.sum()),
                    "f_d": _f_d(i1[w], i2[w], i3[w], io[w]),
                }
            )
    return out


def _f_d(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, po: np.ndarray) -> float:
    abba, baba = _abba_baba(p1, p2, p3, po)
    num = float((abba - baba).sum())
    donor = np.maximum(p2, p3)  # per-site dynamic donor
    abba_d, baba_d = _abba_baba(p1, donor, donor, po)
    den = float((abba_d - baba_d).sum())
    return round(num / den, 6) if den != 0 else float("nan")


def _d_terms(
    p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, po: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    abba, baba = _abba_baba(p1, p2, p3, po)
    return abba - baba, abba + baba


def _ratio_with_blocks(
    arrays: tuple[np.ndarray, ...],
    terms: _TermFunction,
    *,
    n_blocks: int,
    chrom: np.ndarray | None,
    pos: np.ndarray | None,
    block_mode: BlockMode,
    block_size_bp: int,
) -> JackknifeResult:
    """Evaluate a ratio statistic and estimate uncertainty with the requested blocks.

    The historical site-count mode retains the original behaviour. Coordinate-aware modes
    deliberately work contig-by-contig, retaining only numerator and denominator sums for
    each jackknife block.
    """
    vectors = _validate_vectors(arrays)
    if block_mode not in {"site_count", "chromosome", "fixed_bp"}:
        raise ValueError(f"unknown block_mode: {block_mode!r}")

    if block_mode == "site_count":
        mask = _finite_mask(*vectors)
        num, den = terms(*(v[mask] for v in vectors))
        return _jackknife_ratio(num, den, n_blocks=n_blocks)

    if chrom is None:
        raise ValueError(f"block_mode={block_mode!r} requires chrom coordinates")
    chrom_arr = np.asarray(chrom, dtype=str)
    pos_arr = None if pos is None else np.asarray(pos, dtype=np.int64)
    if block_mode == "fixed_bp" and pos_arr is None:
        raise ValueError("block_mode='fixed_bp' requires positions")
    _validate_coordinates(chrom_arr, pos_arr, vectors[0].size)
    if block_mode == "fixed_bp" and block_size_bp <= 0:
        raise ValueError("block_size_bp must be positive")

    block_num: list[float] = []
    block_den: list[float] = []
    n_sites = 0
    for contig in dict.fromkeys(chrom_arr.tolist()):
        on_contig = chrom_arr == contig
        contig_arrays = tuple(v[on_contig] for v in vectors)
        finite = _finite_mask(*contig_arrays)
        if not np.any(finite):
            continue
        num, den = terms(*(v[finite] for v in contig_arrays))
        n_sites += int(num.size)

        if block_mode == "chromosome":
            block_num.append(float(num.sum()))
            block_den.append(float(den.sum()))
            continue

        assert pos_arr is not None
        block_index = (pos_arr[on_contig][finite] - 1) // block_size_bp
        for index in np.unique(block_index):
            in_block = block_index == index
            block_num.append(float(num[in_block].sum()))
            block_den.append(float(den[in_block].sum()))

    return _jackknife_from_block_sums(block_num, block_den, n_sites=n_sites)


def _validate_vectors(arrays: tuple[np.ndarray, ...]) -> tuple[np.ndarray, ...]:
    vectors = tuple(np.asarray(array, dtype=float) for array in arrays)
    if not vectors:
        raise ValueError("at least one statistic vector is required")
    length = vectors[0].size
    if any(vector.ndim != 1 or vector.size != length for vector in vectors):
        raise ValueError("statistic vectors must be one-dimensional and have matching lengths")
    return vectors


def _validate_coordinates(
    chrom: np.ndarray,
    pos: np.ndarray | None,
    length: int,
) -> None:
    if chrom.ndim != 1 or chrom.size != length:
        raise ValueError("chrom must be one-dimensional and match statistic vector length")
    if pos is not None:
        if pos.ndim != 1 or pos.size != length:
            raise ValueError("pos must be one-dimensional and match statistic vector length")
        if np.any(pos <= 0):
            raise ValueError("positions must be positive and 1-based")


def _jackknife_ratio(num: np.ndarray, den: np.ndarray, *, n_blocks: int) -> JackknifeResult:
    """Delete-one block jackknife for a ratio-of-sums estimator theta = Σnum / Σden."""
    n = num.shape[0]
    if n_blocks < 1:
        raise ValueError("n_blocks must be >= 1")
    blocks = np.array_split(np.arange(n), min(n_blocks, n)) if n else []
    return _jackknife_from_block_sums(
        [float(num[block].sum()) for block in blocks],
        [float(den[block].sum()) for block in blocks],
        n_sites=n,
    )


def _jackknife_from_block_sums(
    block_num: list[float],
    block_den: list[float],
    *,
    n_sites: int,
) -> JackknifeResult:
    """Delete one chromosome/fixed block at a time from a ratio-of-sums statistic."""
    numerator: np.ndarray = np.asarray(block_num, dtype=float)
    denominator: np.ndarray = np.asarray(block_den, dtype=float)
    n_blocks = int(numerator.size)
    total_num = float(numerator.sum())
    total_den = float(denominator.sum())
    if n_sites == 0 or n_blocks == 0 or total_den == 0:
        return JackknifeResult(float("nan"), float("nan"), float("nan"), n_sites, n_blocks)
    theta = total_num / total_den

    partials = []
    for num_block, den_block in zip(numerator, denominator, strict=True):
        d_left = total_den - den_block
        if d_left == 0:
            continue
        partials.append((total_num - num_block) / d_left)
    thetas = np.asarray(partials, dtype=float)
    g = thetas.size
    if g < 2:
        return JackknifeResult(theta, float("nan"), float("nan"), n_sites, n_blocks)
    var = (g - 1) / g * np.sum((thetas - thetas.mean()) ** 2)
    se = float(np.sqrt(var))
    z = float(theta / se) if se > 0 else float("nan")
    return JackknifeResult(theta, se, z, n_sites, n_blocks)
