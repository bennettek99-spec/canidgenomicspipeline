"""Materialized genotype store shared by the population-genomics stages.

Rather than re-parse the VCF in every analysis, one ``load_genotypes`` stage reads it once
into a compact on-disk matrix; PCA, F_ST, and diversity all read that matrix back. This is
the small-scale stand-in for the chunked Zarr/sgkit store used at the thousands-of-genomes
scale — same contract (load once, analyse many), different backing format.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import allel
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

# Variants per block for the chunk-iterating API. Sized so one chunk of a
# few-hundred-sample cohort stays in the low tens of megabytes.
DEFAULT_CHUNK_VARIANTS = 50_000


@dataclass(slots=True)
class Genotypes:
    """An in-memory genotype matrix plus its coordinate and sample context."""

    calls: allel.GenotypeArray  # (n_variants, n_samples, ploidy)
    pos: np.ndarray  # (n_variants,) int
    chrom: np.ndarray  # (n_variants,) str
    samples: np.ndarray  # (n_samples,) str

    @property
    def n_variants(self) -> int:
        return self.calls.shape[0]

    @property
    def n_samples(self) -> int:
        return self.calls.shape[1]

    def allele_counts(self, subpop: list[int] | None = None) -> allel.AlleleCountsArray:
        return self.calls.count_alleles(subpop=subpop)

    def sample_index(self) -> dict[str, int]:
        return {s: i for i, s in enumerate(self.samples)}


def save_genotypes(path: Path, genotypes: Genotypes, *, backend: str = "npz") -> Path:
    """Persist genotypes as compact NPZ or an out-of-core memory-mapped directory."""
    path = Path(path)
    if backend == "npy_mmap":
        path.mkdir(parents=True, exist_ok=True)
        np.save(path / "gt.npy", np.asarray(genotypes.calls, dtype=np.int8), allow_pickle=False)
        np.save(path / "pos.npy", np.asarray(genotypes.pos, dtype=np.int64), allow_pickle=False)
        np.save(path / "chrom.npy", np.asarray(genotypes.chrom, dtype="U32"), allow_pickle=False)
        np.save(
            path / "samples.npy", np.asarray(genotypes.samples, dtype="U64"), allow_pickle=False
        )
        (path / "metadata.json").write_text(
            json.dumps({"schema_version": 1, "backend": "npy_mmap"}, indent=2) + "\n",
            encoding="utf-8",
        )
        return path
    if backend != "npz":
        raise ValueError(f"unsupported genotype storage backend: {backend}")
    np.savez_compressed(
        path,
        gt=np.asarray(genotypes.calls, dtype=np.int8),
        pos=np.asarray(genotypes.pos, dtype=np.int64),
        chrom=np.asarray(genotypes.chrom, dtype="U32"),
        samples=np.asarray(genotypes.samples, dtype="U64"),
    )
    # np.savez appends .npz if absent; normalize the returned path to the real file.
    return path if path.suffix == ".npz" else path.with_suffix(".npz")


def load_genotypes(path: Path) -> Genotypes:
    """Load a :class:`Genotypes` previously written by :func:`save_genotypes`."""
    path = Path(path)
    if path.is_dir():
        return Genotypes(
            calls=allel.GenotypeArray(np.load(path / "gt.npy", mmap_mode="r")),
            pos=np.load(path / "pos.npy", mmap_mode="r"),
            chrom=np.load(path / "chrom.npy", mmap_mode="r").astype(str),
            samples=np.load(path / "samples.npy", mmap_mode="r").astype(str),
        )
    with np.load(path, allow_pickle=False) as data:
        return Genotypes(
            calls=allel.GenotypeArray(data["gt"]),
            pos=data["pos"],
            chrom=data["chrom"].astype(str),
            samples=data["samples"].astype(str),
        )


@dataclass(slots=True)
class GenotypeChunk:
    """A contiguous block of variants, with the coordinates that describe it.

    ``start``/``stop`` are indices into the full variant axis, so a caller can
    write results back into a whole-cohort array without tracking its own offset.
    """

    start: int
    stop: int
    calls: allel.GenotypeArray  # (stop - start, n_samples, ploidy)
    pos: np.ndarray
    chrom: np.ndarray
    samples: np.ndarray

    @property
    def n_variants(self) -> int:
        return self.stop - self.start


def iter_genotype_chunks(
    source: Path | str | Genotypes,
    *,
    chunk_variants: int = DEFAULT_CHUNK_VARIANTS,
) -> Iterator[GenotypeChunk]:
    """Yield blocks of variants without holding the whole matrix in memory.

    For the ``npy_mmap`` backend each block is read from the memory map on
    demand, so peak resident memory scales with ``chunk_variants`` rather than
    with the cohort. NPZ archives cannot be memory-mapped, so that path decodes
    once and then slices; the iteration contract is identical either way.

    This is for algorithms that only need a variant window at a time — allele
    counts, per-site diversity, windowed scans. Whole-matrix methods (PCA,
    distance, NMF admixture) still use :func:`load_genotypes`.
    """
    if chunk_variants < 1:
        raise ValueError(f"chunk_variants must be positive, got {chunk_variants}")

    if isinstance(source, Genotypes):
        gt, pos, chrom, samples = (source.calls, source.pos, source.chrom, source.samples)
    else:
        path = Path(source)
        if path.is_dir():
            gt = np.load(path / "gt.npy", mmap_mode="r")
            pos = np.load(path / "pos.npy", mmap_mode="r")
            chrom = np.load(path / "chrom.npy", mmap_mode="r")
            samples = np.load(path / "samples.npy", mmap_mode="r").astype(str)
        else:
            genotypes = load_genotypes(path)
            gt, pos, chrom, samples = (
                genotypes.calls,
                genotypes.pos,
                genotypes.chrom,
                genotypes.samples,
            )

    n_variants = gt.shape[0]
    for start in range(0, n_variants, chunk_variants):
        stop = min(start + chunk_variants, n_variants)
        yield GenotypeChunk(
            start=start,
            stop=stop,
            calls=allel.GenotypeArray(np.asarray(gt[start:stop])),
            pos=np.asarray(pos[start:stop]),
            chrom=np.asarray(chrom[start:stop]).astype(str),
            samples=np.asarray(samples).astype(str),
        )


def chunked_allele_counts(
    source: Path | str | Genotypes,
    *,
    subpop: list[int] | None = None,
    chunk_variants: int = DEFAULT_CHUNK_VARIANTS,
) -> allel.AlleleCountsArray:
    """Allele counts over all variants, accumulated one chunk at a time.

    Equivalent to ``load_genotypes(path).allele_counts(subpop)`` but never holds
    more than ``chunk_variants`` rows of the genotype matrix at once.
    """
    blocks = [
        np.asarray(chunk.calls.count_alleles(subpop=subpop))
        for chunk in iter_genotype_chunks(source, chunk_variants=chunk_variants)
    ]
    if not blocks:
        return allel.AlleleCountsArray(np.empty((0, 2), dtype="i4"))
    width = max(block.shape[1] for block in blocks)
    padded = [
        block if block.shape[1] == width else np.pad(block, ((0, 0), (0, width - block.shape[1])))
        for block in blocks
    ]
    return allel.AlleleCountsArray(np.concatenate(padded, axis=0))


def chunked_alt_frequency(
    source: Path | str | Genotypes,
    *,
    subpop: list[int] | None = None,
    chunk_variants: int = DEFAULT_CHUNK_VARIANTS,
) -> np.ndarray:
    """Per-variant ALT allele frequency; NaN where no allele was called."""
    counts = np.asarray(chunked_allele_counts(source, subpop=subpop, chunk_variants=chunk_variants))
    if counts.size == 0:
        return np.empty(0, dtype=float)
    total = counts.sum(axis=1)
    alt = counts[:, 1:].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        freqs = np.where(total > 0, alt / total, np.nan)
    return freqs.astype(float)


def allele_difference_matrix(gn: np.ndarray) -> np.ndarray:
    """Pairwise allele-difference distance in [0, 1] from an ALT-count matrix.

    ``gn`` is (n_sites, n_samples); returns an (n_samples, n_samples) matrix equal to the
    mean per-site absolute ALT-count difference divided by 2.
    """
    X = np.asarray(gn, dtype=float).T
    n_sites = max(X.shape[1], 1)
    return squareform(pdist(X, metric="cityblock") / (2.0 * n_sites))


def load_sample_labels(sample_sheet_csv: Path) -> pd.DataFrame:
    """Read the normalized sample sheet, indexed by ``sample_id``.

    Guaranteed columns: ``taxon``, ``population`` (plus any optional metadata columns).
    """
    df = pd.read_csv(sample_sheet_csv, dtype=str)
    return df.set_index("sample_id")


def align_labels(genotypes: Genotypes, labels: pd.DataFrame) -> pd.DataFrame:
    """Return per-sample labels aligned to the genotype-matrix sample order.

    Samples absent from the sample sheet are filled with ``"unknown"`` so downstream code
    never key-errors on an unlabeled sample.
    """
    rows = []
    for sample_id in genotypes.samples:
        if sample_id in labels.index:
            row = labels.loc[sample_id]
            rows.append(
                {
                    "sample_id": sample_id,
                    "taxon": str(row.get("taxon", "unknown")),
                    "population": str(row.get("population", "unknown")),
                }
            )
        else:
            rows.append({"sample_id": sample_id, "taxon": "unknown", "population": "unknown"})
    return pd.DataFrame(rows)


def population_indices(
    genotypes: Genotypes, labels: pd.DataFrame, *, column: str = "population"
) -> dict[str, list[int]]:
    """Map each population label to the genotype-matrix column indices of its samples.

    Only samples present in *both* the genotype matrix and the sample sheet are included;
    the mapping is keyed to the genotype sample order so allele counts align.
    """
    idx = genotypes.sample_index()
    groups: dict[str, list[int]] = {}
    for sample_id, row in labels.iterrows():
        col = idx.get(str(sample_id))
        if col is None:
            continue
        groups.setdefault(str(row[column]), []).append(col)
    return {k: sorted(v) for k, v in sorted(groups.items())}
