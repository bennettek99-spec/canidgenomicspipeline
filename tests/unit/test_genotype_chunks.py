"""Tests for the chunk-iterating genotype API over the memory-mapped store.

The contract has two halves: chunked results must equal the whole-matrix
results exactly, and the memory-mapped path must not materialize the whole
genotype matrix to produce them.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import allel
import numpy as np
import pytest

from canidae.stages.popgen.store import (
    DEFAULT_CHUNK_VARIANTS,
    Genotypes,
    chunked_allele_counts,
    chunked_alt_frequency,
    iter_genotype_chunks,
    load_genotypes,
    save_genotypes,
)


def _cohort(n_variants: int = 250, n_samples: int = 8, seed: int = 3) -> Genotypes:
    rng = np.random.default_rng(seed)
    gt = rng.integers(0, 2, size=(n_variants, n_samples, 2), dtype="i1")
    # A fully uncalled site, so the missing-data paths are covered.
    gt[7] = -1
    return Genotypes(
        calls=allel.GenotypeArray(gt),
        pos=np.arange(1, n_variants + 1, dtype=np.int64) * 100,
        chrom=np.array(["1"] * n_variants),
        samples=np.array([f"S{i:02d}" for i in range(n_samples)]),
    )


@pytest.fixture()
def mmap_store(tmp_path: Path) -> tuple[Path, Genotypes]:
    geno = _cohort()
    path = save_genotypes(tmp_path / "store", geno, backend="npy_mmap")
    return path, geno


def test_chunks_tile_the_variant_axis_exactly(mmap_store) -> None:
    path, geno = mmap_store
    chunks = list(iter_genotype_chunks(path, chunk_variants=60))
    assert [c.n_variants for c in chunks] == [60, 60, 60, 60, 10]
    assert chunks[0].start == 0
    assert chunks[-1].stop == geno.n_variants
    # Boundaries meet with no gap and no overlap.
    for previous, following in itertools.pairwise(chunks):
        assert previous.stop == following.start


def test_chunk_contents_match_the_whole_matrix(mmap_store) -> None:
    path, geno = mmap_store
    rebuilt = np.concatenate(
        [np.asarray(c.calls) for c in iter_genotype_chunks(path, chunk_variants=37)]
    )
    assert np.array_equal(rebuilt, np.asarray(geno.calls))

    positions = np.concatenate([c.pos for c in iter_genotype_chunks(path, chunk_variants=37)])
    assert np.array_equal(positions, geno.pos)
    for chunk in iter_genotype_chunks(path, chunk_variants=37):
        assert chunk.samples.tolist() == geno.samples.tolist()
        assert chunk.chrom.tolist() == ["1"] * chunk.n_variants


def test_iteration_does_not_materialize_the_whole_matrix(mmap_store, monkeypatch) -> None:
    """The mmap backend must stream; a chunk may not be larger than requested."""
    path, geno = mmap_store
    real_load = np.load
    full_loads: list[str] = []

    def spy(file, *args, **kwargs):
        if kwargs.get("mmap_mode") is None and str(file).endswith("gt.npy"):
            full_loads.append(str(file))
        return real_load(file, *args, **kwargs)

    monkeypatch.setattr(np, "load", spy)
    chunks = list(iter_genotype_chunks(path, chunk_variants=25))
    assert not full_loads, "genotype matrix was read without mmap_mode"
    assert max(c.n_variants for c in chunks) == 25
    assert sum(c.n_variants for c in chunks) == geno.n_variants


def test_chunked_allele_counts_match_the_whole_matrix(mmap_store) -> None:
    path, geno = mmap_store
    expected = np.asarray(geno.allele_counts())
    observed = np.asarray(chunked_allele_counts(path, chunk_variants=33))
    assert observed.shape == expected.shape
    assert np.array_equal(observed, expected)


def test_chunked_allele_counts_honour_a_subpopulation(mmap_store) -> None:
    path, geno = mmap_store
    subpop = [0, 2, 4]
    expected = np.asarray(geno.allele_counts(subpop=subpop))
    observed = np.asarray(chunked_allele_counts(path, subpop=subpop, chunk_variants=40))
    assert np.array_equal(observed, expected)


def test_chunked_alt_frequency_matches_a_direct_computation(mmap_store) -> None:
    path, geno = mmap_store
    counts = np.asarray(geno.allele_counts())
    total = counts.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        expected = np.where(total > 0, counts[:, 1:].sum(axis=1) / total, np.nan)
    observed = chunked_alt_frequency(path, chunk_variants=64)
    assert np.allclose(observed, expected, equal_nan=True)
    # The all-missing site has no called alleles and must stay NaN.
    assert np.isnan(observed[7])


def test_npz_archives_iterate_with_the_same_contract(tmp_path: Path) -> None:
    geno = _cohort(n_variants=120)
    path = save_genotypes(tmp_path / "g.npz", geno)
    chunks = list(iter_genotype_chunks(path, chunk_variants=50))
    assert [c.n_variants for c in chunks] == [50, 50, 20]
    assert np.array_equal(
        np.asarray(chunked_allele_counts(path, chunk_variants=50)),
        np.asarray(load_genotypes(path).allele_counts()),
    )


def test_an_in_memory_cohort_can_be_chunked_directly() -> None:
    geno = _cohort(n_variants=30)
    chunks = list(iter_genotype_chunks(geno, chunk_variants=12))
    assert [c.n_variants for c in chunks] == [12, 12, 6]
    assert np.array_equal(
        np.asarray(chunked_allele_counts(geno, chunk_variants=12)),
        np.asarray(geno.allele_counts()),
    )


def test_a_single_chunk_covers_a_cohort_smaller_than_the_chunk_size(mmap_store) -> None:
    path, geno = mmap_store
    chunks = list(iter_genotype_chunks(path, chunk_variants=DEFAULT_CHUNK_VARIANTS))
    assert len(chunks) == 1
    assert chunks[0].n_variants == geno.n_variants


def test_a_non_positive_chunk_size_is_rejected(mmap_store) -> None:
    path, _ = mmap_store
    with pytest.raises(ValueError, match="must be positive"):
        list(iter_genotype_chunks(path, chunk_variants=0))
