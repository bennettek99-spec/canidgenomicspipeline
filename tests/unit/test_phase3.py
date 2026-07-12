from __future__ import annotations

import gzip
from pathlib import Path

import allel
import numpy as np
import pytest

from canidae.stages.introgression.fstats import d_statistic, f4
from canidae.stages.phylogenetics.ml_tree import write_phylip
from canidae.stages.phylogenetics.neighbor_joining import (
    bipartitions,
    neighbor_joining,
    to_newick,
)
from canidae.stages.phylogenetics.treemix import write_treemix_input
from canidae.stages.popgen.store import Genotypes


def _abba_baba_freqs(n_abba: int, n_baba: int, *, seed: int = 0):
    """P1,P2,P3,O frequencies with n_abba ABBA and n_baba BABA sites, interleaved (as in
    real data) so the block jackknife sees homogeneous blocks."""
    n = n_abba + n_baba
    rng = np.random.default_rng(seed)
    is_abba = np.zeros(n, dtype=bool)
    is_abba[rng.permutation(n)[:n_abba]] = True
    p1 = np.where(is_abba, 0.0, 1.0)
    p2 = np.where(is_abba, 1.0, 0.0)
    p3 = np.ones(n)
    po = np.zeros(n)
    return p1, p2, p3, po


def test_d_statistic_detects_abba_excess() -> None:
    d = d_statistic(*_abba_baba_freqs(600, 200), n_blocks=20)
    assert 0.4 < d.estimate < 0.6      # (600-200)/(600+200) = 0.5
    assert abs(d.z) > 3                 # significant


def test_d_statistic_symmetric_is_zero() -> None:
    d = d_statistic(*_abba_baba_freqs(300, 300), n_blocks=20)
    assert abs(d.estimate) < 0.05


def test_f4_nonzero_under_gene_flow() -> None:
    stat = f4(*_abba_baba_freqs(600, 200), n_blocks=20)
    assert abs(stat.estimate) > 0.3


def test_coordinate_aware_jackknife_keeps_chromosome_blocks_distinct() -> None:
    p1, p2, p3, po = _abba_baba_freqs(4, 4)
    chrom = np.array(["1", "1", "1", "1", "2", "2", "2", "2"])
    # The same local coordinates occur on both chromosomes. A correct fixed-block
    # jackknife must still retain four independent blocks, not merge them into two.
    pos = np.array([100, 200, 1_000_100, 1_000_200] * 2)

    fixed = d_statistic(
        p1, p2, p3, po,
        chrom=chrom,
        pos=pos,
        block_mode="fixed_bp",
        block_size_bp=1_000_000,
    )
    by_chromosome = d_statistic(
        p1, p2, p3, po,
        chrom=chrom,
        block_mode="chromosome",
    )
    stat_f4 = f4(
        p1, p2, p3, po,
        chrom=chrom,
        block_mode="chromosome",
    )

    assert fixed.n_blocks == 4
    assert by_chromosome.n_blocks == 2
    assert stat_f4.n_blocks == 2
    assert np.isclose(fixed.estimate, by_chromosome.estimate)


def test_fixed_bp_blocks_require_positions() -> None:
    p1, p2, p3, po = _abba_baba_freqs(4, 4)
    with pytest.raises(ValueError, match="requires positions"):
        d_statistic(
            p1, p2, p3, po,
            chrom=np.array(["1"] * 8),
            block_mode="fixed_bp",
        )


def test_neighbor_joining_recovers_topology() -> None:
    # additive distances for the tree ((A,B),(C,D))
    labels = ["A", "B", "C", "D"]
    d = np.array([
        [0, 2, 3, 3],
        [2, 0, 3, 3],
        [3, 3, 0, 2],
        [3, 3, 2, 0],
    ], dtype=float)
    tree = neighbor_joining(d, labels)
    splits = bipartitions(tree)
    assert frozenset({"C", "D"}) in splits    # the A,B | C,D split (canonicalized vs ref A)


def test_to_newick_is_parseable() -> None:
    from io import StringIO

    from Bio import Phylo

    labels = ["A", "B", "C", "D"]
    d = np.array([[0, 2, 3, 3], [2, 0, 3, 3], [3, 3, 0, 2], [3, 3, 2, 0]], dtype=float)
    newick = to_newick(neighbor_joining(d, labels))
    tree = Phylo.read(StringIO(newick), "newick")
    assert tree.count_terminals() == 4


def _demo_genotypes() -> Genotypes:
    gt = np.array(
        [
            [[0, 0], [0, 1], [1, 1], [-1, -1]],
            [[1, 1], [0, 0], [0, 1], [0, 0]],
        ],
        dtype="i1",
    )
    return Genotypes(
        calls=allel.GenotypeArray(gt),
        pos=np.array([100, 200], dtype=np.int64),
        chrom=np.array(["1", "1"]),
        samples=np.array(["A", "B", "C", "D"]),
    )


def test_write_phylip_iupac(tmp_path: Path) -> None:
    geno = _demo_genotypes()
    path = write_phylip(geno, tmp_path / "aln.phy", segregating_only=False)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "4 2"                 # 4 taxa, 2 sites
    seqs = {ln.split()[0]: ln.split()[1] for ln in lines[1:]}
    assert seqs["B"] == "WA"                 # het then hom-ref
    assert seqs["D"][0] == "N"               # missing genotype -> N


def test_write_treemix_input(tmp_path: Path) -> None:
    geno = _demo_genotypes()
    groups = {"pop1": [0, 1], "pop2": [2, 3]}
    path = write_treemix_input(geno, groups, tmp_path / "tm.gz")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().split()
        first = fh.readline().split()
    assert header == ["pop1", "pop2"]
    assert all("," in field for field in first)   # "ref,alt" per population
