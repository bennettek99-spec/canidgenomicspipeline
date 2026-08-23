"""Property tests for genotype encodings and small file-format boundaries."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import allel
import numpy as np
from hypothesis import given
from hypothesis import strategies as st

from canidae.analysis.genotypes import allele_count, dosage, parse_gt_field
from canidae.stages.popgen.plink import decode_bed, write_plink_bed
from canidae.stages.popgen.store import Genotypes
from canidae.stages.processing.vcf_io import read_biallelic_snps, write_minimal_vcf


def _gt(count: int) -> list[int]:
    return {
        -1: [-1, -1],
        0: [0, 0],
        1: [0, 1],
        2: [1, 1],
    }[count]


@given(st.lists(st.integers(min_value=0, max_value=2), min_size=1, max_size=64))
def test_string_genotype_scales_are_inverse(counts: list[int]) -> None:
    strings = [("0/0", "0/1", "1/1")[count] for count in counts]
    assert [allele_count(gt) for gt in strings] == counts
    assert [dosage(gt) for gt in strings] == [count / 2 for count in counts]
    assert all(parse_gt_field(f"{gt}:40", 0) == gt for gt in strings)


@given(st.lists(st.integers(min_value=-1, max_value=2), min_size=9, max_size=9))
def test_vcf_round_trip_preserves_alt_counts(values: list[int]) -> None:
    calls: np.ndarray = np.asarray([_gt(value) for value in values], dtype="i1").reshape(3, 3, 2)
    with TemporaryDirectory() as directory:
        path = write_minimal_vcf(
            Path(directory) / "roundtrip.vcf",
            np.array(["1", "1", "2"]),
            np.array([100, 200, 100]),
            np.array(["A", "C", "G"]),
            np.array(["G", "T", "A"]),
            np.array(["S1", "S2", "S3"]),
            calls,
        )
        observed = read_biallelic_snps(path)
    np.testing.assert_array_equal(observed.gt, calls)


@given(st.lists(st.integers(min_value=-1, max_value=2), min_size=20, max_size=20))
def test_plink_round_trip_preserves_alt_counts(values: list[int]) -> None:
    calls: np.ndarray = np.asarray([_gt(value) for value in values], dtype="i1").reshape(4, 5, 2)
    geno = Genotypes(
        calls=allel.GenotypeArray(calls),
        pos=np.arange(100, 500, 100, dtype=np.int64),
        chrom=np.array(["1", "1", "1", "2"]),
        samples=np.array(["S1", "S2", "S3", "S4", "S5"]),
    )
    with TemporaryDirectory() as directory:
        fileset = write_plink_bed(geno, Path(directory) / "roundtrip")
        expected: np.ndarray = np.asarray(values, dtype="i1").reshape(4, 5)
        np.testing.assert_array_equal(decode_bed(fileset, n_samples=5), expected)
