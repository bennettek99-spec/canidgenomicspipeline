from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from canidae.core.errors import StageInputError
from canidae.stages.processing.harmonize import _orientation
from canidae.stages.processing.vcf_io import read_biallelic_snps, write_minimal_vcf
from canidae.stages.processing.wgs import (
    align_command,
    bcftools_joint_command,
    call_command,
    joint_command,
    read_fastq_manifest,
)


def test_vcf_write_read_roundtrip(tmp_path: Path) -> None:
    chrom = np.array(["1", "1", "1"])
    pos = np.array([100, 200, 300], dtype=np.int64)
    ref = np.array(["A", "C", "G"])
    alt = np.array(["T", "G", "A"])
    samples = np.array(["s1", "s2"])
    gt = np.array([[[0, 0], [0, 1]], [[1, 1], [-1, -1]], [[0, 1], [1, 1]]], dtype="i1")
    path = write_minimal_vcf(tmp_path / "x.vcf", chrom, pos, ref, alt, samples, gt)
    back = read_biallelic_snps(path)
    assert list(back.samples) == ["s1", "s2"]
    assert back.pos.tolist() == [100, 200, 300]
    np.testing.assert_array_equal(back.gt, gt)


def test_align_command_pipeline() -> None:
    cmd = align_command("W1", "r1.fq.gz", "r2.fq.gz", "ref.fa", "W1.cram", 8)
    assert "bwa-mem2 mem -t 8" in cmd
    assert "SM:W1" in cmd
    assert "samtools markdup" in cmd and "W1.cram" in cmd


def test_align_command_supports_interleaved_and_single_end_reads() -> None:
    interleaved = align_command(
        "I1", "reads.fq.gz", "", "ref.fa", "I1.cram", 2, layout="interleaved"
    )
    single = align_command("S1", "reads.fq.gz", "", "ref.fa", "S1.cram", 2, layout="single")
    assert " -p reads.fq.gz " in interleaved
    assert "reads.fq.gz " in single and " -p " not in single


def test_call_command_backends() -> None:
    assert "HaplotypeCaller" in call_command("s", "s.cram", "ref", "s.g.vcf", "gatk", 4)
    assert "mpileup" in call_command("s", "s.cram", "ref", "s.g.vcf", "bcftools", 4)
    with pytest.raises(StageInputError):
        call_command("s", "s.cram", "ref", "s.g.vcf", "nonsense", 4)


def test_joint_command_backends() -> None:
    assert "glnexus_cli" in joint_command(["a.g.vcf", "b.g.vcf"], "ref", "out.vcf", "glnexus", 4)
    assert "GenotypeGVCFs" in joint_command(["a.g.vcf"], "ref", "out.vcf", "gatk", 4)


def test_bcftools_joint_command() -> None:
    cmd = bcftools_joint_command(["A.cram", "B.cram", "C.cram"], "ref.fa", "joint.vcf.gz", 8)
    assert "bcftools mpileup -f ref.fa" in cmd
    assert "A.cram B.cram C.cram" in cmd  # all samples in one call
    assert "bcftools call -m -v" in cmd and "joint.vcf.gz" in cmd


def test_fastq_manifest_validation(tmp_path: Path) -> None:
    bad = tmp_path / "m.csv"
    bad.write_text("sample_id,fastq1\nW1,r1.fq\n", encoding="utf-8")
    with pytest.raises(StageInputError):
        read_fastq_manifest(bad)


def test_harmonization_checks_swap_complement_and_palindromic_orientation() -> None:
    assert _orientation("A", "C", "A", "C", True) == ("same", False)
    assert _orientation("A", "C", "C", "A", True) == ("swap", True)
    assert _orientation("A", "C", "T", "G", True) == ("complement", False)
    # Strand-only A/T and C/G matching is deliberately not guessed on a laptop.
    assert _orientation("A", "T", "T", "A", True) == ("swap", True)
    assert _orientation("A", "T", "A", "T", True) == ("same", False)
