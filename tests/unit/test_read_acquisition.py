"""Streaming raw-read acquisition safeguards and layout-preserving downsampling."""

from __future__ import annotations

import gzip
import hashlib
import io
import random
from pathlib import Path

import pytest

from canidae.stages.acquisition import reads
from canidae.stages.acquisition.reads import (
    ReadAcquisitionError,
    ReadArchive,
    download_read_archive,
    enforce_read_transfer_policy,
    estimate_archives,
)


def _gz_fastq(headers: list[str]) -> bytes:
    data = io.BytesIO()
    with gzip.GzipFile(fileobj=data, mode="wb") as handle:
        for header in headers:
            handle.write(f"@{header}\nACGT\n+\n!!!!\n".encode())
    return data.getvalue()


def test_archive_layout_validation_and_transfer_policy() -> None:
    paired = ReadArchive(
        sample_id="wolf", layout="paired", urls=["https://x/r1", "https://x/r2"],
        expected_bytes=[300, 400],
    )
    estimates = estimate_archives([paired])
    assert estimates[0].compressed_bytes == 700
    assert enforce_read_transfer_policy(
        estimates, max_download_bytes=900, confirmation_threshold_bytes=500, confirmed=True
    ) == 700
    with pytest.raises(ReadAcquisitionError, match="explicit confirmation"):
        enforce_read_transfer_policy(
            estimates, max_download_bytes=900, confirmation_threshold_bytes=500, confirmed=False
        )
    with pytest.raises(ValueError, match="exactly 1 URL"):
        ReadArchive(sample_id="bad", layout="single", urls=["a", "b"])


def test_single_and_interleaved_streaming_verify_compressed_checksum(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    single_bytes = _gz_fastq(["r1", "r2"])
    interleaved_bytes = _gz_fastq(["pair1/1", "pair1/2", "pair2/1", "pair2/2"])
    payloads = {"single": single_bytes, "interleaved": interleaved_bytes}

    def fake_open(url: str, _timeout: int) -> io.BytesIO:
        return io.BytesIO(payloads[url])

    monkeypatch.setattr(reads, "_open_stream", fake_open)
    single = ReadArchive(
        sample_id="s", layout="single", urls=["single"],
        checksums=["sha256:" + hashlib.sha256(single_bytes).hexdigest()],
    )
    outputs, _records = download_read_archive(
        single, tmp_path, rng=random.Random(1), timeout_seconds=1
    )
    with gzip.open(outputs[0], "rt") as handle:
        assert handle.read().count("@r") == 2

    interleaved = ReadArchive(sample_id="i", layout="interleaved", urls=["interleaved"],
                              downsample_fraction=0.5)
    outputs, _records = download_read_archive(
        interleaved, tmp_path, rng=random.Random(1), timeout_seconds=1
    )
    with gzip.open(outputs[0], "rt") as handle:
        text = handle.read()
    assert "@pair1/1" in text and "@pair1/2" in text
    assert "@pair2/1" not in text and "@pair2/2" not in text


def test_paired_downsampling_keeps_mates_and_cleans_failed_partial(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payloads = {"r1": _gz_fastq(["pair1/1", "pair2/1"]), "r2": _gz_fastq(["pair1/2", "pair2/2"])}

    def fake_open(url: str, _timeout: int) -> io.BytesIO:
        return io.BytesIO(payloads[url])

    monkeypatch.setattr(reads, "_open_stream", fake_open)
    paired = ReadArchive(
        sample_id="p", layout="paired", urls=["r1", "r2"], downsample_fraction=0.5
    )
    outputs, _records = download_read_archive(
        paired, tmp_path, rng=random.Random(1), timeout_seconds=1
    )
    with gzip.open(outputs[0], "rt") as r1, gzip.open(outputs[1], "rt") as r2:
        assert "@pair1/1" in r1.read()
        assert "@pair1/2" in r2.read()

    mismatch = ReadArchive(
        sample_id="bad", layout="single", urls=["r1"], checksums=["md5:deadbeef"]
    )
    with pytest.raises(ReadAcquisitionError, match="checksum mismatch"):
        download_read_archive(mismatch, tmp_path, rng=random.Random(1), timeout_seconds=1)
    assert not (tmp_path / "bad_single.fastq.gz").exists()
    assert not (tmp_path / "bad_single.fastq.gz.part").exists()
