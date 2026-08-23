"""Laptop-safe streamed acquisition of paired, interleaved, and single-end FASTQ archives.

This stage is deliberately conservative. It estimates every archive before transferring data,
refuses a total above its hard sub-10-GB ceiling, requires an explicit confirmation above a
lower threshold, streams gzip FASTQ records directly into the retained output, and preserves
paired/interleaved read coherence during downsampling. Temporary ``.part`` files are removed
on every failure.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import os
import random
import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, BinaryIO, Literal
from urllib.request import Request, urlopen

import pandas as pd
from pydantic import BaseModel, Field, field_validator, model_validator

from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.acquisition.reduced_panel import (
    DEFAULT_CONFIRMATION_THRESHOLD_BYTES,
    DEFAULT_DOWNLOAD_CEILING_BYTES,
    MAX_DOWNLOAD_CEILING_BYTES,
    format_bytes,
)


class ReadAcquisitionError(RuntimeError):
    """A safe acquisition failure with a plain-language recovery instruction."""

    def __init__(self, message: str, *, recovery: str) -> None:
        super().__init__(message)
        self.recovery = recovery


class ReadArchive(BaseModel):
    """One sample's public FASTQ archive layout and publisher checksums."""

    sample_id: str
    layout: Literal["paired", "interleaved", "single"]
    urls: list[str] = Field(min_length=1)
    checksums: list[str] = Field(default_factory=list)
    expected_bytes: list[int] = Field(default_factory=list)
    downsample_fraction: float = Field(default=1.0, gt=0.0, le=1.0)

    @field_validator("sample_id")
    @classmethod
    def _sample_id_not_blank(cls, value: str) -> str:
        value = value.strip()
        if not value:
            raise ValueError("sample_id cannot be blank")
        return value

    @model_validator(mode="after")
    def _validate_layout(self) -> ReadArchive:
        expected_urls = 2 if self.layout == "paired" else 1
        if len(self.urls) != expected_urls:
            raise ValueError(f"{self.layout} layout requires exactly {expected_urls} URL(s)")
        if self.checksums and len(self.checksums) != len(self.urls):
            raise ValueError("checksums must be empty or have one value per URL")
        if self.expected_bytes and len(self.expected_bytes) != len(self.urls):
            raise ValueError("expected_bytes must be empty or have one value per URL")
        if any(size <= 0 for size in self.expected_bytes):
            raise ValueError("expected_bytes values must be positive")
        return self


class AcquireReadsConfig(StageConfig):
    archives: list[ReadArchive] = Field(default_factory=list)
    max_download_bytes: int = Field(
        default=DEFAULT_DOWNLOAD_CEILING_BYTES, ge=1, le=MAX_DOWNLOAD_CEILING_BYTES
    )
    confirmation_threshold_bytes: int = Field(
        default=DEFAULT_CONFIRMATION_THRESHOLD_BYTES, ge=0, le=MAX_DOWNLOAD_CEILING_BYTES
    )
    confirm_large_transfer: bool = False
    random_seed: int = 1234
    timeout_seconds: int = Field(default=120, ge=5, le=600)
    minimum_free_overhead_bytes: int = Field(default=268_435_456, ge=0)

    @model_validator(mode="after")
    def _validate_threshold(self) -> AcquireReadsConfig:
        if not self.archives:
            raise ValueError("acquire_reads needs at least one archive")
        ids = [archive.sample_id for archive in self.archives]
        if len(ids) != len(set(ids)):
            raise ValueError("acquire_reads archive sample_id values must be unique")
        if self.confirmation_threshold_bytes > self.max_download_bytes:
            raise ValueError("confirmation_threshold_bytes cannot exceed max_download_bytes")
        return self


@dataclass(frozen=True, slots=True)
class ArchiveEstimate:
    sample_id: str
    layout: str
    compressed_bytes: int


def estimate_archives(
    archives: list[ReadArchive],
    *,
    timeout_seconds: int = 30,
    content_length: dict[str, int] | None = None,
) -> list[ArchiveEstimate]:
    """Estimate compressed transfer size using supplied metadata or HEAD/range metadata."""
    estimates: list[ArchiveEstimate] = []
    for archive in archives:
        sizes = archive.expected_bytes or [
            (content_length or {}).get(url) or _remote_size(url, timeout_seconds)
            for url in archive.urls
        ]
        estimates.append(
            ArchiveEstimate(
                sample_id=archive.sample_id,
                layout=archive.layout,
                compressed_bytes=sum(sizes),
            )
        )
    return estimates


def enforce_read_transfer_policy(
    estimates: list[ArchiveEstimate],
    *,
    max_download_bytes: int,
    confirmation_threshold_bytes: int,
    confirmed: bool,
) -> int:
    """Return total estimated bytes or raise before any archive body is requested."""
    total = sum(estimate.compressed_bytes for estimate in estimates)
    if total > max_download_bytes:
        raise ReadAcquisitionError(
            f"estimated archive download {format_bytes(total)} exceeds hard ceiling "
            f"{format_bytes(max_download_bytes)}",
            recovery=(
                "Choose fewer samples, use a lower downsampled/public preset, or prepare a "
                "smaller archive set. The hard ceiling cannot be raised to 10 GB or more."
            ),
        )
    if total > confirmation_threshold_bytes and not confirmed:
        raise ReadAcquisitionError(
            f"estimated archive download is {format_bytes(total)}; "
            "explicit confirmation is required",
            recovery=(
                "Review the listed archive sizes, then set confirm_large_transfer: true only "
                "if the transfer and disk use are acceptable."
            ),
        )
    return total


@STAGES.register("acquire_reads")
class AcquireReadsStage(Stage):
    """Stream selected public FASTQ archives into a manifest consumable by ``align``."""

    name = "acquire_reads"
    config_model = AcquireReadsConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet", optional=True)]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.RAW_READS, "fastq_manifest")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: AcquireReadsConfig = self.config  # type: ignore[assignment]
        if ctx.datastore.has(ArtifactKind.SAMPLE_SHEET, "sample_sheet"):
            sheet = pd.read_csv(
                ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path, dtype=str
            )
            known = set(sheet.get("sample_id", pd.Series(dtype=str)).dropna().astype(str))
            requested = {archive.sample_id for archive in cfg.archives}
            missing = sorted(requested - known)
            if missing:
                raise ReadAcquisitionError(
                    "archive sample_id values are absent from the registered sample sheet: "
                    + ", ".join(missing),
                    recovery=(
                        "Correct the sample sheet or archive sample_id values before retrying."
                    ),
                )
        estimates = estimate_archives(cfg.archives, timeout_seconds=cfg.timeout_seconds)
        total = enforce_read_transfer_policy(
            estimates,
            max_download_bytes=cfg.max_download_bytes,
            confirmation_threshold_bytes=cfg.confirmation_threshold_bytes,
            confirmed=cfg.confirm_large_transfer,
        )
        stage_dir = ctx.datastore.stage_dir(self.name)
        _require_disk_space(stage_dir, total + cfg.minimum_free_overhead_bytes)

        rows: list[dict[str, object]] = []
        checksums: list[dict[str, object]] = []
        for index, archive in enumerate(cfg.archives):
            sample_rng = random.Random(cfg.random_seed + index)
            outputs, records = download_read_archive(
                archive,
                stage_dir,
                rng=sample_rng,
                timeout_seconds=cfg.timeout_seconds,
            )
            checksums.extend(records)
            row: dict[str, object] = {
                "sample_id": archive.sample_id,
                "layout": archive.layout,
                "fastq1": str(outputs[0]),
                "fastq2": str(outputs[1]) if archive.layout == "paired" else "",
                "downsample_fraction": archive.downsample_fraction,
                "source_urls": ";".join(archive.urls),
            }
            rows.append(row)

        manifest = pd.DataFrame(rows)
        manifest_path = stage_dir / "fastq_manifest.csv"
        _atomic_csv(manifest, manifest_path)
        provenance_path = stage_dir / "read_acquisition_manifest.json"
        _atomic_json(
            provenance_path,
            {
                "schema_version": 1,
                "created_at": datetime.now(UTC).isoformat(timespec="seconds"),
                "estimated_download_bytes": total,
                "hard_ceiling_bytes": cfg.max_download_bytes,
                "confirmation_threshold_bytes": cfg.confirmation_threshold_bytes,
                "confirmed_large_transfer": cfg.confirm_large_transfer,
                "archives": checksums,
                "cleanup": "temporary .part files are removed after success or failure",
            },
        )
        artifact = ctx.datastore.add(
            ArtifactKind.RAW_READS,
            "fastq_manifest",
            manifest_path,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "n_samples": len(rows),
                "estimated_download_bytes": total,
                "manifest": str(provenance_path),
                "layouts": sorted({archive.layout for archive in cfg.archives}),
            },
        )
        return StageResult(
            artifacts=[artifact],
            metrics={
                "n_samples": len(rows),
                "estimated_download_bytes": total,
                "layouts": sorted({archive.layout for archive in cfg.archives}),
            },
        )


def download_read_archive(
    archive: ReadArchive, output_dir: Path, *, rng: random.Random, timeout_seconds: int
) -> tuple[list[Path], list[dict[str, object]]]:
    """Stream one archive set, preserving paired/interleaved units during downsampling."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if archive.layout == "paired":
        paths = [
            output_dir / f"{archive.sample_id}_R1.fastq.gz",
            output_dir / f"{archive.sample_id}_R2.fastq.gz",
        ]
        records = _download_paired(archive, paths, rng=rng, timeout_seconds=timeout_seconds)
        return paths, records

    suffix = "interleaved" if archive.layout == "interleaved" else "single"
    path = output_dir / f"{archive.sample_id}_{suffix}.fastq.gz"
    records_per_unit = 2 if archive.layout == "interleaved" else 1
    record = _download_one(
        archive.urls[0],
        path,
        checksum=archive.checksums[0] if archive.checksums else "",
        fraction=archive.downsample_fraction,
        rng=rng,
        records_per_unit=records_per_unit,
        timeout_seconds=timeout_seconds,
    )
    return [path], [record]


def _download_paired(
    archive: ReadArchive, paths: list[Path], *, rng: random.Random, timeout_seconds: int
) -> list[dict[str, object]]:
    temporary = [path.with_suffix(path.suffix + ".part") for path in paths]
    responses: list[Any] = []
    try:
        responses = [_open_stream(url, timeout_seconds) for url in archive.urls]
        raw = [_HashingReader(response) for response in responses]
        readers = [gzip.GzipFile(fileobj=item, mode="rb") for item in raw]
        with gzip.open(temporary[0], "wb") as out1, gzip.open(temporary[1], "wb") as out2:
            while True:
                pair = [_fastq_record(reader) for reader in readers]
                if pair == [None, None]:
                    break
                if any(item is None for item in pair):
                    raise ReadAcquisitionError(
                        f"paired archive for {archive.sample_id} has unmatched FASTQ records",
                        recovery=(
                            "Verify the publisher paired-read files and retry with matching "
                            "R1/R2 URLs."
                        ),
                    )
                if rng.random() <= archive.downsample_fraction:
                    out1.writelines(pair[0] or [])
                    out2.writelines(pair[1] or [])
        for reader in readers:
            reader.close()
        records: list[dict[str, object]] = []
        for url, item, expected in zip(
            archive.urls, raw, archive.checksums or ["", ""], strict=True
        ):
            digests = item.digests()
            _verify_expected_digest(digests, expected, url)
            records.append({"url": url, "source_digests": digests, "expected_checksum": expected})
        for source, destination in zip(temporary, paths, strict=True):
            os.replace(source, destination)
        return records
    except Exception:
        for path in temporary:
            path.unlink(missing_ok=True)
        raise
    finally:
        for response in responses:
            close = getattr(response, "close", None)
            if close:
                close()


def _download_one(
    url: str,
    destination: Path,
    *,
    checksum: str,
    fraction: float,
    rng: random.Random,
    records_per_unit: int,
    timeout_seconds: int,
) -> dict[str, object]:
    temporary = destination.with_suffix(destination.suffix + ".part")
    response = None
    try:
        response = _open_stream(url, timeout_seconds)
        raw = _HashingReader(response)
        reader = gzip.GzipFile(fileobj=raw, mode="rb")
        with gzip.open(temporary, "wb") as output:
            while True:
                unit = [_fastq_record(reader) for _ in range(records_per_unit)]
                if all(record is None for record in unit):
                    break
                if any(record is None for record in unit):
                    raise ReadAcquisitionError(
                        f"interleaved/single archive {url} ended mid-read unit",
                        recovery=(
                            "Verify the archive is valid gzip FASTQ and retry the acquisition."
                        ),
                    )
                if rng.random() <= fraction:
                    for record in unit:
                        output.writelines(record or [])
        reader.close()
        digests = raw.digests()
        _verify_expected_digest(digests, checksum, url)
        os.replace(temporary, destination)
        return {"url": url, "source_digests": digests, "expected_checksum": checksum}
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        if response is not None:
            response.close()


def _fastq_record(reader: gzip.GzipFile | BinaryIO) -> list[bytes] | None:
    first = reader.readline()
    if not first:
        return None
    rest = [reader.readline() for _ in range(3)]
    if any(not line for line in rest):
        raise ReadAcquisitionError(
            "FASTQ ended partway through a record",
            recovery="Verify the source archive integrity and retry the acquisition.",
        )
    return [first, *rest]


class _HashingReader(io.RawIOBase):
    """A non-buffering wrapper that hashes compressed source bytes as gzip consumes them."""

    def __init__(self, source: BinaryIO) -> None:
        self.source = source
        self._sha256 = hashlib.sha256()
        self._md5 = hashlib.md5()

    def readable(self) -> bool:
        return True

    def read(self, size: int = -1) -> bytes:
        data = self.source.read(size)
        if data:
            self._sha256.update(data)
            self._md5.update(data)
        return data

    def readinto(self, buffer: bytearray) -> int:  # type: ignore[override]
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)

    def digests(self) -> dict[str, str]:
        return {"sha256": self._sha256.hexdigest(), "md5": self._md5.hexdigest()}


def _verify_expected_digest(digests: dict[str, str], expected: str, url: str) -> None:
    if not expected:
        return
    algorithm, separator, value = expected.lower().partition(":")
    if not separator:
        if len(expected) == 32:
            algorithm, value = "md5", expected.lower()
        elif len(expected) == 64:
            algorithm, value = "sha256", expected.lower()
        else:
            raise ReadAcquisitionError(
                f"unrecognized checksum format for {url}: {expected!r}",
                recovery="Use md5:<digest> or sha256:<digest> from the archive publisher.",
            )
    if algorithm not in digests or digests[algorithm] != value:
        raise ReadAcquisitionError(
            f"checksum mismatch for {url}",
            recovery="Do not use the downloaded archive. Check the publisher checksum and retry.",
        )


def _remote_size(url: str, timeout_seconds: int) -> int:
    try:
        response = urlopen(
            Request(url, method="HEAD", headers={"User-Agent": "CANIS/0.2"}),
            timeout=timeout_seconds,
        )
        with response:
            length = response.headers.get("Content-Length")
            if length and int(length) > 0:
                return int(length)
    except Exception:
        pass
    try:
        request = Request(url, headers={"Range": "bytes=0-0", "User-Agent": "CANIS/0.2"})
        with urlopen(request, timeout=timeout_seconds) as response:
            content_range = response.headers.get("Content-Range", "")
            if "/" in content_range:
                return int(content_range.rsplit("/", 1)[1])
    except Exception as exc:
        raise ReadAcquisitionError(
            f"could not estimate archive size for {url}: {exc}",
            recovery=(
                "Provide expected_bytes in the archive config or use an HTTP source exposing "
                "Content-Length/Content-Range before retrying."
            ),
        ) from exc
    raise ReadAcquisitionError(
        f"archive size is not available for {url}",
        recovery="Provide expected_bytes explicitly; CANIS will not download an unbounded archive.",
    )


def _open_stream(url: str, timeout_seconds: int):
    return urlopen(Request(url, headers={"User-Agent": "CANIS/0.2"}), timeout=timeout_seconds)


def _require_disk_space(path: Path, required_bytes: int) -> None:
    free = shutil.disk_usage(path).free
    if free < required_bytes:
        raise ReadAcquisitionError(
            f"only {format_bytes(free)} free at {path}; acquisition needs at least "
            f"{format_bytes(required_bytes)}",
            recovery=(
                "Free disk space or reduce the archive selection/downsample fraction before "
                "retrying."
            ),
        )


def _atomic_csv(table: pd.DataFrame, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    try:
        table.to_csv(temporary, index=False)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(destination: Path, payload: dict) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


__all__ = [
    "AcquireReadsConfig",
    "AcquireReadsStage",
    "ArchiveEstimate",
    "ReadAcquisitionError",
    "ReadArchive",
    "download_read_archive",
    "enforce_read_transfer_policy",
    "estimate_archives",
]
