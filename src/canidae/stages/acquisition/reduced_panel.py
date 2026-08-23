"""Laptop-safe extraction of a compact genotype panel from a remote indexed VCF.

``ReducedPanelStage`` is an acquisition entry point for published, bgzip-compressed VCFs
with a tabix ``.tbi`` index.  It deliberately downloads only byte ranges that overlap a
fixed marker panel, selects the requested samples, validates biallelic complete GT calls,
and writes a compact local VCF plus a machine-readable extraction manifest.

The stage is intentionally dependency-light: it uses the Python standard library rather
than requiring ``bcftools``.  It is suitable for the 2K, 10K, and 25K laptop presets, not
for silently acquiring full genomes.  A transfer estimate is calculated before source VCF
ranges are requested; the default hard ceiling is below 10 GB and estimates above the
confirmation threshold require an explicit opt-in in the configuration.
"""

from __future__ import annotations

import gzip
import hashlib
import http.client
import json
import os
import shutil
import struct
import tempfile
import threading
import time
import urllib.error
import urllib.request
from collections import Counter
from collections.abc import Callable, Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import pandas as pd
from pydantic import Field, field_validator, model_validator

from canidae.core.atomic import atomic_write_text
from canidae.core.errors import IntegrityError, StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.acquisition.sample_sheet import read_sample_sheet

# Named presets are deliberately small enough for a laptop analysis.  The values are the
# target number of fixed marker-panel positions before call validation removes unavailable
# or incomplete variants.
PANEL_PRESETS: dict[str, int] = {"2k": 2_000, "10k": 10_000, "25k": 25_000}
DEFAULT_DOWNLOAD_CEILING_BYTES = 9_000_000_000
MAX_DOWNLOAD_CEILING_BYTES = 9_999_999_999
DEFAULT_CONFIRMATION_THRESHOLD_BYTES = 1_000_000_000

DEFAULT_SOURCE_VCF_URL = (
    "https://research.nhgri.nih.gov/dog_genome/downloads/datasets/WGS/"
    "722g.990.SNP.INDEL.chrAll.vcf.gz"
)
DEFAULT_PANEL_URL = "https://download.cncb.ac.cn/dogsd/dog10k/variations/CFA31_IlluminaHD.vcf.gz"
DEFAULT_PANEL_MD5 = "4a87088b17631fb6237210044ac099a6"

_BGZF_MAX_BLOCK = 65_536
_HEADER_ALLOWANCE_BYTES = 16 * 1024 * 1024
_MAX_ATTEMPTS = 5
_BACKOFF_SECONDS = 2.0
_RETRYABLE_HTTP_STATUS = frozenset({429, 500, 502, 503, 504})


class ReducedPanelError(StageInputError):
    """An expected, actionable reduced-panel acquisition failure."""


class TransferConfirmationRequiredError(ReducedPanelError):
    """The estimate is safe but needs explicit user confirmation before transfer."""


class TransferLimitExceededError(ReducedPanelError):
    """The estimate or an actual response exceeds the configured safety ceiling."""


@dataclass(frozen=True, slots=True)
class ReferenceIndex:
    """The bins and linear index for one tabix reference sequence."""

    bins: dict[int, list[tuple[int, int]]]
    linear: list[int]


@dataclass(slots=True)
class TransferBudget:
    """Tracks actual bytes and refuses transfers over the configured ceiling."""

    limit: int
    used: int = 0
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    def reserve(self, amount: int) -> None:
        with self._lock:
            self._reserve_unlocked(amount)

    def add(self, amount: int) -> None:
        with self._lock:
            self._reserve_unlocked(amount)
            self.used += amount

    def _reserve_unlocked(self, amount: int) -> None:
        if amount < 0 or self.used + amount > self.limit:
            raise TransferLimitExceededError(
                "transfer would exceed laptop safety ceiling: "
                f"{self.used + amount:,} > {self.limit:,} bytes"
            )


@dataclass(slots=True)
class RangeCache:
    """Persistent validated cache for remote BGZF byte ranges."""

    root: Path
    namespace: str
    hits: int = 0
    reused_bytes: int = 0

    @property
    def directory(self) -> Path:
        return self.root / self.namespace

    def get(self, start: int, end: int) -> bytes | None:
        path = self._path(start, end)
        digest_path = path.with_suffix(".sha256")
        try:
            data = path.read_bytes()
            expected = digest_path.read_text(encoding="ascii").strip()
        except OSError:
            return None
        if len(data) != end - start + 1 or hashlib.sha256(data).hexdigest() != expected:
            path.unlink(missing_ok=True)
            digest_path.unlink(missing_ok=True)
            return None
        self.hits += 1
        self.reused_bytes += len(data)
        return data

    def put(self, start: int, end: int, data: bytes) -> None:
        if len(data) != end - start + 1:
            return
        path = self._path(start, end)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(f".tmp-{threading.get_ident()}")
        digest_tmp = path.with_suffix(f".sha256.tmp-{threading.get_ident()}")
        tmp.write_bytes(data)
        digest_tmp.write_text(hashlib.sha256(data).hexdigest(), encoding="ascii")
        os.replace(tmp, path)
        os.replace(digest_tmp, path.with_suffix(".sha256"))

    def _path(self, start: int, end: int) -> Path:
        return self.directory / f"{start}-{end}.bin"


@dataclass(frozen=True, slots=True)
class TransferEstimate:
    """A conservative estimate made before source VCF ranges are fetched."""

    prerequisite_bytes: int
    range_bytes: int
    header_allowance_bytes: int
    n_ranges: int
    n_target_sites: int

    @property
    def total_bytes(self) -> int:
        return self.prerequisite_bytes + self.range_bytes + self.header_allowance_bytes

    def as_dict(self) -> dict[str, int]:
        return {**asdict(self), "total_bytes": self.total_bytes}


@dataclass(frozen=True, slots=True)
class ReducedPanelResult:
    """The paths and accounting emitted by one successful extraction."""

    vcf_path: Path
    manifest_path: Path
    retained_sites: int
    transfer_estimate: TransferEstimate
    downloaded_bytes: int
    skipped_records: dict[str, int]
    output_sha256: str
    cached_ranges: int = 0
    reused_bytes: int = 0


def preset_site_count(preset: str) -> int:
    """Return the fixed marker target for a named laptop preset."""
    try:
        return PANEL_PRESETS[preset]
    except KeyError as exc:
        choices = ", ".join(PANEL_PRESETS)
        raise ReducedPanelError(
            f"unknown reduced-panel preset '{preset}'; choose {choices}"
        ) from exc


def format_bytes(value: float) -> str:
    """Render a byte count for an error message or user-facing progress output."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1000 or unit == "TB":
            return f"{value:.2f} {unit}" if unit != "B" else f"{value} B"
        value /= 1000
    return f"{value:.2f} TB"  # pragma: no cover - loop always returns


def enforce_transfer_policy(
    estimate: TransferEstimate,
    *,
    ceiling_bytes: int,
    confirmation_threshold_bytes: int,
    confirmed: bool,
) -> None:
    """Reject unsafe work, or require an explicit confirmation for a large safe transfer."""
    if estimate.total_bytes > ceiling_bytes:
        raise TransferLimitExceededError(
            "indexed extraction estimate exceeds the configured safety ceiling: "
            f"{format_bytes(estimate.total_bytes)} > {format_bytes(ceiling_bytes)}"
        )
    if estimate.total_bytes > confirmation_threshold_bytes and not confirmed:
        raise TransferConfirmationRequiredError(
            "indexed extraction is estimated at "
            f"{format_bytes(estimate.total_bytes)} across {estimate.n_ranges:,} byte ranges. "
            "Set stages.reduced_panel.confirm_large_transfer=true after reviewing the estimate."
        )


def verify_checksum(data: bytes, expected: str, *, label: str) -> str:
    """Verify a ``md5`` or ``sha256`` checksum and return the normalized digest string.

    A bare 32-character value is treated as MD5 for compatibility with common public data
    catalogues.  An empty expected value records the SHA-256 without asserting a publisher
    checksum, which is appropriate for sources that do not provide one.
    """
    sha256 = hashlib.sha256(data).hexdigest()
    if not expected:
        return f"sha256:{sha256}"
    algorithm, _, wanted = expected.partition(":")
    if not wanted:
        algorithm, wanted = ("md5", algorithm) if len(algorithm) == 32 else ("sha256", algorithm)
    algorithm = algorithm.lower()
    if algorithm == "md5":
        actual = hashlib.md5(data).hexdigest()
    elif algorithm == "sha256":
        actual = sha256
    else:
        raise ReducedPanelError(f"unsupported checksum algorithm for {label}: {algorithm}")
    if actual.lower() != wanted.lower():
        raise IntegrityError(
            f"{label} checksum mismatch: expected {algorithm}:{wanted}, got {algorithm}:{actual}"
        )
    return f"{algorithm}:{actual}"


def choose_samples(sample_sheet: pd.DataFrame, requested: list[str]) -> pd.DataFrame:
    """Select a validated ordered subset from a normalized sample sheet."""
    if not requested:
        return sample_sheet.copy().reset_index(drop=True)
    if len(requested) != len(set(requested)):
        raise ReducedPanelError("selected_samples contains duplicate sample IDs")
    by_id = sample_sheet.set_index("sample_id", drop=False)
    missing = [sample_id for sample_id in requested if sample_id not in by_id.index]
    if missing:
        raise ReducedPanelError(
            "selected_samples must be present in sample_sheet; missing: " + ", ".join(missing)
        )
    return by_id.loc[requested].reset_index(drop=True)


def _request(
    url: str,
    budget: TransferBudget,
    byte_range: tuple[int, int] | None = None,
    *,
    timeout_seconds: int,
    progress: Callable[[str], None] | None = None,
) -> bytes:
    """Fetch bytes with bounded retries and strict HTTP range-response validation."""
    headers = {"User-Agent": "CANIS-reduced-panel/0.1"}
    if byte_range is not None:
        start, end = byte_range
        headers["Range"] = f"bytes={start}-{end}"
    request = urllib.request.Request(url, headers=headers)
    last_error: Exception | None = None
    for attempt in range(1, _MAX_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                status = getattr(response, "status", response.getcode())
                if byte_range is not None and status != 206:
                    raise TransferLimitExceededError(
                        f"server ignored byte-range request for {url} (HTTP {status}); "
                        "aborting to prevent a full source download"
                    )
                declared = response.headers.get("Content-Length")
                if declared is not None:
                    budget.reserve(int(declared))
                data = response.read()
            budget.add(len(data))
            return data
        except TransferLimitExceededError:
            raise
        except urllib.error.HTTPError as exc:
            if exc.code not in _RETRYABLE_HTTP_STATUS:
                raise ReducedPanelError(f"download failed for {url}: HTTP {exc.code}") from exc
            last_error = exc
        except (urllib.error.URLError, http.client.IncompleteRead, OSError) as exc:
            last_error = exc
        if attempt < _MAX_ATTEMPTS:
            delay = _BACKOFF_SECONDS * (2 ** (attempt - 1))
            if progress:
                progress(
                    f"transfer attempt {attempt}/{_MAX_ATTEMPTS} failed ({last_error}); "
                    f"retrying in {delay:.0f}s"
                )
            time.sleep(delay)
    raise ReducedPanelError(
        f"transfer failed after {_MAX_ATTEMPTS} attempts for {url}: {last_error}"
    )


def _parse_tabix(raw_bgzf: bytes) -> tuple[list[str], list[ReferenceIndex]]:
    """Parse the small gzipped tabix index required for indexed byte-range selection."""
    raw = gzip.decompress(raw_bgzf)
    if raw[:4] != b"TBI\x01":
        raise ReducedPanelError("source index is not a tabix TBI index")
    offset = 4

    def i32() -> int:
        nonlocal offset
        value = struct.unpack_from("<i", raw, offset)[0]
        offset += 4
        return value

    n_ref = i32()
    for _ in range(6):
        i32()  # format, sequence/begin/end columns, meta character, skipped lines
    names_len = i32()
    names = raw[offset : offset + names_len].rstrip(b"\x00").decode().split("\x00")
    offset += names_len
    refs: list[ReferenceIndex] = []
    for _ in range(n_ref):
        n_bin = i32()
        bins: dict[int, list[tuple[int, int]]] = {}
        for _ in range(n_bin):
            bin_id = struct.unpack_from("<I", raw, offset)[0]
            offset += 4
            n_chunk = i32()
            chunks = []
            for _ in range(n_chunk):
                chunks.append(struct.unpack_from("<QQ", raw, offset))
                offset += 16
            bins[bin_id] = chunks
        n_linear = i32()
        linear = list(struct.unpack_from(f"<{n_linear}Q", raw, offset)) if n_linear else []
        offset += 8 * n_linear
        refs.append(ReferenceIndex(bins=bins, linear=linear))
    if len(names) != len(refs):
        raise ReducedPanelError(
            f"tabix reference mismatch: {len(names)} names, {len(refs)} indexes"
        )
    return names, refs


def _reg2bins(start: int, end: int) -> list[int]:
    end -= 1
    bins = [0]
    bins.extend(range(1 + (start >> 26), 1 + (end >> 26) + 1))
    bins.extend(range(9 + (start >> 23), 9 + (end >> 23) + 1))
    bins.extend(range(73 + (start >> 20), 73 + (end >> 20) + 1))
    bins.extend(range(585 + (start >> 17), 585 + (end >> 17) + 1))
    bins.extend(range(4681 + (start >> 14), 4681 + (end >> 14) + 1))
    return bins


def _chunks_for(index: ReferenceIndex, pos_1based: int) -> list[tuple[int, int]]:
    start = pos_1based - 1
    linear_i = start >> 14
    min_offset = index.linear[linear_i] if linear_i < len(index.linear) else 0
    chunks: list[tuple[int, int]] = []
    for bin_id in _reg2bins(start, start + 1):
        chunks.extend(chunk for chunk in index.bins.get(bin_id, []) if chunk[1] > min_offset)
    return chunks


def _merge_chunks(chunks: set[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(chunks)
    merged: list[list[int]] = []
    for begin, end in ordered:
        if not merged:
            merged.append([begin, end])
            continue
        previous = merged[-1]
        previous_byte_end = (previous[1] >> 16) + _BGZF_MAX_BLOCK
        current_byte_start = begin >> 16
        if current_byte_start <= previous_byte_end:
            previous[1] = max(previous[1], end)
        else:
            merged.append([begin, end])
    return [(begin, end) for begin, end in merged]


def _decompress_bgzf(data: bytes, first_uncompressed_offset: int = 0) -> bytes:
    """Decompress complete BGZF blocks returned by an HTTP byte-range request."""
    output = bytearray()
    cursor = 0
    while cursor + 18 <= len(data):
        if data[cursor : cursor + 2] != b"\x1f\x8b":
            break
        block_size = struct.unpack_from("<H", data, cursor + 16)[0] + 1
        if cursor + block_size > len(data):
            break
        output.extend(gzip.decompress(data[cursor : cursor + block_size]))
        cursor += block_size
    return bytes(output[first_uncompressed_offset:])


def _panel_positions(panel_data: bytes, target_sites: int) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    for line in gzip.decompress(panel_data).decode().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t", 3)
        if fields[0].isdigit() and 1 <= int(fields[0]) <= 38:
            rows.append((fields[0], int(fields[1])))
    if not rows:
        raise ReducedPanelError("no autosomal positions found in marker panel")
    stride = max(1, (len(rows) + target_sites - 1) // target_sites)
    return rows[::stride]


def _resolve_contig(panel_contig: str, source_contigs: set[str]) -> str | None:
    """Map common numeric/``chr`` autosomal naming styles without invoking liftover."""
    candidates = (panel_contig, f"chr{panel_contig}")
    for candidate in candidates:
        if candidate in source_contigs:
            return candidate
    return None


def _source_header(
    url: str,
    budget: TransferBudget,
    *,
    timeout_seconds: int,
    progress: Callable[[str], None] | None,
) -> tuple[list[str], list[str]]:
    fetch_size = 1_048_576
    while fetch_size <= _HEADER_ALLOWANCE_BYTES:
        compressed = _request(
            url,
            budget,
            (0, fetch_size - 1),
            timeout_seconds=timeout_seconds,
            progress=progress,
        )
        text = _decompress_bgzf(compressed).decode(errors="replace")
        lines = text.splitlines()
        column_line = next((line for line in lines if line.startswith("#CHROM\t")), None)
        if column_line is not None:
            return [line for line in lines if line.startswith("##")], column_line.split("\t")
        fetch_size *= 2
    raise ReducedPanelError("VCF header exceeded 16 MB compressed; refusing to continue")


def _genotype(
    sample_field: str,
    formats: list[str],
    *,
    min_gq: float,
    min_dp: int,
) -> str | None:
    values = sample_field.split(":")
    gt_index = formats.index("GT")
    if gt_index >= len(values):
        return None
    gt = values[gt_index]
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    for name, threshold in (("GQ", min_gq), ("DP", min_dp)):
        if threshold <= 0:
            continue
        if name not in formats:
            return None
        index = formats.index(name)
        try:
            value = float(values[index])
        except (IndexError, ValueError):
            return None
        if value < threshold:
            return None
    return gt


def _contig_sort_key(contig: str) -> tuple[int, int | str]:
    raw = contig.removeprefix("chr")
    return (0, int(raw)) if raw.isdigit() else (1, contig)


def _atomic_gzip_vcf(
    output_path: Path,
    meta_lines: list[str],
    columns: list[str],
    sample_ids: list[str],
    records: dict[tuple[str, int], str],
    *,
    preset: str,
) -> None:
    """Write the compact VCF to a neighboring temporary file before promotion."""
    with tempfile.NamedTemporaryFile(
        mode="wb",
        suffix=".vcf.gz.tmp",
        prefix=f".{output_path.stem}.",
        dir=output_path.parent,
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        with gzip.open(temporary, "wt", encoding="utf-8", compresslevel=6) as out:
            for line in meta_lines:
                out.write(f"{line}\n")
            out.write(f"##CANIS_reduced_panel=preset:{preset}\n")
            out.write("\t".join(columns[:9] + sample_ids) + "\n")
            for key in sorted(records, key=lambda item: (_contig_sort_key(item[0]), item[1])):
                out.write(f"{records[key]}\n")
        temporary.replace(output_path)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _indexed_plan(
    panel_data: bytes,
    index_data: bytes,
    target_sites: int,
    prerequisite_bytes: int,
) -> tuple[set[tuple[str, int]], list[tuple[int, int]], TransferEstimate]:
    names, indexes = _parse_tabix(index_data)
    name_to_index = {name: index for name, index in zip(names, indexes, strict=True)}
    positions = _panel_positions(panel_data, target_sites)
    source_contigs = set(name_to_index)
    targets = {
        (mapped, pos)
        for chrom, pos in positions
        if (mapped := _resolve_contig(chrom, source_contigs)) is not None
    }
    if not targets:
        raise ReducedPanelError("marker-panel chromosome names do not match the source index")
    chunks: set[tuple[int, int]] = set()
    for chrom, pos in targets:
        chunks.update(_chunks_for(name_to_index[chrom], pos))
    merged = _merge_chunks(chunks)
    range_bytes = sum(((end >> 16) + _BGZF_MAX_BLOCK) - (begin >> 16) for begin, end in merged)
    estimate = TransferEstimate(
        prerequisite_bytes=prerequisite_bytes,
        range_bytes=range_bytes,
        header_allowance_bytes=_HEADER_ALLOWANCE_BYTES,
        n_ranges=len(merged),
        n_target_sites=len(targets),
    )
    return targets, merged, estimate


def preflight_indexed_panel(
    *,
    target_sites: int,
    source_vcf_url: str = DEFAULT_SOURCE_VCF_URL,
    panel_url: str = DEFAULT_PANEL_URL,
    panel_checksum: str = DEFAULT_PANEL_MD5,
    index_checksum: str = "",
    max_download_bytes: int = DEFAULT_DOWNLOAD_CEILING_BYTES,
    timeout_seconds: int = 120,
) -> TransferEstimate:
    """Calculate the exact indexed transfer plan without fetching source VCF ranges."""
    budget = TransferBudget(max_download_bytes)
    panel_data = _request(panel_url, budget, timeout_seconds=timeout_seconds)
    verify_checksum(panel_data, panel_checksum, label="marker panel")
    index_data = _request(f"{source_vcf_url}.tbi", budget, timeout_seconds=timeout_seconds)
    verify_checksum(index_data, index_checksum, label="tabix index")
    return _indexed_plan(panel_data, index_data, target_sites, budget.used)[2]


def extract_indexed_panel(
    *,
    source_vcf_url: str,
    panel_url: str,
    panel_checksum: str,
    index_checksum: str,
    sample_ids: list[str],
    target_sites: int,
    preset: str,
    output_path: Path,
    manifest_path: Path,
    max_download_bytes: int = DEFAULT_DOWNLOAD_CEILING_BYTES,
    confirmation_threshold_bytes: int = DEFAULT_CONFIRMATION_THRESHOLD_BYTES,
    confirm_large_transfer: bool = False,
    min_retained_sites: int | None = None,
    timeout_seconds: int = 120,
    temporary_dir: Path | None = None,
    range_cache_dir: Path | None = None,
    range_workers: int = 2,
    min_genotype_quality: float = 0.0,
    min_genotype_depth: int = 0,
    progress: Callable[[str], None] | None = None,
) -> ReducedPanelResult:
    """Extract a compact VCF from remote tabix-indexed data.

    The function is public to make the network-specific portion usable from a future UI and
    independently testable.  It never downloads the full source VCF: range responses must
    be HTTP 206 and every actual response is accounted against ``max_download_bytes``.
    """
    if not sample_ids:
        raise ReducedPanelError("at least one selected sample is required")
    if len(sample_ids) != len(set(sample_ids)):
        raise ReducedPanelError("selected sample IDs must be unique")
    if target_sites <= 0:
        raise ReducedPanelError("target_sites must be positive")
    if not 1 <= max_download_bytes <= MAX_DOWNLOAD_CEILING_BYTES:
        raise ReducedPanelError(
            "max_download_bytes must be positive and remain below 10 GB "
            f"(<= {MAX_DOWNLOAD_CEILING_BYTES:,})"
        )
    if confirmation_threshold_bytes < 0:
        raise ReducedPanelError("confirmation_threshold_bytes must be >= 0")
    if not 1 <= range_workers <= 4:
        raise ReducedPanelError("range_workers must be between 1 and 4")

    output_path = Path(output_path)
    manifest_path = Path(manifest_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_root = Path(temporary_dir) if temporary_dir else output_path.parent
    temp_root.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".reduced-panel-", dir=temp_root))
    budget = TransferBudget(max_download_bytes)
    try:
        if progress:
            progress("downloading marker panel and tabix index")
        panel_data = _request(panel_url, budget, timeout_seconds=timeout_seconds, progress=progress)
        panel_digest = verify_checksum(panel_data, panel_checksum, label="marker panel")
        (work / "marker_panel.vcf.gz").write_bytes(panel_data)

        index_url = f"{source_vcf_url}.tbi"
        index_data = _request(index_url, budget, timeout_seconds=timeout_seconds, progress=progress)
        index_digest = verify_checksum(index_data, index_checksum, label="tabix index")
        (work / "source.vcf.gz.tbi").write_bytes(index_data)

        targets, merged, estimate = _indexed_plan(panel_data, index_data, target_sites, budget.used)
        if progress:
            progress(
                f"estimated transfer {format_bytes(estimate.total_bytes)} for "
                f"{len(targets):,} targets in {len(merged):,} byte ranges"
            )
        enforce_transfer_policy(
            estimate,
            ceiling_bytes=max_download_bytes,
            confirmation_threshold_bytes=confirmation_threshold_bytes,
            confirmed=confirm_large_transfer,
        )

        meta_lines, columns = _source_header(
            source_vcf_url, budget, timeout_seconds=timeout_seconds, progress=progress
        )
        missing_samples = [sample for sample in sample_ids if sample not in columns]
        if missing_samples:
            raise ReducedPanelError(
                "selected samples are absent from the remote VCF header: "
                + ", ".join(missing_samples)
            )
        sample_columns = [columns.index(sample) for sample in sample_ids]

        records: dict[tuple[str, int], str] = {}
        skipped: Counter[str] = Counter()
        cache = None
        if range_cache_dir is not None:
            namespace = hashlib.sha256(f"{source_vcf_url}|{index_digest}".encode()).hexdigest()[:24]
            cache = RangeCache(Path(range_cache_dir), namespace)

        def fetch(item: tuple[int, tuple[int, int]]) -> tuple[int, int, bytes]:
            number, (virtual_begin, virtual_end) = item
            byte_start = virtual_begin >> 16
            byte_end = (virtual_end >> 16) + _BGZF_MAX_BLOCK - 1
            data = cache.get(byte_start, byte_end) if cache else None
            if data is None:
                data = _request(
                    source_vcf_url,
                    budget,
                    (byte_start, byte_end),
                    timeout_seconds=timeout_seconds,
                    progress=progress,
                )
                if cache:
                    cache.put(byte_start, byte_end, data)
            return number, virtual_begin, data

        range_started = time.monotonic()

        def bounded_fetches() -> Iterator[tuple[int, int, bytes]]:
            # Executor.map submits its entire iterable eagerly on supported Python
            # versions. Feed it small batches so completed BGZF ranges cannot build
            # up into a multi-gigabyte in-memory queue behind one slow response.
            indexed_ranges = list(enumerate(merged, start=1))
            batch_size = max(range_workers * 2, 1)
            with ThreadPoolExecutor(max_workers=range_workers) as pool:
                for offset in range(0, len(indexed_ranges), batch_size):
                    yield from pool.map(fetch, indexed_ranges[offset : offset + batch_size])

        for number, virtual_begin, data in bounded_fetches():
            text = _decompress_bgzf(data, virtual_begin & 0xFFFF).decode(errors="replace")
            for line in text.splitlines():
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) <= max(sample_columns):
                    skipped["truncated_record"] += 1
                    continue
                try:
                    key = (fields[0], int(fields[1]))
                except (ValueError, IndexError):
                    skipped["malformed_position"] += 1
                    continue
                if key not in targets or key in records:
                    continue
                if len(fields) < 10:
                    skipped["missing_format_or_samples"] += 1
                    continue
                ref, alt = fields[3], fields[4]
                if len(ref) != 1 or len(alt) != 1 or "," in alt:
                    skipped["not_biallelic_snp"] += 1
                    continue
                formats = fields[8].split(":")
                if "GT" not in formats:
                    skipped["missing_gt_field"] += 1
                    continue
                genotypes = [
                    _genotype(
                        fields[index],
                        formats,
                        min_gq=min_genotype_quality,
                        min_dp=min_genotype_depth,
                    )
                    for index in sample_columns
                ]
                if any(gt is None for gt in genotypes):
                    skipped["missing_or_non_diploid_gt"] += 1
                    continue
                records[key] = "\t".join(fields[:8] + ["GT"] + [str(gt) for gt in genotypes])
            if progress and (number % 100 == 0 or number == len(merged)):
                reused = cache.reused_bytes if cache else 0
                elapsed = max(time.monotonic() - range_started, 0.001)
                remaining_s = (len(merged) - number) / max(number / elapsed, 1e-9)
                progress(
                    f"ranges {number:,}/{len(merged):,}; downloaded "
                    f"{format_bytes(budget.used)}; reused {format_bytes(reused)}; "
                    f"retained {len(records):,} SNPs; ETA {remaining_s / 60:.1f} min"
                )

        minimum = min_retained_sites if min_retained_sites is not None else min(1_000, target_sites)
        if len(records) < minimum:
            raise ReducedPanelError(
                f"only {len(records):,} complete biallelic SNPs were retained; expected at least "
                f"{minimum:,}. Check the marker panel, selected samples, and source assembly."
            )
        _atomic_gzip_vcf(output_path, meta_lines, columns, sample_ids, records, preset=preset)
        output_sha256 = hashlib.sha256(output_path.read_bytes()).hexdigest()
        manifest = {
            "schema_version": 2,
            "kind": "canis_reduced_panel",
            "created_at": datetime.now(UTC).isoformat(timespec="seconds"),
            "source": {"vcf_url": source_vcf_url, "tabix_index_url": index_url},
            "selection": {"sample_ids": sample_ids, "target_sites": target_sites, "preset": preset},
            "marker_panel": {
                "url": panel_url,
                "checksum": panel_digest,
                "n_targeted_sites": len(targets),
            },
            "tabix_index": {"checksum": index_digest},
            "transfer": {
                "estimate": estimate.as_dict(),
                "downloaded_bytes": budget.used,
                "cached_ranges": cache.hits if cache else 0,
                "reused_bytes": cache.reused_bytes if cache else 0,
                "hard_ceiling_bytes": max_download_bytes,
                "confirmation_threshold_bytes": confirmation_threshold_bytes,
                "confirmed_large_transfer": confirm_large_transfer,
            },
            "validation": {
                "retained_sites": len(records),
                "minimum_retained_sites": minimum,
                "skipped_records": dict(sorted(skipped.items())),
                "required_calls": "complete diploid biallelic SNP GT for every selected sample",
                "minimum_genotype_quality": min_genotype_quality,
                "minimum_genotype_depth": min_genotype_depth,
            },
            "cleanup": {"temporary_marker_panel_and_index": "removed after extraction"},
            "output": {"path": str(output_path), "sha256": output_sha256},
        }
        _atomic_json(manifest_path, manifest)
        return ReducedPanelResult(
            vcf_path=output_path,
            manifest_path=manifest_path,
            retained_sites=len(records),
            transfer_estimate=estimate,
            downloaded_bytes=budget.used,
            skipped_records=dict(sorted(skipped.items())),
            output_sha256=output_sha256,
            cached_ranges=cache.hits if cache else 0,
            reused_bytes=cache.reused_bytes if cache else 0,
        )
    finally:
        shutil.rmtree(work, ignore_errors=True)


class ReducedPanelConfig(StageConfig):
    """Configuration for a reduced-panel acquisition run.

    ``sample_sheet`` is the user-facing sample selection.  Set ``selected_samples`` to a
    subset when a project sheet has more candidates than the laptop run should include.
    """

    sample_sheet: Path = Field(..., description="CSV containing samples to extract.")
    selected_samples: list[str] = Field(default_factory=list)
    source_vcf_url: str = DEFAULT_SOURCE_VCF_URL
    panel_url: str = DEFAULT_PANEL_URL
    panel_checksum: str = DEFAULT_PANEL_MD5
    index_checksum: str = ""
    preset: Literal["2k", "10k", "25k"] = "25k"
    max_download_bytes: int = Field(
        default=DEFAULT_DOWNLOAD_CEILING_BYTES, ge=1, le=MAX_DOWNLOAD_CEILING_BYTES
    )
    confirmation_threshold_bytes: int = Field(
        default=DEFAULT_CONFIRMATION_THRESHOLD_BYTES, ge=0, le=MAX_DOWNLOAD_CEILING_BYTES
    )
    confirm_large_transfer: bool = False
    min_retained_sites: int | None = Field(default=None, ge=1)
    timeout_seconds: int = Field(default=120, ge=5, le=600)
    range_workers: int = Field(default=2, ge=1, le=4)
    min_genotype_quality: float = Field(default=0.0, ge=0.0)
    min_genotype_depth: int = Field(default=0, ge=0)
    dataset_id: str = "remote_reduced_panel"
    reference_id: str = "CanFam3.1"

    @field_validator("selected_samples")
    @classmethod
    def _no_empty_sample_ids(cls, value: list[str]) -> list[str]:
        cleaned = [sample.strip() for sample in value]
        if any(not sample for sample in cleaned):
            raise ValueError("selected_samples cannot contain empty IDs")
        if len(cleaned) != len(set(cleaned)):
            raise ValueError("selected_samples cannot contain duplicate IDs")
        return cleaned

    @model_validator(mode="after")
    def _threshold_within_ceiling(self) -> ReducedPanelConfig:
        if self.confirmation_threshold_bytes > self.max_download_bytes:
            raise ValueError("confirmation_threshold_bytes cannot exceed max_download_bytes")
        return self


@STAGES.register("reduced_panel")
class ReducedPanelStage(Stage):
    """Acquire remote indexed VCF calls as normal pipeline artifacts.

    It produces the same ``callset`` and ``sample_sheet`` roles as ``ingest`` so a recipe can
    begin with ``reduced_panel -> qc -> load_genotypes -> ...`` without an out-of-band script.
    Do not include both ``ingest`` and ``reduced_panel`` in the same recipe because both are
    acquisition entry points that declare those roles.
    """

    name = "reduced_panel"
    config_model = ReducedPanelConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
            ArtifactSpec(ArtifactKind.CALLSET, "callset"),
        ]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: ReducedPanelConfig = self.config  # type: ignore[assignment]
        sheet_path = _resolve(cfg.sample_sheet, ctx.config.paths.root)
        selected = choose_samples(read_sample_sheet(sheet_path), cfg.selected_samples)
        sample_ids = selected["sample_id"].astype(str).tolist()
        target_sites = preset_site_count(cfg.preset)

        stage_dir = ctx.datastore.stage_dir(self.name)
        output_path = stage_dir / f"reduced_panel_{cfg.preset}.vcf.gz"
        manifest_path = stage_dir / f"reduced_panel_{cfg.preset}.manifest.json"

        def progress(message: str) -> None:
            # Stages should not print directly: a normal logger can later surface this in the UI.
            from canidae.core.logging import get_logger

            get_logger("stage.reduced_panel").info("%s", message)
            if ctx.run_dir is not None:
                atomic_write_text(
                    Path(ctx.run_dir) / "progress.json",
                    json.dumps(
                        {
                            "stage": self.name,
                            "message": message,
                            "updated_at": datetime.now(UTC).isoformat(timespec="seconds"),
                        },
                        indent=2,
                    )
                    + "\n",
                )

        result = extract_indexed_panel(
            source_vcf_url=cfg.source_vcf_url,
            panel_url=cfg.panel_url,
            panel_checksum=cfg.panel_checksum,
            index_checksum=cfg.index_checksum,
            sample_ids=sample_ids,
            target_sites=target_sites,
            preset=cfg.preset,
            output_path=output_path,
            manifest_path=manifest_path,
            max_download_bytes=cfg.max_download_bytes,
            confirmation_threshold_bytes=cfg.confirmation_threshold_bytes,
            confirm_large_transfer=cfg.confirm_large_transfer,
            min_retained_sites=cfg.min_retained_sites,
            timeout_seconds=cfg.timeout_seconds,
            temporary_dir=stage_dir,
            range_cache_dir=ctx.config.paths.cache_root / "reduced-panel-ranges",
            range_workers=cfg.range_workers,
            min_genotype_quality=cfg.min_genotype_quality,
            min_genotype_depth=cfg.min_genotype_depth,
            progress=progress,
        )
        normalized_sheet = stage_dir / "sample_sheet.csv"
        selected.to_csv(normalized_sheet, index=False)

        sheet_art = ctx.datastore.add(
            ArtifactKind.SAMPLE_SHEET,
            "sample_sheet",
            normalized_sheet,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"n_samples": len(selected), "source": str(sheet_path)},
        )
        callset_art = ctx.datastore.add(
            ArtifactKind.CALLSET,
            "callset",
            result.vcf_path,
            fmt=FileFormat.VCF,
            produced_by=self.name,
            metadata={
                "dataset_id": cfg.dataset_id,
                "source": cfg.source_vcf_url,
                "reference_id": cfg.reference_id,
                "reference_build": cfg.reference_id,
                "reduced_panel": True,
                "panel_relative": True,
                "panel_name": f"fixed_snp_panel_{cfg.preset}",
                "callable_sites": None,
                "preset": cfg.preset,
                "target_sites": target_sites,
                "retained_sites": result.retained_sites,
                "manifest": str(result.manifest_path),
                "sha256": result.output_sha256,
                "hard_call_filters": {
                    "min_genotype_quality": cfg.min_genotype_quality,
                    "min_genotype_depth": cfg.min_genotype_depth,
                },
            },
        )
        return StageResult(
            artifacts=[sheet_art, callset_art],
            metrics={
                "preset": cfg.preset,
                "target_sites": target_sites,
                "retained_sites": result.retained_sites,
                "n_samples": len(selected),
                "estimated_transfer_bytes": result.transfer_estimate.total_bytes,
                "downloaded_bytes": result.downloaded_bytes,
                "cached_ranges": result.cached_ranges,
                "reused_bytes": result.reused_bytes,
            },
        )


def _resolve(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else root / path


__all__ = [
    "DEFAULT_CONFIRMATION_THRESHOLD_BYTES",
    "DEFAULT_DOWNLOAD_CEILING_BYTES",
    "DEFAULT_PANEL_MD5",
    "DEFAULT_PANEL_URL",
    "DEFAULT_SOURCE_VCF_URL",
    "PANEL_PRESETS",
    "RangeCache",
    "ReducedPanelConfig",
    "ReducedPanelError",
    "ReducedPanelResult",
    "ReducedPanelStage",
    "TransferConfirmationRequiredError",
    "TransferEstimate",
    "TransferLimitExceededError",
    "choose_samples",
    "enforce_transfer_policy",
    "extract_indexed_panel",
    "format_bytes",
    "preflight_indexed_panel",
    "preset_site_count",
    "verify_checksum",
]
