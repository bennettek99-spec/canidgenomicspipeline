"""Indexed remote VCF byte-range helpers (tabix + BGZF).

Pure transport/index utilities used by laptop-safe panel extraction and the
hybrid-canid diagnostic scripts. Stage-level reduced-panel acquisition keeps its
own cache/progress layer on top of the same algorithms.
"""

from __future__ import annotations

import gzip
import http.client
import struct
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

SOURCE_VCF_URL = (
    "https://research.nhgri.nih.gov/dog_genome/downloads/datasets/WGS/"
    "722g.990.SNP.INDEL.chrAll.vcf.gz"
)
DEFAULT_LIMIT = 9_000_000_000
BGZF_MAX_BLOCK = 65_536
MAX_TRANSFER_ATTEMPTS = 5
TRANSFER_BACKOFF_SECONDS = 2.0
RETRYABLE_HTTP_STATUS = frozenset({429, 500, 502, 503, 504})
_USER_AGENT = "CANIS-indexed-vcf/1.0"


@dataclass
class ReferenceIndex:
    bins: dict[int, list[tuple[int, int]]]
    linear: list[int]


@dataclass
class TransferBudget:
    limit: int
    used: int = 0

    def reserve(self, amount: int) -> None:
        if amount < 0 or self.used + amount > self.limit:
            raise RuntimeError(
                f"transfer would exceed safety limit: {self.used + amount:,} > {self.limit:,} bytes"
            )

    def add(self, amount: int) -> None:
        self.reserve(amount)
        self.used += amount


def request(
    url: str,
    budget: TransferBudget,
    byte_range: tuple[int, int] | None = None,
    *,
    user_agent: str = _USER_AGENT,
) -> bytes:
    """HTTP GET with optional byte range, budget tracking, and transient retries."""
    headers = {"User-Agent": user_agent}
    if byte_range is not None:
        start, end = byte_range
        headers["Range"] = f"bytes={start}-{end}"
    http_request = urllib.request.Request(url, headers=headers)

    # The NIH host intermittently drops or times out a byte-range request. A single such
    # failure must not abort a long run, so transient network/HTTP errors are retried
    # with exponential backoff. Hard failures (byte-range ignored, transfer-budget
    # exceeded, non-retryable HTTP status) still propagate immediately.
    last_error: Exception | None = None
    for attempt in range(1, MAX_TRANSFER_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(http_request, timeout=120) as response:
                status = getattr(response, "status", response.getcode())
                if byte_range is not None and status != 206:
                    raise RuntimeError(
                        f"server ignored byte-range request for {url} (HTTP {status}); "
                        "aborting to prevent a full source download"
                    )
                declared = response.headers.get("Content-Length")
                if declared is not None:
                    budget.reserve(int(declared))
                data = response.read()
            budget.add(len(data))
            return data
        except urllib.error.HTTPError as exc:
            if exc.code not in RETRYABLE_HTTP_STATUS:
                raise
            last_error = exc
        except (urllib.error.URLError, http.client.IncompleteRead, OSError) as exc:
            last_error = exc
        if attempt < MAX_TRANSFER_ATTEMPTS:
            delay = TRANSFER_BACKOFF_SECONDS * (2 ** (attempt - 1))
            print(
                f"   ! transfer attempt {attempt}/{MAX_TRANSFER_ATTEMPTS} failed "
                f"({last_error}); retrying in {delay:.0f}s",
                flush=True,
            )
            time.sleep(delay)
    raise RuntimeError(
        f"transfer failed after {MAX_TRANSFER_ATTEMPTS} attempts for {url}: {last_error}"
    )


def download_small(url: str, path: Path, budget: TransferBudget) -> bytes:
    if path.exists() and path.stat().st_size:
        return path.read_bytes()
    data = request(url, budget)
    path.write_bytes(data)
    return data


def parse_tabix(raw_bgzf: bytes) -> tuple[list[str], list[ReferenceIndex]]:
    raw = gzip.decompress(raw_bgzf)
    if raw[:4] != b"TBI\x01":
        raise RuntimeError("source index is not a tabix TBI index")
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
        raise RuntimeError(f"tabix reference mismatch: {len(names)} names, {len(refs)} indexes")
    return names, refs


def reg2bins(start: int, end: int) -> list[int]:
    end -= 1
    bins = [0]
    bins.extend(range(1 + (start >> 26), 1 + (end >> 26) + 1))
    bins.extend(range(9 + (start >> 23), 9 + (end >> 23) + 1))
    bins.extend(range(73 + (start >> 20), 73 + (end >> 20) + 1))
    bins.extend(range(585 + (start >> 17), 585 + (end >> 17) + 1))
    bins.extend(range(4681 + (start >> 14), 4681 + (end >> 14) + 1))
    return bins


def chunks_for(index: ReferenceIndex, pos_1based: int) -> list[tuple[int, int]]:
    start = pos_1based - 1
    linear_i = start >> 14
    min_offset = index.linear[linear_i] if linear_i < len(index.linear) else 0
    chunks: list[tuple[int, int]] = []
    for bin_id in reg2bins(start, start + 1):
        chunks.extend(chunk for chunk in index.bins.get(bin_id, []) if chunk[1] > min_offset)
    return chunks


def merge_chunks(chunks: set[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(chunks)
    merged: list[list[int]] = []
    for begin, end in ordered:
        if not merged:
            merged.append([begin, end])
            continue
        previous = merged[-1]
        previous_byte_end = (previous[1] >> 16) + BGZF_MAX_BLOCK
        current_byte_start = begin >> 16
        if current_byte_start <= previous_byte_end:
            previous[1] = max(previous[1], end)
        else:
            merged.append([begin, end])
    return [(begin, end) for begin, end in merged]


def decompress_bgzf(data: bytes, first_uncompressed_offset: int = 0) -> bytes:
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


def source_header(url: str, budget: TransferBudget) -> tuple[list[str], list[str]]:
    fetch_size = 1_048_576
    while fetch_size <= 16_777_216:
        compressed = request(url, budget, (0, fetch_size - 1))
        text = decompress_bgzf(compressed).decode(errors="replace")
        lines = text.splitlines()
        column_line = next((line for line in lines if line.startswith("#CHROM\t")), None)
        if column_line is not None:
            return [line for line in lines if line.startswith("##")], column_line.split("\t")
        fetch_size *= 2
    raise RuntimeError("VCF header exceeded 16 MB compressed; refusing to continue")


def genotype(sample_field: str, gt_index: int) -> str | None:
    values = sample_field.split(":")
    if gt_index >= len(values):
        return None
    gt = values[gt_index]
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    return gt


# Backward-compatible private aliases used by older script call sites / tests.
_request = request
_parse_tabix = parse_tabix
_chunks_for = chunks_for
_merge_chunks = merge_chunks
_decompress_bgzf = decompress_bgzf
_source_header = source_header
_genotype = genotype
_download_small = download_small
_reg2bins = reg2bins
