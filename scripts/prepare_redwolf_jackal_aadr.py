"""Prepare the laptop Red Wolf/jackal panel without bcftools or a full VCF download.

The source VCF is bgzip-compressed and tabix-indexed. This script reads its small public
index, resolves the compressed chunks overlapping an evenly spaced CanineHD marker panel,
and downloads only those byte ranges. A hard transfer ceiling prevents an accidental full
309 GB response or unexpectedly large indexed extraction.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import http.client
import json
import struct
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

SOURCE_VCF_URL = (
    "https://research.nhgri.nih.gov/dog_genome/downloads/datasets/WGS/"
    "722g.990.SNP.INDEL.chrAll.vcf.gz"
)
PANEL_URL = (
    "https://download.cncb.ac.cn/dogsd/dog10k/variations/"
    "CFA31_IlluminaHD.vcf.gz"
)
PANEL_MD5 = "4a87088b17631fb6237210044ac099a6"
SAMPLES = ("Wolf25", "Wolf26", "GoldenJackal01")
DEFAULT_LIMIT = 9_000_000_000
BGZF_MAX_BLOCK = 65_536
MAX_TRANSFER_ATTEMPTS = 5
TRANSFER_BACKOFF_SECONDS = 2.0
RETRYABLE_HTTP_STATUS = frozenset({429, 500, 502, 503, 504})


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
                f"transfer would exceed safety limit: {self.used + amount:,} > "
                f"{self.limit:,} bytes"
            )

    def add(self, amount: int) -> None:
        self.reserve(amount)
        self.used += amount


def _request(url: str, budget: TransferBudget, byte_range: tuple[int, int] | None = None) -> bytes:
    headers = {"User-Agent": "CANIS-redwolf-jackal-panel/1.0"}
    if byte_range is not None:
        start, end = byte_range
        headers["Range"] = f"bytes={start}-{end}"
    request = urllib.request.Request(url, headers=headers)

    # The NIH host intermittently drops or times out a byte-range request. A single such
    # failure must not abort an ~80-minute run, so transient network/HTTP errors are retried
    # with exponential backoff. Hard failures (byte-range ignored, transfer-budget exceeded,
    # non-retryable HTTP status) still propagate immediately.
    last_error: Exception | None = None
    for attempt in range(1, MAX_TRANSFER_ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                status = getattr(response, "status", response.getcode())
                if byte_range is not None and status != 206:
                    raise RuntimeError(
                        f"server ignored byte-range request for {url} (HTTP {status}); aborting "
                        "to prevent a full source download"
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


def _download_small(url: str, path: Path, budget: TransferBudget) -> bytes:
    if path.exists() and path.stat().st_size:
        return path.read_bytes()
    data = _request(url, budget)
    path.write_bytes(data)
    return data


def _parse_tabix(raw_bgzf: bytes) -> tuple[list[str], list[ReferenceIndex]]:
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
    chunks = []
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
        previous_byte_end = (previous[1] >> 16) + BGZF_MAX_BLOCK
        current_byte_start = begin >> 16
        if current_byte_start <= previous_byte_end:
            previous[1] = max(previous[1], end)
        else:
            merged.append([begin, end])
    return [(begin, end) for begin, end in merged]


def _decompress_bgzf(data: bytes, first_uncompressed_offset: int = 0) -> bytes:
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
        raise RuntimeError("no CanFam3.1 autosomal positions found in CanineHD panel")
    stride = max(1, (len(rows) + target_sites - 1) // target_sites)
    return rows[::stride]


def _source_header(url: str, budget: TransferBudget) -> tuple[list[str], list[str]]:
    fetch_size = 1_048_576
    while fetch_size <= 16_777_216:
        compressed = _request(url, budget, (0, fetch_size - 1))
        text = _decompress_bgzf(compressed).decode(errors="replace")
        lines = text.splitlines()
        column_line = next((line for line in lines if line.startswith("#CHROM\t")), None)
        if column_line is not None:
            return [line for line in lines if line.startswith("##")], column_line.split("\t")
        fetch_size *= 2
    raise RuntimeError("VCF header exceeded 16 MB compressed; refusing to continue")


def _genotype(sample_field: str, gt_index: int) -> str | None:
    values = sample_field.split(":")
    if gt_index >= len(values):
        return None
    gt = values[gt_index]
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    return gt


def prepare(args: argparse.Namespace) -> Path:
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    budget = TransferBudget(args.max_download_bytes)
    panel_path = out_dir / "CFA31_IlluminaHD.vcf.gz"
    index_path = out_dir / "722g.990.SNP.INDEL.chrAll.vcf.gz.tbi"
    output_path = out_dir / f"redwolf_jackal_caninehd{args.target_sites}.vcf.gz"
    manifest_path = output_path.with_suffix("").with_suffix(".manifest.json")

    print(">> downloading small marker panel and tabix index", flush=True)
    panel_data = _download_small(PANEL_URL, panel_path, budget)
    # MD5 is the publisher-provided file-integrity value, not a security primitive.
    if hashlib.md5(panel_data).hexdigest() != PANEL_MD5:
        raise RuntimeError("CanineHD panel MD5 mismatch")
    index_data = _download_small(f"{SOURCE_VCF_URL}.tbi", index_path, budget)
    names, indexes = _parse_tabix(index_data)
    name_to_index = {name: index for name, index in zip(names, indexes, strict=True)}

    positions = _panel_positions(panel_data, args.target_sites)
    mapped_positions: list[tuple[str, int]] = []
    for chrom, pos in positions:
        source_chrom = chrom if chrom in name_to_index else f"chr{chrom}"
        if source_chrom in name_to_index:
            mapped_positions.append((source_chrom, pos))
    targets = set(mapped_positions)
    if not targets:
        raise RuntimeError("CanineHD chromosome names do not match the 722-genome index")

    chunks: set[tuple[int, int]] = set()
    for chrom, pos in targets:
        chunks.update(_chunks_for(name_to_index[chrom], pos))
    merged = _merge_chunks(chunks)
    estimated = sum(((end >> 16) + BGZF_MAX_BLOCK) - (begin >> 16) for begin, end in merged)
    print(
        f">> {len(targets):,} targets resolve to {len(merged):,} byte ranges; "
        f"estimated transfer {estimated / 1_000_000_000:.2f} GB",
        flush=True,
    )
    if budget.used + estimated > budget.limit:
        raise RuntimeError(
            f"indexed extraction estimate exceeds safety limit: "
            f"{budget.used + estimated:,} > {budget.limit:,} bytes"
        )

    meta_lines, columns = _source_header(SOURCE_VCF_URL, budget)
    missing_samples = [sample for sample in SAMPLES if sample not in columns]
    if missing_samples:
        raise RuntimeError(f"samples absent from source VCF: {missing_samples}")
    sample_columns = [columns.index(sample) for sample in SAMPLES]

    records: dict[tuple[str, int], str] = {}
    for number, (virtual_begin, virtual_end) in enumerate(merged, start=1):
        byte_start = virtual_begin >> 16
        byte_end = (virtual_end >> 16) + BGZF_MAX_BLOCK - 1
        data = _request(SOURCE_VCF_URL, budget, (byte_start, byte_end))
        text = _decompress_bgzf(data, virtual_begin & 0xFFFF).decode(errors="replace")
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) <= max(sample_columns):
                continue
            try:
                key = (fields[0], int(fields[1]))
            except (ValueError, IndexError):
                continue
            if key not in targets or key in records:
                continue
            ref, alt = fields[3], fields[4]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            formats = fields[8].split(":")
            if "GT" not in formats:
                continue
            gt_index = formats.index("GT")
            genotypes = [_genotype(fields[index], gt_index) for index in sample_columns]
            if any(gt is None for gt in genotypes):
                continue
            output = fields[:8] + ["GT"] + [str(gt) for gt in genotypes]
            records[key] = "\t".join(output)
        if number % 100 == 0 or number == len(merged):
            print(
                f"   ranges {number:,}/{len(merged):,}; "
                f"downloaded {budget.used / 1_000_000_000:.2f} GB; "
                f"retained {len(records):,} SNPs",
                flush=True,
            )

    if len(records) < 1000:
        raise RuntimeError(f"only {len(records)} complete biallelic SNPs retained; expected >=1000")
    with gzip.open(output_path, "wt", encoding="utf-8", compresslevel=6) as handle:
        for line in meta_lines:
            handle.write(f"{line}\n")
        handle.write("##CANIS_reduced_panel=CFA31_IlluminaHD_even_subset\n")
        handle.write("\t".join(columns[:9] + list(SAMPLES)) + "\n")
        for key in sorted(records, key=lambda item: (int(item[0].removeprefix("chr")), item[1])):
            handle.write(f"{records[key]}\n")

    digest = hashlib.sha256(output_path.read_bytes()).hexdigest()
    manifest = {
        "dataset_id": "redwolf_jackal_722g_caninehd",
        "description": "AADR-style reduced CanFam3.1 SNP panel; not raw genomes.",
        "created_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "source": {
            "bioproject": "PRJNA448733",
            "analysis_accession": "SRZ189891",
            "vcf_url": SOURCE_VCF_URL,
            "reference": "CanFam3.1",
        },
        "samples": [
            {"sample_id": "Wolf25", "taxon": "red_wolf", "biosample": "SAMN02921317"},
            {"sample_id": "Wolf26", "taxon": "red_wolf", "biosample": "SAMN02921318"},
            {
                "sample_id": "GoldenJackal01",
                "taxon": "golden_jackal",
                "biosample": "SAMN03366713",
            },
        ],
        "marker_panel": {
            "name": "CFA31_IlluminaHD",
            "requested_sites": args.target_sites,
            "targeted_sites": len(targets),
            "retained_sites": len(records),
        },
        "transfer": {"downloaded_bytes": budget.used, "hard_limit_bytes": budget.limit},
        "output": {"path": str(output_path), "sha256": digest},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f">> dataset ready: {output_path}", flush=True)
    print(f">> total downloaded: {budget.used / 1_000_000_000:.2f} GB", flush=True)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target-sites", type=int, default=25_000)
    parser.add_argument("--out-dir", type=Path, default=Path("data/redwolf_jackal_aadr"))
    parser.add_argument("--max-download-bytes", type=int, default=DEFAULT_LIMIT)
    args = parser.parse_args()
    if args.target_sites < 1000:
        parser.error("--target-sites must be at least 1000")
    if not 1 <= args.max_download_bytes < 10_000_000_000:
        parser.error("--max-download-bytes must be positive and below 10 GB")
    return args


if __name__ == "__main__":
    try:
        prepare(parse_args())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
