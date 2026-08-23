"""Prepare the laptop Red Wolf/jackal panel without bcftools or a full VCF download.

The source VCF is bgzip-compressed and tabix-indexed. This script reads its small public
index, resolves the compressed chunks overlapping an evenly spaced CanineHD marker panel,
and downloads only those byte ranges. A hard transfer ceiling prevents an accidental full
309 GB response or unexpectedly large indexed extraction.

Transport helpers live in :mod:`canidae.io.indexed_vcf`; this module keeps the
panel-selection orchestration and CLI.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

from canidae.io import indexed_vcf
from canidae.io.indexed_vcf import (
    BGZF_MAX_BLOCK,
    DEFAULT_LIMIT,
    SOURCE_VCF_URL,
    TransferBudget,
    chunks_for,
    decompress_bgzf,
    download_small,
    genotype,
    merge_chunks,
    parse_tabix,
    request,
    source_header,
)

# Re-export private names so older call sites / tests that import this module as
# ``indexed_vcf`` keep working during the library migration.
_request = request
_parse_tabix = parse_tabix
_chunks_for = chunks_for
_merge_chunks = merge_chunks
_decompress_bgzf = decompress_bgzf
_source_header = source_header
_genotype = genotype
_download_small = download_small
ReferenceIndex = indexed_vcf.ReferenceIndex

PANEL_URL = "https://download.cncb.ac.cn/dogsd/dog10k/variations/CFA31_IlluminaHD.vcf.gz"
PANEL_MD5 = "4a87088b17631fb6237210044ac099a6"
SAMPLES = ("Wolf25", "Wolf26", "GoldenJackal01")


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


def prepare(args: argparse.Namespace) -> Path:
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    budget = TransferBudget(args.max_download_bytes)
    panel_path = out_dir / "CFA31_IlluminaHD.vcf.gz"
    index_path = out_dir / "722g.990.SNP.INDEL.chrAll.vcf.gz.tbi"
    output_path = out_dir / f"redwolf_jackal_caninehd{args.target_sites}.vcf.gz"
    manifest_path = output_path.with_suffix("").with_suffix(".manifest.json")

    print(">> downloading small marker panel and tabix index", flush=True)
    panel_data = download_small(PANEL_URL, panel_path, budget)
    # MD5 is the publisher-provided file-integrity value, not a security primitive.
    if hashlib.md5(panel_data).hexdigest() != PANEL_MD5:
        raise RuntimeError("CanineHD panel MD5 mismatch")
    index_data = download_small(f"{SOURCE_VCF_URL}.tbi", index_path, budget)
    names, indexes = parse_tabix(index_data)
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
        chunks.update(chunks_for(name_to_index[chrom], pos))
    merged = merge_chunks(chunks)
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

    meta_lines, columns = source_header(SOURCE_VCF_URL, budget)
    missing_samples = [sample for sample in SAMPLES if sample not in columns]
    if missing_samples:
        raise RuntimeError(f"samples absent from source VCF: {missing_samples}")
    sample_columns = [columns.index(sample) for sample in SAMPLES]

    records: dict[tuple[str, int], str] = {}
    for number, (virtual_begin, virtual_end) in enumerate(merged, start=1):
        byte_start = virtual_begin >> 16
        byte_end = (virtual_end >> 16) + BGZF_MAX_BLOCK - 1
        data = request(SOURCE_VCF_URL, budget, (byte_start, byte_end))
        text = decompress_bgzf(data, virtual_begin & 0xFFFF).decode(errors="replace")
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
            genotypes = [genotype(fields[index], gt_index) for index in sample_columns]
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
