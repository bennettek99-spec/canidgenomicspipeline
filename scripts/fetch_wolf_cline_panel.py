"""Fetch an unascertained window panel for the red / Great Lakes / eastern wolf study.

Rather than a fixed dog-array SNP list (which is ascertained in dogs), this reads
*contiguous* stretches of the indexed NHGRI 722-genome VCF at evenly spaced points
across the 38 autosomes, so every variant called in those windows is kept. That
makes the panel suitable for f-statistics and per-window heterozygosity.

Only byte ranges are transferred; the budget is hard-capped. Output (gitignored):

    data/wolf_cline/panel.npz      genotypes (alt-allele count, -1 missing), GQ,
                                   chrom/pos, window id and window spans

Run:
    python scripts/fetch_wolf_cline_panel.py --windows 400 --window-bytes 3500000
"""

from __future__ import annotations

import argparse
import csv
import time
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from canidae.io.indexed_vcf import (
    SOURCE_VCF_URL,
    TransferBudget,
    chunks_for,
    decompress_bgzf,
    download_small,
    parse_tabix,
    request,
    source_header,
)

ROOT = Path(__file__).resolve().parents[1]
SAMPLES_CSV = ROOT / "configs" / "examples" / "wolf_cline_samples.csv"
OUT_DIR = ROOT / "data" / "wolf_cline"
AUTOSOMES = [f"chr{i}" for i in range(1, 39)]


def load_sample_ids(path: Path = SAMPLES_CSV) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [row["sample_id"] for row in csv.DictReader(handle)]


def contig_lengths(meta_lines: list[str]) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line in meta_lines:
        if line.startswith("##contig=<"):
            fields = dict(item.split("=", 1) for item in line[10:-1].split(",") if "=" in item)
            if "ID" in fields and "length" in fields:
                lengths[fields["ID"]] = int(fields["length"])
    return lengths


def window_starts(lengths: dict[str, int], n_windows: int) -> list[tuple[str, int]]:
    """Evenly spaced genomic anchor points, proportional to autosome length."""
    total = sum(lengths[c] for c in AUTOSOMES)
    step = total / n_windows
    anchors: list[tuple[str, int]] = []
    offset = step / 2
    for contig in AUTOSOMES:
        length = lengths[contig]
        while offset < length:
            anchors.append((contig, int(offset) + 1))
            offset += step
        offset -= length
    return anchors


def parse_window(
    text: str, contig: str, columns: list[int]
) -> tuple[list[tuple[int, list[int], list[int]]], int, int]:
    """Biallelic PASS SNPs on *contig* plus the span of all records seen."""
    records: list[tuple[int, list[int], list[int]]] = []
    first = last = 0
    lines = text.split("\n")[:-1]  # the final line is truncated by the byte range
    for line in lines:
        fields = line.split("\t")
        if len(fields) < 10 or fields[0] != contig:
            continue
        pos = int(fields[1])
        first = first or pos
        last = pos
        if fields[6] != "PASS" or len(fields[3]) != 1 or len(fields[4]) != 1:
            continue
        fmt = fields[8].split(":")
        gt_i, gq_i = fmt.index("GT"), fmt.index("GQ") if "GQ" in fmt else -1
        gts: list[int] = []
        gqs: list[int] = []
        for column in columns:
            values = fields[column].split(":")
            alleles = values[gt_i].replace("|", "/").split("/")
            if len(alleles) == 2 and all(a in ("0", "1") for a in alleles):
                gts.append(int(alleles[0]) + int(alleles[1]))
            else:
                gts.append(-1)
            try:
                gqs.append(min(int(values[gq_i]), 99) if gq_i >= 0 else 0)
            except (IndexError, ValueError):
                gqs.append(0)
        if any(g > 0 for g in gts):  # drop sites monomorphic-reference in this cohort
            records.append((pos, gts, gqs))
    return records, first, last


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", type=int, default=400)
    parser.add_argument("--window-bytes", type=int, default=3_500_000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-bytes", type=int, default=2_500_000_000)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    budget = TransferBudget(args.max_bytes)
    names, refs = parse_tabix(
        download_small(SOURCE_VCF_URL + ".tbi", OUT_DIR / "source.vcf.gz.tbi", budget)
    )
    meta, header = source_header(SOURCE_VCF_URL, budget)
    samples = load_sample_ids()
    missing = [s for s in samples if s not in header]
    if missing:
        raise SystemExit(f"samples absent from source VCF: {missing}")
    columns = [header.index(s) for s in samples]

    anchors = window_starts(contig_lengths(meta), args.windows)
    plans: list[tuple[int, str, int]] = []
    for window_id, (contig, pos) in enumerate(anchors):
        chunks = chunks_for(refs[names.index(contig)], pos)
        if chunks:
            plans.append((window_id, contig, min(c[0] for c in chunks)))
    print(f"{len(plans)} windows x {args.window_bytes / 1e6:.1f} MB", flush=True)

    Parsed = tuple[int, str, list[tuple[int, list[int], list[int]]], int, int]

    def fetch(plan: tuple[int, str, int]) -> Parsed:
        # Parse inside the worker so only small per-window SNP lists are held in memory.
        window_id, contig, virtual = plan
        start = virtual >> 16
        data = request(SOURCE_VCF_URL, budget, (start, start + args.window_bytes - 1))
        text = decompress_bgzf(data, virtual & 0xFFFF).decode()
        return (window_id, contig, *parse_window(text, contig, columns))

    def results() -> Iterator[Parsed]:
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            yield from pool.map(fetch, plans)

    chrom: list[int] = []
    pos_out: list[int] = []
    win_out: list[int] = []
    gt_rows: list[list[int]] = []
    gq_rows: list[list[int]] = []
    spans: list[tuple[int, int, int, int]] = []
    started = time.time()
    for done, (window_id, contig, records, first, last) in enumerate(results(), 1):
        chrom_i = AUTOSOMES.index(contig) + 1
        spans.append((window_id, chrom_i, first, last))
        for pos, gts, gqs in records:
            chrom.append(chrom_i)
            pos_out.append(pos)
            win_out.append(window_id)
            gt_rows.append(gts)
            gq_rows.append(gqs)
        if done % 20 == 0 or done == len(plans):
            rate = budget.used / (time.time() - started) / 1e6
            print(
                f"  {done}/{len(plans)} windows, {len(pos_out):,} SNPs, "
                f"{budget.used / 1e9:.2f} GB at {rate:.1f} MB/s",
                flush=True,
            )

    np.savez_compressed(
        OUT_DIR / "panel.npz",
        samples=np.array(samples),
        chrom=np.array(chrom, dtype=np.int8),
        pos=np.array(pos_out, dtype=np.int64),
        window=np.array(win_out, dtype=np.int32),
        gt=np.array(gt_rows, dtype=np.int8),
        gq=np.array(gq_rows, dtype=np.uint8),
        spans=np.array(spans, dtype=np.int64),
        bytes_transferred=np.array(budget.used),
        source=np.array(SOURCE_VCF_URL),
    )
    print(f"wrote {OUT_DIR / 'panel.npz'} ({len(pos_out):,} SNPs)")


if __name__ == "__main__":
    main()
