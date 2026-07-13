"""Build a laptop-safe, cross-source bridge panel for the eastern-coyote report.

The eastern-coyote study deposit is a RADseq PED/MAP matrix, while the earlier
proxy analysis uses the NHGRI 722-canid WGS VCF.  These formats must not be
concatenated directly.  This tool instead selects RADseq loci, obtains the exact
WGS records at those coordinates through tabix byte ranges, verifies each allele
pair, reorients the RADseq calls to the WGS REF/ALT alleles, and writes a merged
VCF for exploratory population-level comparison.

The output is a *cross-platform bridge panel*, not a whole-genome joint callset.
It is appropriate for an exploratory PCA/FST/distance/report comparison, but not
for a standalone estimate of wolf ancestry or reference-allele interpretation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import http.client
import json
import socket
import time
import urllib.parse
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from canidae.stages.acquisition import reduced_panel as rp
from canidae.stages.processing.vcf_io import read_biallelic_snps


SOURCE_VCF_URL = (
    "https://research.nhgri.nih.gov/dog_genome/downloads/datasets/WGS/"
    "722g.990.SNP.INDEL.chrAll.vcf.gz"
)
WGS_SAMPLES = (
    ("Coyote01", "coyote", "coyote_wgs", "north_america"),
    ("Coyote02", "coyote", "coyote_wgs", "north_america"),
    ("AlaskanWolf", "gray_wolf", "gray_wolf_wgs", "north_america"),
    ("AlgonquinWolf13467", "eastern_wolf", "eastern_wolf_wgs", "north_america"),
    ("AlgonquinWolf13470", "eastern_wolf", "eastern_wolf_wgs", "north_america"),
    ("GoldenJackal01", "golden_jackal", "golden_jackal_wgs", "",),
)


def prefer_ipv4() -> None:
    """Avoid a known stalled IPv6 route to the public NHGRI host on this laptop."""
    original = socket.getaddrinfo

    def ipv4_only(host, port, family=0, type=0, proto=0, flags=0):
        answers = original(host, port, family, type, proto, flags)
        ipv4 = [answer for answer in answers if answer[0] == socket.AF_INET]
        return ipv4 or answers

    socket.getaddrinfo = ipv4_only


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--eastern-vcf",
        type=Path,
        default=Path(
            "data/eastern_coyote_radseq_workspace/store/qc/qc_filtered.vcf"
        ),
        help="QC-filtered eastern-coyote VCF created by the RADseq run.",
    )
    parser.add_argument(
        "--eastern-sample-sheet",
        type=Path,
        default=Path(
            "data/eastern_coyote_radseq_workspace/store/qc/qc_sample_sheet.csv"
        ),
        help="QC-filtered sample sheet that matches --eastern-vcf.",
    )
    parser.add_argument("--source-vcf-url", default=SOURCE_VCF_URL)
    parser.add_argument(
        "--source-header-vcf",
        type=Path,
        default=Path(
            "data/coywolf_laptop_workspace/store/reduced_panel/reduced_panel_2k.vcf.gz"
        ),
        help=(
            "Existing local subset of the same source VCF, used only to obtain the "
            "#CHROM sample header without an extra remote request."
        ),
    )
    parser.add_argument(
        "--index-cache",
        type=Path,
        default=Path("data/eastern_coyote_wgs_bridge/722g.990.SNP.INDEL.chrAll.vcf.gz.tbi"),
        help="Small cached tabix index. It is downloaded once if absent.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/eastern_coyote_wgs_bridge"),
        help="Directory for the bridge VCF, sample sheet, and manifest.",
    )
    parser.add_argument(
        "--target-sites",
        type=int,
        default=1000,
        help="Uniformly sampled autosomal RADseq loci to query from WGS.",
    )
    parser.add_argument("--min-retained-sites", type=int, default=250)
    parser.add_argument("--max-download-bytes", type=int, default=9_000_000_000)
    parser.add_argument("--confirmation-threshold-bytes", type=int, default=1_000_000_000)
    parser.add_argument(
        "--confirm-large-transfer",
        action="store_true",
        help="Permit a preflighted transfer above the confirmation threshold.",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Calculate the exact indexed-range estimate but do not download WGS records.",
    )
    return parser.parse_args()


def _choose_targets(eastern, target_sites: int):
    autosomal = [
        index
        for index, chrom in enumerate(eastern.chrom)
        if str(chrom).isdigit() and 1 <= int(str(chrom)) <= 38
    ]
    if not autosomal:
        raise ValueError("no numeric autosomal loci found in the eastern-coyote VCF")
    if not 1 <= target_sites <= len(autosomal):
        raise ValueError(
            f"--target-sites must be between 1 and {len(autosomal):,}; got {target_sites:,}"
        )
    return [autosomal[(number * len(autosomal)) // target_sites] for number in range(target_sites)]


def _load_or_fetch_index(args: argparse.Namespace) -> tuple[bytes, int]:
    if args.index_cache.exists():
        raw = args.index_cache.read_bytes()
        return raw, 0
    args.index_cache.parent.mkdir(parents=True, exist_ok=True)
    budget = rp.TransferBudget(args.max_download_bytes)
    raw = rp._request(
        f"{args.source_vcf_url}.tbi",
        budget,
        timeout_seconds=120,
    )
    args.index_cache.write_bytes(raw)
    return raw, budget.used


def _estimate(index_data: bytes, selected, *, target_sites: int) -> dict:
    names, indexes = rp._parse_tabix(index_data)
    by_contig = dict(zip(names, indexes, strict=True))
    source_contigs = set(by_contig)
    targets: dict[tuple[str, int], int] = {}
    for index in selected.indices:
        source_chrom = rp._resolve_contig(str(selected.eastern.chrom[index]), source_contigs)
        if source_chrom is not None:
            targets[(source_chrom, int(selected.eastern.pos[index]))] = index
    if not targets:
        raise ValueError("none of the selected RADseq contigs match the WGS tabix index")
    chunks: set[tuple[int, int]] = set()
    for chrom, position in targets:
        chunks.update(rp._chunks_for(by_contig[chrom], position))
    merged = rp._merge_chunks(chunks)
    range_bytes = sum(
        ((end >> 16) + rp._BGZF_MAX_BLOCK) - (begin >> 16) for begin, end in merged
    )
    return {
        "names": names,
        "by_contig": by_contig,
        "targets": targets,
        "merged": merged,
        "transfer": rp.TransferEstimate(
            prerequisite_bytes=len(index_data),
            range_bytes=range_bytes,
            header_allowance_bytes=16 * 1024 * 1024,
            n_ranges=len(merged),
            n_target_sites=len(targets),
        ),
        "requested_target_sites": target_sites,
    }


class SelectedEastern:
    def __init__(self, eastern, indices: list[int]):
        self.eastern = eastern
        self.indices = indices


class PersistentRangeClient:
    """Strict byte-range client that reuses one HTTPS connection across tabix chunks."""

    def __init__(self, url: str, budget: rp.TransferBudget, *, timeout_seconds: int = 120):
        parsed = urllib.parse.urlsplit(url)
        if parsed.scheme != "https" or not parsed.hostname:
            raise ValueError(f"bridge source must be an HTTPS URL: {url}")
        self.host = parsed.hostname
        self.path = urllib.parse.urlunsplit(("", "", parsed.path, parsed.query, ""))
        self.port = parsed.port or 443
        self.budget = budget
        self.timeout_seconds = timeout_seconds
        self.connection: http.client.HTTPSConnection | None = None

    def close(self) -> None:
        if self.connection is not None:
            self.connection.close()
            self.connection = None

    def get(self, byte_range: tuple[int, int], progress) -> bytes:
        start, end = byte_range
        expected_max = end - start + 1
        last_error: Exception | None = None
        for attempt in range(1, 6):
            try:
                if self.connection is None:
                    self.connection = http.client.HTTPSConnection(
                        self.host, self.port, timeout=self.timeout_seconds
                    )
                self.connection.request(
                    "GET",
                    self.path,
                    headers={
                        "Range": f"bytes={start}-{end}",
                        "User-Agent": "CANIS-eastern-coyote-WGS-bridge/0.1",
                        "Connection": "keep-alive",
                    },
                )
                response = self.connection.getresponse()
                if response.status != 206:
                    raise RuntimeError(
                        f"source ignored byte range {start}-{end} (HTTP {response.status})"
                    )
                declared = response.getheader("Content-Length")
                if declared is not None:
                    declared_bytes = int(declared)
                    if declared_bytes > expected_max:
                        raise RuntimeError(
                            f"source response exceeds requested byte range: {declared_bytes} > "
                            f"{expected_max}"
                        )
                    self.budget.reserve(declared_bytes)
                content_range = response.getheader("Content-Range") or ""
                if not content_range.startswith(f"bytes {start}-"):
                    raise RuntimeError(f"unexpected Content-Range for request {start}-{end}: {content_range}")
                data = response.read()
                if len(data) > expected_max:
                    raise RuntimeError(
                        f"source body exceeds requested byte range: {len(data)} > {expected_max}"
                    )
                self.budget.add(len(data))
                return data
            except Exception as exc:
                last_error = exc
                self.close()
                if attempt < 5:
                    delay = 2 ** (attempt - 1)
                    progress(
                        f"range request attempt {attempt}/5 failed ({exc}); retrying in {delay}s"
                    )
                    time.sleep(delay)
        raise RuntimeError(f"range request failed after 5 attempts: {last_error}")


def _source_records(
    args: argparse.Namespace,
    plan: dict,
    index_downloaded_bytes: int,
) -> tuple[dict[tuple[str, int], tuple[list[str], list[str]]], list[str], int, Counter[str]]:
    def progress(message: str) -> None:
        print(message, flush=True)

    budget = rp.TransferBudget(args.max_download_bytes)
    if index_downloaded_bytes:
        budget.add(index_downloaded_bytes)
    columns = _source_columns(args, budget)
    sample_ids = [row[0] for row in WGS_SAMPLES]
    missing = [sample for sample in sample_ids if sample not in columns]
    if missing:
        raise ValueError("WGS samples absent from source header: " + ", ".join(missing))
    sample_columns = [columns.index(sample) for sample in sample_ids]
    records: dict[tuple[str, int], tuple[list[str], list[str]]] = {}
    skipped: Counter[str] = Counter()
    client = PersistentRangeClient(args.source_vcf_url, budget)
    try:
        for number, (virtual_begin, virtual_end) in enumerate(plan["merged"], start=1):
            byte_start = virtual_begin >> 16
            byte_end = (virtual_end >> 16) + rp._BGZF_MAX_BLOCK - 1
            data = client.get((byte_start, byte_end), progress)
            text = rp._decompress_bgzf(data, virtual_begin & 0xFFFF).decode(errors="replace")
            for line in text.splitlines():
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) <= max(sample_columns):
                    skipped["truncated_record"] += 1
                    continue
                try:
                    key = (fields[0], int(fields[1]))
                except (IndexError, ValueError):
                    skipped["malformed_position"] += 1
                    continue
                if key not in plan["targets"] or key in records:
                    continue
                if len(fields) < 10 or len(fields[3]) != 1 or len(fields[4]) != 1 or "," in fields[4]:
                    skipped["not_biallelic_snp"] += 1
                    continue
                formats = fields[8].split(":")
                if "GT" not in formats:
                    skipped["missing_gt_field"] += 1
                    continue
                gt_index = formats.index("GT")
                genotypes = [rp._genotype(fields[column], gt_index) for column in sample_columns]
                if any(genotype is None for genotype in genotypes):
                    skipped["missing_wgs_gt"] += 1
                    continue
                records[key] = (fields[:5], [str(genotype) for genotype in genotypes])
            if number % 100 == 0 or number == len(plan["merged"]):
                print(
                    f"Downloaded WGS ranges {number:,}/{len(plan['merged']):,}; "
                    f"{rp.format_bytes(budget.used)}; exact-target records {len(records):,}",
                    flush=True,
                )
    finally:
        client.close()
    return records, columns, budget.used, skipped


def _source_columns(args: argparse.Namespace, budget: rp.TransferBudget) -> list[str]:
    """Use the local proxy subset header when available; otherwise fetch the source header."""
    if args.source_header_vcf.exists():
        with gzip.open(args.source_header_vcf, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#CHROM\t"):
                    return line.rstrip("\n").split("\t")
        raise ValueError(f"#CHROM header not found in {args.source_header_vcf}")
    _meta, columns = rp._source_header(
        args.source_vcf_url, budget, timeout_seconds=120, progress=None
    )
    return columns


def _source_cache_path(args: argparse.Namespace) -> Path:
    return args.out_dir / f"eastern_coyote_wgs_bridge_{args.target_sites}.source_records.json"


def _write_source_cache(
    path: Path,
    args: argparse.Namespace,
    source_records: dict[tuple[str, int], tuple[list[str], list[str]]],
    source_skips: Counter[str],
) -> None:
    payload = {
        "schema_version": 1,
        "source_vcf_url": args.source_vcf_url,
        "target_sites": args.target_sites,
        "records": [
            {
                "chrom": chrom,
                "position": position,
                "fields": fields,
                "genotypes": genotypes,
            }
            for (chrom, position), (fields, genotypes) in sorted(source_records.items())
        ],
        "source_record_skips": dict(sorted(source_skips.items())),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def _load_source_cache(
    path: Path, args: argparse.Namespace
) -> tuple[dict[tuple[str, int], tuple[list[str], list[str]]], Counter[str]] | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1 or payload.get("source_vcf_url") != args.source_vcf_url:
        return None
    if payload.get("target_sites") != args.target_sites:
        return None
    records: dict[tuple[str, int], tuple[list[str], list[str]]] = {}
    for row in payload.get("records", []):
        records[(str(row["chrom"]), int(row["position"]))] = (
            list(row["fields"]), list(row["genotypes"])
        )
    return records, Counter(payload.get("source_record_skips", {}))


def _remap_eastern_gt(genotypes, east_ref: str, east_alt: str, wgs_ref: str, wgs_alt: str) -> list[str]:
    allele_map = {
        0: "0" if east_ref == wgs_ref else "1",
        1: "0" if east_alt == wgs_ref else "1",
    }
    output: list[str] = []
    for call in genotypes:
        left, right = int(call[0]), int(call[1])
        if left < 0 or right < 0:
            output.append("./.")
        else:
            output.append(f"{allele_map[left]}/{allele_map[right]}")
    return output


def _write_bridge_vcf(
    destination: Path,
    selected: SelectedEastern,
    plan: dict,
    source_records: dict[tuple[str, int], tuple[list[str], list[str]]],
) -> tuple[int, Counter[str]]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    skips: Counter[str] = Counter()
    retained = 0
    sample_ids = [str(sample) for sample in selected.eastern.samples] + [row[0] for row in WGS_SAMPLES]
    with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=CANIS_eastern_coyote_WGS_bridge_panel\n")
        handle.write(f"##source_wgs={SOURCE_VCF_URL}\n")
        handle.write("##source_radseq_doi=10.5061/dryad.7f9q2cd\n")
        handle.write("##CANIS_bridge_panel=RADseq loci reoriented to matching NHGRI WGS alleles\n")
        handle.write(
            "##CANIS_warning=Cross-platform and ascertainment-limited exploratory panel; "
            "not a whole-genome joint callset\n"
        )
        handle.write(
            "##INFO=<ID=BRIDGE,Number=1,Type=String,Description=\"Coordinate and "
            "unordered-allele match between RADseq and WGS\">\n"
        )
        handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Unphased genotype\">\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(sample_ids) + "\n")
        for key in sorted(plan["targets"], key=lambda item: (int(item[0].removeprefix("chr")), item[1])):
            east_index = plan["targets"][key]
            source = source_records.get(key)
            if source is None:
                skips["no_complete_wgs_record"] += 1
                continue
            fields, source_gt = source
            east_ref = str(selected.eastern.ref[east_index])
            east_alt = str(selected.eastern.alt[east_index])
            wgs_ref, wgs_alt = fields[3], fields[4]
            if {east_ref, east_alt} != {wgs_ref, wgs_alt}:
                skips["allele_pair_mismatch"] += 1
                continue
            eastern_gt = _remap_eastern_gt(
                selected.eastern.gt[east_index], east_ref, east_alt, wgs_ref, wgs_alt
            )
            handle.write(
                f"{key[0]}\t{key[1]}\t{fields[2]}\t{wgs_ref}\t{wgs_alt}\t.\tPASS\t"
                "BRIDGE=coordinate_and_allele_pair_match\tGT\t"
            )
            handle.write("\t".join(eastern_gt + source_gt) + "\n")
            retained += 1
    temporary.replace(destination)
    return retained, skips


def _write_sample_sheet(destination: Path, eastern_sheet: Path, eastern_sample_ids: list[str]) -> None:
    with eastern_sheet.open(newline="", encoding="utf-8") as handle:
        eastern_rows = list(csv.DictReader(handle))
    rows = [row for row in eastern_rows if row["sample_id"] in set(eastern_sample_ids)]
    found = {row["sample_id"] for row in rows}
    missing = [sample for sample in eastern_sample_ids if sample not in found]
    if missing:
        raise ValueError("eastern sample metadata missing: " + ", ".join(missing))
    fields = ["sample_id", "taxon", "population", "region", "country", "source_study", "accession"]
    for sample_id, taxon, population, region in WGS_SAMPLES:
        rows.append({
            "sample_id": sample_id,
            "taxon": taxon,
            "population": population,
            "region": region,
            "country": "",
            "source_study": "Plassais2019_PRJNA448733_NHGRI722",
            "accession": "PRJNA448733",
        })
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{field: row.get(field, "") for field in fields} for row in rows])


def main() -> int:
    args = parse_args()
    prefer_ipv4()
    if args.min_retained_sites < 1:
        raise ValueError("--min-retained-sites must be positive")
    eastern = read_biallelic_snps(args.eastern_vcf)
    selected = SelectedEastern(eastern, _choose_targets(eastern, args.target_sites))
    index_data, index_downloaded_bytes = _load_or_fetch_index(args)
    plan = _estimate(index_data, selected, target_sites=args.target_sites)
    estimate = plan["transfer"]
    print(
        "Preflight: "
        f"{estimate.n_target_sites:,} target loci, {estimate.n_ranges:,} byte ranges, "
        f"estimated {rp.format_bytes(estimate.total_bytes)}.",
        flush=True,
    )
    rp.enforce_transfer_policy(
        estimate,
        ceiling_bytes=args.max_download_bytes,
        confirmation_threshold_bytes=args.confirmation_threshold_bytes,
        confirmed=args.confirm_large_transfer,
    )
    if args.preflight_only:
        return 0

    source_cache = _source_cache_path(args)
    cached = _load_source_cache(source_cache, args)
    if cached is None:
        source_records, _columns, downloaded_bytes, source_skips = _source_records(
            args, plan, index_downloaded_bytes
        )
        _write_source_cache(source_cache, args, source_records, source_skips)
    else:
        source_records, source_skips = cached
        downloaded_bytes = 0
        print(f"Reused {len(source_records):,} cached WGS target records from {source_cache}")
    vcf_path = args.out_dir / f"eastern_coyote_wgs_bridge_{args.target_sites}.vcf.gz"
    sample_path = args.out_dir / "eastern_coyote_wgs_bridge_samples.csv"
    retained, merge_skips = _write_bridge_vcf(vcf_path, selected, plan, source_records)
    if retained < args.min_retained_sites:
        vcf_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"only {retained:,} bridge loci passed coordinate and allele checks; "
            f"expected at least {args.min_retained_sites:,}. No bridge VCF was kept."
        )
    _write_sample_sheet(sample_path, args.eastern_sample_sheet, [str(s) for s in eastern.samples])
    manifest = {
        "created_utc": datetime.now(UTC).isoformat(),
        "source_wgs": {"url": args.source_vcf_url, "samples": [row[0] for row in WGS_SAMPLES]},
        "source_radseq": {
            "doi": "10.5061/dryad.7f9q2cd",
            "qc_vcf": str(args.eastern_vcf),
            "n_samples": len(eastern.samples),
        },
        "selection": {
            "requested_target_sites": args.target_sites,
            "autosomal_radseq_sites_available": len(_choose_targets(eastern, len([i for i, c in enumerate(eastern.chrom) if str(c).isdigit() and 1 <= int(str(c)) <= 38]))),
            "n_ranges": estimate.n_ranges,
        },
        "transfer": {
            "estimated_bytes": estimate.total_bytes,
            "downloaded_bytes": downloaded_bytes,
            "hard_ceiling_bytes": args.max_download_bytes,
        },
        "validation": {
            "n_wgs_coordinate_records": len(source_records),
            "n_retained_bridge_sites": retained,
            "source_record_skips": dict(sorted(source_skips.items())),
            "merge_skips": dict(sorted(merge_skips.items())),
            "criterion": "same coordinate and unordered biallelic allele pair; RADseq calls remapped to WGS REF/ALT",
        },
        "outputs": {
            "vcf": str(vcf_path),
            "sample_sheet": str(sample_path),
            "source_record_cache": str(source_cache),
        },
        "limitations": [
            "Cross-platform RADseq/WGS bridge panel, not a whole-genome joint callset.",
            "Loci are selected from an eastern-coyote RADseq panel and are not genome-wide callable sites.",
            "Use for exploratory population comparison; do not treat as standalone wolf-ancestry evidence.",
        ],
    }
    manifest_path = args.out_dir / f"eastern_coyote_wgs_bridge_{args.target_sites}.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {vcf_path}")
    print(f"Wrote {sample_path}")
    print(f"Wrote {manifest_path}")
    print(f"Retained {retained:,} cross-source bridge SNPs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
