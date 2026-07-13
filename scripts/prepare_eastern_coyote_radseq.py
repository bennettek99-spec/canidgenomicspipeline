"""Create a small, auditable VCF from the Heppenheimer et al. RADseq PLINK files.

The Dryad deposit provides a PLINK PED/MAP genotype matrix rather than a VCF.  This
utility selects the specified samples and materializes its biallelic SNPs as a VCF
that CANIS can use on a laptop.  The result is intentionally marked as *derived*:
the PED representation does not contain an assembly reference allele, so REF/ALT
are deterministic cohort allele labels and must not be used for reference-oriented
annotation or locus-level allele interpretation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable


SOURCE_DOI = "10.5061/dryad.7f9q2cd"
SOURCE_URL = "https://datadryad.org/dataset/doi:10.5061/dryad.7f9q2cd"
SOURCE_CITATION = "Heppenheimer et al. (2018), Ecology and Evolution"
VALID_BASES = frozenset({"A", "C", "G", "T"})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--map",
        type=Path,
        default=Path(
            "data/eastern_coyote_radseq_source/dryad_2018/"
            "Heppenheimer_et_al_Ecol_Evol_2018.map"
        ),
        help="Source PLINK MAP file.",
    )
    parser.add_argument(
        "--ped",
        type=Path,
        default=Path(
            "data/eastern_coyote_radseq_source/dryad_2018/"
            "Heppenheimer_et_al_Ecol_Evol_2018.ped"
        ),
        help="Source PLINK PED file.",
    )
    parser.add_argument(
        "--sample-sheet",
        type=Path,
        default=Path("configs/examples/eastern_coyote_radseq_samples.csv"),
        help="CANIS sample sheet. Its sample_id values select PED records.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/eastern_coyote_radseq"),
        help="Directory for the derived VCF and manifest.",
    )
    parser.add_argument(
        "--min-site-call-rate",
        type=float,
        default=0.95,
        help="Minimum complete diploid call rate across the selected cohort.",
    )
    parser.add_argument(
        "--min-retained-sites",
        type=int,
        default=1000,
        help="Fail if fewer sites survive conversion, guarding against a bad selection.",
    )
    return parser.parse_args()


def read_sample_ids(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "sample_id" not in rows[0]:
        raise ValueError(f"sample sheet must contain at least one sample_id: {path}")
    sample_ids = [row["sample_id"].strip() for row in rows]
    if any(not sample_id for sample_id in sample_ids):
        raise ValueError(f"sample sheet contains an empty sample_id: {path}")
    duplicates = sorted({sample_id for sample_id in sample_ids if sample_ids.count(sample_id) > 1})
    if duplicates:
        raise ValueError(f"sample sheet has duplicate sample IDs: {duplicates}")
    return sample_ids


def read_map(path: Path) -> list[tuple[str, str, int]]:
    markers: list[tuple[str, str, int]] = []
    seen_coordinates: set[tuple[str, int]] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.split()
            if len(fields) != 4:
                raise ValueError(f"{path}:{line_number}: expected 4 MAP columns, got {len(fields)}")
            chrom, marker_id, _genetic_distance, raw_position = fields
            try:
                position = int(raw_position)
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid position {raw_position!r}") from exc
            if position < 1:
                raise ValueError(f"{path}:{line_number}: position must be positive, got {position}")
            coordinate = (chrom, position)
            if coordinate in seen_coordinates:
                raise ValueError(
                    f"{path}:{line_number}: duplicate chromosome/position {chrom}:{position}; "
                    "cannot make a standards-compliant VCF"
                )
            seen_coordinates.add(coordinate)
            markers.append((chrom, marker_id, position))
    if not markers:
        raise ValueError(f"MAP file is empty: {path}")
    return markers


def read_selected_ped(path: Path, sample_ids: Iterable[str], n_markers: int) -> dict[str, list[str]]:
    wanted = set(sample_ids)
    selected: dict[str, list[str]] = {}
    expected_fields = 6 + 2 * n_markers
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed PED row")
            sample_id = fields[1]
            if sample_id not in wanted:
                continue
            if len(fields) != expected_fields:
                raise ValueError(
                    f"{path}:{line_number}: expected {expected_fields} fields for {n_markers} markers, "
                    f"got {len(fields)}"
                )
            if sample_id in selected:
                raise ValueError(f"{path}:{line_number}: duplicate PED sample {sample_id}")
            selected[sample_id] = fields[6:]
    missing = [sample_id for sample_id in sample_ids if sample_id not in selected]
    if missing:
        raise ValueError(f"selected samples are absent from PED: {', '.join(missing)}")
    return selected


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def format_genotype(alleles: tuple[str, str], ref: str, alt: str) -> str:
    left, right = alleles
    if left == "0" or right == "0":
        return "./."
    if left not in VALID_BASES or right not in VALID_BASES:
        return "./."
    allele_to_index = {ref: "0", alt: "1"}
    if left not in allele_to_index or right not in allele_to_index:
        return "./."
    return f"{allele_to_index[left]}/{allele_to_index[right]}"


def write_vcf(
    destination: Path,
    markers: list[tuple[str, str, int]],
    selected: dict[str, list[str]],
    sample_ids: list[str],
    min_site_call_rate: float,
) -> tuple[int, Counter[str]]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    counters: Counter[str] = Counter()
    retained = 0
    with gzip.open(destination, "wt", encoding="utf-8", newline="") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=Heppenheimer_et_al_Ecol_Evol_2018_Dryad_PLINK_PED_MAP\n")
        handle.write(f"##source_doi={SOURCE_DOI}\n")
        handle.write("##CANIS_origin=locally derived from PLINK PED/MAP; not a primary called VCF\n")
        handle.write(
            "##CANIS_allele_note=REF and ALT are sorted cohort allele labels, not "
            "reference-genome-oriented alleles\n"
        )
        handle.write(
            "##CANIS_coordinate_note=coordinates are copied from the deposited MAP file; "
            "verify assembly before locus-level interpretation\n"
        )
        handle.write(
            "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source data DOI\">\n"
        )
        handle.write(
            "##INFO=<ID=CALLRATE,Number=1,Type=Float,Description=\"Complete diploid "
            "call rate in the selected cohort\">\n"
        )
        handle.write(
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Unphased genotype\">\n"
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        handle.write("\t".join(sample_ids))
        handle.write("\n")

        for marker_index, (chrom, marker_id, position) in enumerate(markers):
            genotype_alleles = [
                (selected[sample_id][2 * marker_index], selected[sample_id][2 * marker_index + 1])
                for sample_id in sample_ids
            ]
            called = [
                pair for pair in genotype_alleles
                if pair[0] in VALID_BASES and pair[1] in VALID_BASES
            ]
            call_rate = len(called) / len(sample_ids)
            if call_rate < min_site_call_rate:
                counters["low_call_rate"] += 1
                continue

            observed = sorted(
                {allele for pair in genotype_alleles for allele in pair if allele in VALID_BASES}
            )
            if len(observed) < 2:
                counters["monomorphic"] += 1
                continue
            if len(observed) > 2:
                counters["multiallelic"] += 1
                continue
            ref, alt = observed
            formatted = [format_genotype(pair, ref, alt) for pair in genotype_alleles]
            if any(value == "./." and pair[0] != "0" and pair[1] != "0" for value, pair in zip(formatted, genotype_alleles, strict=True)):
                counters["invalid_alleles"] += 1
                continue

            handle.write(
                f"{chrom}\t{position}\t{marker_id}\t{ref}\t{alt}\t.\tPASS\t"
                f"SOURCE={SOURCE_DOI};CALLRATE={call_rate:.4f}\tGT\t"
            )
            handle.write("\t".join(formatted))
            handle.write("\n")
            retained += 1
    return retained, counters


def main() -> int:
    args = parse_args()
    if not 0 < args.min_site_call_rate <= 1:
        raise ValueError("--min-site-call-rate must be in (0, 1]")
    if args.min_retained_sites < 1:
        raise ValueError("--min-retained-sites must be positive")

    sample_ids = read_sample_ids(args.sample_sheet)
    markers = read_map(args.map)
    selected = read_selected_ped(args.ped, sample_ids, len(markers))

    vcf_path = args.out_dir / "eastern_coyote_ontario_arizona_radseq.vcf.gz"
    retained, counters = write_vcf(
        vcf_path, markers, selected, sample_ids, args.min_site_call_rate
    )
    if retained < args.min_retained_sites:
        vcf_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"only {retained} sites survived conversion; expected at least "
            f"{args.min_retained_sites}. No output VCF was kept."
        )

    manifest = {
        "created_utc": datetime.now(UTC).isoformat(),
        "source": {
            "citation": SOURCE_CITATION,
            "doi": SOURCE_DOI,
            "url": SOURCE_URL,
            "format": "PLINK PED/MAP genotype calls",
            "map": str(args.map),
            "ped": str(args.ped),
            "map_sha256": sha256(args.map),
            "ped_sha256": sha256(args.ped),
        },
        "derived_callset": {
            "path": str(vcf_path),
            "format": "VCFv4.2 gzipped",
            "origin": "locally derived from the source PED/MAP files",
            "samples": sample_ids,
            "n_samples": len(sample_ids),
            "n_source_markers": len(markers),
            "n_retained_biallelic_snps": retained,
            "min_site_call_rate": args.min_site_call_rate,
            "skipped": dict(sorted(counters.items())),
            "allele_orientation": (
                "REF/ALT are deterministic sorted cohort allele labels; they are not "
                "reference-genome oriented."
            ),
            "appropriate_use": (
                "Cohort-level PCA, FST, diversity, distance, and NMF exploratory analysis. "
                "Not reference-allele annotation, whole-genome inference, or standalone "
                "wolf-ancestry estimation."
            ),
        },
    }
    manifest_path = args.out_dir / "eastern_coyote_ontario_arizona_radseq.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {vcf_path}")
    print(f"Wrote {manifest_path}")
    print(
        f"Selected {len(sample_ids)} samples; retained {retained:,} biallelic SNPs "
        f"from {len(markers):,} source loci."
    )
    if counters:
        print("Skipped loci: " + ", ".join(f"{name}={count:,}" for name, count in sorted(counters.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
