"""Assign the dog ancestry of the NYC coydog family to domestic-dog breeds.

The NYC coydog validation scored NY01/NY04/NY05/T211 against only 10 dogs from
3 breeds.  The source WGS VCF (NHGRI 722-genome panel, ~144 breeds) carries
every sample in the same records, so the same ~337 indexed byte ranges already
transferred once can be re-fetched to recover genotype calls for ALL panel
samples at the same bridge loci.

This remains a 252-locus cross-platform diagnostic, not a genomic ancestry
estimate. Pure helpers live in :mod:`canidae.analysis` and
:mod:`canidae.io.indexed_vcf`.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

import validate_nyc_coydog as nyc

from canidae.analysis.breed_panel import (
    INFERRED_VILLAGE_CODES,
    breed_of,
    classify_group,
)
from canidae.analysis.bridge_panel import bridge_sites
from canidae.analysis.genotypes import dosage
from canidae.analysis.reference_mixture import (
    MIN_EFFECTIVE_P,
    SINGLE_BREED_GAP,
    allele_frequencies,
    binom2_logpmf,
    calls_log_likelihood,
    leave_one_out,
    log_likelihood,
)
from canidae.io.indexed_vcf import (
    BGZF_MAX_BLOCK,
    SOURCE_VCF_URL,
    TransferBudget,
    chunks_for,
    decompress_bgzf,
    genotype,
    merge_chunks,
    parse_tabix,
    request,
    source_header,
)


def fetch_all_sample_genotypes(
    sites: dict[tuple[str, int], tuple[str, str]],
    out_dir: Path,
    max_bytes: int,
    known_index: Path | None = None,
) -> tuple[list[str], dict[tuple[str, int], list[str | None]], dict[str, object]]:
    """Fetch every VCF sample at the bridge loci, resuming per byte range."""
    out_dir.mkdir(parents=True, exist_ok=True)
    index_path = out_dir / "722g.990.SNP.INDEL.chrAll.vcf.gz.tbi"
    if not index_path.exists() and known_index is not None and known_index.exists():
        index_path = known_index
    budget = TransferBudget(max_bytes)
    if index_path.exists():
        raw_index = index_path.read_bytes()
    else:
        raw_index = request(f"{SOURCE_VCF_URL}.tbi", budget)
        index_path.write_bytes(raw_index)
    names, indexes = parse_tabix(raw_index)
    by_chrom = dict(zip(names, indexes, strict=True))
    missing_chromosomes = sorted({chrom for chrom, _ in sites if chrom not in by_chrom})
    if missing_chromosomes:
        raise RuntimeError(f"bridge chromosomes absent from WGS index: {missing_chromosomes}")

    chunks: set[tuple[int, int]] = set()
    for chrom, pos in sites:
        chunks.update(chunks_for(by_chrom[chrom], pos))
    merged = merge_chunks(chunks)
    estimated = sum(((end >> 16) + BGZF_MAX_BLOCK) - (begin >> 16) for begin, end in merged)
    if budget.used + estimated > max_bytes:
        raise RuntimeError(
            f"indexed extraction estimate exceeds safety limit: "
            f"{budget.used + estimated:,} > {max_bytes:,} bytes"
        )

    _, columns = source_header(SOURCE_VCF_URL, budget)
    samples = columns[9:]
    if not samples:
        raise RuntimeError("source VCF exposes no genotype columns")

    checkpoint = out_dir / "fetch_checkpoint.jsonl"
    records: dict[tuple[str, int], list[str | None]] = {}
    done: set[int] = set()
    if checkpoint.exists():
        with checkpoint.open(encoding="utf-8") as handle:
            for line in handle:
                entry = json.loads(line)
                done.add(entry["range"])
                for chrom, position, genotypes in entry["records"]:
                    records[(chrom, position)] = genotypes

    for ordinal, (virtual_begin, virtual_end) in enumerate(merged):
        if ordinal in done:
            continue
        start = virtual_begin >> 16
        end = (virtual_end >> 16) + BGZF_MAX_BLOCK - 1
        raw = request(SOURCE_VCF_URL, budget, (start, end))
        text = decompress_bgzf(raw, virtual_begin & 0xFFFF).decode(errors="replace")
        kept: list[list[object]] = []
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[8] == "":
                continue
            try:
                key = (fields[0], int(fields[1]))
            except ValueError:
                continue
            if key not in sites or key in records:
                continue
            ref, alt = sites[key]
            if fields[3].upper() != ref or fields[4].upper() != alt:
                continue
            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                continue
            gt_index = format_fields.index("GT")
            genotypes = [genotype(field, gt_index) for field in fields[9:]]
            records[key] = genotypes
            kept.append([key[0], key[1], genotypes])
        with checkpoint.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps({"range": ordinal, "records": kept}) + "\n")
        if ordinal % 25 == 0 or ordinal == len(merged) - 1:
            print(
                f"WGS ranges {ordinal + 1}/{len(merged)}; "
                f"downloaded {budget.used / 1_000_000:.1f} MB; "
                f"retained {len(records)} loci",
                flush=True,
            )

    metadata = {
        "n_samples": len(samples),
        "n_ranges": len(merged),
        "estimated_bytes": estimated,
        "downloaded_bytes": budget.used,
        "retained_loci": len(records),
    }
    return samples, records, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bridge-vcf",
        type=Path,
        default=Path("data/eastern_coyote_wgs_bridge/eastern_coyote_wgs_bridge_1000.vcf.gz"),
    )
    parser.add_argument(
        "--reference-json",
        type=Path,
        default=Path("data/nyc_coydog_validation/reference_genotypes.json"),
        help="prior 13-sample reference; restricts analysis to its retained loci",
    )
    parser.add_argument(
        "--calls-dir",
        type=Path,
        default=Path("data/nyc_coydog_validation"),
        help="directory with <SAMPLE>_calls.csv from the NYC validation",
    )
    parser.add_argument(
        "--validation-json",
        type=Path,
        default=Path("data/nyc_coydog_validation/validation.json"),
        help="strict manifest carrying per-sample dog_fraction estimates",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("data/nyc_coydog_breeds"))
    parser.add_argument(
        "--known-index",
        type=Path,
        default=Path("data/nyc_coydog_validation/722g.990.SNP.INDEL.chrAll.vcf.gz.tbi"),
    )
    parser.add_argument("--min-depth", type=int, default=8)
    parser.add_argument("--min-group", type=int, default=3)
    parser.add_argument("--max-download-bytes", type=int, default=9_000_000_000)
    parser.add_argument("--top-k", type=int, default=5, help="breeds to report per NYC sample")
    args = parser.parse_args()
    if not 1 <= args.max_download_bytes < 10_000_000_000:
        parser.error("--max-download-bytes must be positive and below 10 GB")
    if args.min_group < 2:
        parser.error("--min-group must be at least 2")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sites = bridge_sites(args.bridge_vcf)
    saved = json.loads(args.reference_json.read_text(encoding="utf-8"))
    reference_keys = sorted(
        ((key.rsplit(":", 1)[0], int(key.rsplit(":", 1)[1])) for key in saved),
        key=lambda item: (int(item[0].removeprefix("chr")), item[1]),
    )
    print(f"Using {len(reference_keys)} previously retained reference loci", flush=True)

    samples, records, fetch_metadata = fetch_all_sample_genotypes(
        sites, args.out_dir, args.max_download_bytes, args.known_index
    )
    records = {key: records[key] for key in reference_keys}
    (args.out_dir / "all_sample_genotypes.json").write_text(
        json.dumps(
            {
                "samples": samples,
                "loci": {
                    f"{chrom}:{pos}": records[(chrom, pos)] for (chrom, pos) in reference_keys
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Fetched {len(samples)} samples at {len(records)} loci", flush=True)

    groups: dict[str, list[int]] = defaultdict(list)
    for index, sample in enumerate(samples):
        groups[breed_of(sample)].append(index)
    categories = {group: classify_group(group) for group in groups}
    breed_groups = {
        group: indices
        for group, indices in groups.items()
        if categories[group] == "breed" and len(indices) >= args.min_group
    }
    village_indices = [
        index
        for group, indices in groups.items()
        if categories[group] == "village"
        for index in indices
    ]
    candidates: dict[str, list[int]] = dict(breed_groups)
    if len(village_indices) >= args.min_group:
        candidates["VillageDog(all regions)"] = village_indices
    for code, region in INFERRED_VILLAGE_CODES.items():
        if len(groups.get(code, [])) >= args.min_group:
            candidates[f"VillageDog({region}, inferred)"] = groups[code]
    wild_groups = {
        group: indices for group, indices in groups.items() if categories[group] == "wild"
    }
    pooled_dogs = [
        index
        for group, indices in groups.items()
        if categories[group] != "wild"
        for index in indices
    ]
    print(
        f"{len(candidates)} candidate groups (>= {args.min_group} members, "
        f"{len(pooled_dogs)} dogs pooled); wild outgroups: "
        f"{sorted(wild_groups) or 'none'}",
        flush=True,
    )
    coyote_indices = next(
        (indices for name, indices in wild_groups.items() if "coyote" in name.lower()), []
    )
    if not coyote_indices:
        raise RuntimeError("no coyote group available for the mixture model")
    coyote_panel = allele_frequencies(reference_keys, records, coyote_indices)
    pooled_panel = allele_frequencies(reference_keys, records, pooled_dogs)
    panels = {
        group: allele_frequencies(reference_keys, records, indices)
        for group, indices in candidates.items()
    }

    loo_rows, loo_summary = leave_one_out(
        reference_keys, records, samples, candidates, args.min_group
    )
    with (args.out_dir / "leave_one_out.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(loo_rows[0]))
        writer.writeheader()
        writer.writerows(loo_rows)
    print(f"Leave-one-out top-1 accuracy: {loo_summary}", flush=True)

    validation = json.loads(args.validation_json.read_text(encoding="utf-8"))
    dog_fractions = {str(row["sample_id"]): row["dog_fraction"] for row in validation["results"]}
    score_rows: list[dict[str, object]] = []
    per_sample_top: dict[str, object] = {}
    for sample_id in nyc.NYC_RUNS:
        calls = nyc.recalculate_calls_from_csv(
            args.calls_dir / f"{sample_id}_calls.csv",
            {key: sites[key] for key in reference_keys},
            args.min_depth,
        )
        dog_fraction = dog_fractions.get(sample_id)
        if dog_fraction is None:
            print(f"Skipping {sample_id}: no dog_fraction estimate", flush=True)
            continue
        entries: list[dict[str, object]] = []
        for breed, panel in panels.items():
            ll, used = calls_log_likelihood(
                calls, reference_keys, panel, coyote_panel, float(dog_fraction)
            )
            entries.append({"breed": breed, "log_likelihood": round(ll, 2), "loci_used": used})
        any_ll, any_used = calls_log_likelihood(
            calls, reference_keys, pooled_panel, coyote_panel, float(dog_fraction)
        )
        entries.sort(key=lambda entry: entry["log_likelihood"], reverse=True)
        for rank, entry in enumerate(entries[: args.top_k], start=1):
            score_rows.append(
                {
                    "sample_id": sample_id,
                    "rank": rank,
                    "breed": entry["breed"],
                    "log_likelihood": entry["log_likelihood"],
                    "gap_to_next": (
                        round(
                            float(entries[rank - 1]["log_likelihood"])
                            - float(entries[rank]["log_likelihood"]),
                            2,
                        )
                        if rank < len(entries)
                        else None
                    ),
                    "loci_used": entry["loci_used"],
                }
            )
        best = entries[0]
        per_sample_top[sample_id] = {
            "dog_fraction": dog_fraction,
            "best_breed": best["breed"],
            "best_log_likelihood": best["log_likelihood"],
            "any_dog_log_likelihood": round(any_ll, 2),
            "single_breed_gap": round(float(best["log_likelihood"]) - any_ll, 2),
            "single_breed_supported": bool(
                float(best["log_likelihood"]) - any_ll > SINGLE_BREED_GAP
            ),
            "loci_used": any_used,
        }
        print(
            f"{sample_id}: best {best['breed']} "
            f"(LL {best['log_likelihood']}, any-dog {any_ll:.2f})",
            flush=True,
        )

    with (args.out_dir / "breed_scores.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(score_rows[0]))
        writer.writeheader()
        writer.writerows(score_rows)

    (args.out_dir / "breed_assignment.json").write_text(
        json.dumps(
            {
                "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
                "sources": {
                    "wgs_bioproject": "PRJNA448733",
                    "n_wgs_samples": len(samples),
                    "n_candidate_groups": len(candidates),
                    "candidate_sizes": {
                        group: len(indices)
                        for group, indices in sorted(
                            candidates.items(), key=lambda item: -len(item[1])
                        )
                    },
                    "group_categories": {group: categories[group] for group in sorted(groups)},
                    "inferred_labels": {
                        code: f"VillageDog({region}, inferred)"
                        for code, region in INFERRED_VILLAGE_CODES.items()
                        if len(groups.get(code, [])) >= args.min_group
                    },
                    "wild_outgroups": {
                        group: len(indices) for group, indices in sorted(wild_groups.items())
                    },
                    "loci": len(reference_keys),
                },
                "transfer": fetch_metadata,
                "calibration": loo_summary,
                "results": per_sample_top,
                "limitations": [
                    "252 cross-platform bridge loci; breed resolution is far below a "
                    "genotype-array or WGS standard.",
                    "The model scores which breed panel best explains the dog-derived "
                    "alleles; it cannot phase or directly observe the parent.",
                    "A mixed-breed or unpanelled dog parent should appear as a small "
                    "single-breed gap versus the pooled any-dog panel.",
                    "F1/backcross mixture likelihood assumes independent allele draws, an "
                    "approximation at this marker count.",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {args.out_dir / 'breed_assignment.json'}", flush=True)


__all__ = [
    "INFERRED_VILLAGE_CODES",
    "MIN_EFFECTIVE_P",
    "SINGLE_BREED_GAP",
    "allele_frequencies",
    "binom2_logpmf",
    "breed_of",
    "calls_log_likelihood",
    "classify_group",
    "dosage",
    "fetch_all_sample_genotypes",
    "leave_one_out",
    "log_likelihood",
    "main",
]


if __name__ == "__main__":
    main()
