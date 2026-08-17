"""Cross-platform WGS/RADseq bridge-panel I/O and known data caveats."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from canidae.analysis.genotypes import allele_count, parse_gt_field

# NOTE: the 2026-07-13 bridge VCF's six WGS sample columns carry mislabeled
# genotypes (verified against live source records: e.g. its "Coyote01" column
# matches panel sample 140447_S11).  Its 36 RADseq columns are unaffected.
# Analyses must never use the bridge VCF's WGS columns; WGS genotypes come from
# independently fetched panel records.
BRIDGE_WGS_COLUMNS_BUG = True
BRIDGE_WGS_IDS = frozenset(
    {
        "Coyote01",
        "Coyote02",
        "AlaskanWolf",
        "AlgonquinWolf13467",
        "AlgonquinWolf13470",
        "GoldenJackal01",
    }
)


def bridge_sites(path: Path) -> dict[tuple[str, int], tuple[str, str]]:
    """Load chrom/pos -> (REF, ALT) from a bridge VCF; requires >=100 loci."""
    sites: dict[tuple[str, int], tuple[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 5:
                continue
            sites[(fields[0], int(fields[1]))] = (fields[3].upper(), fields[4].upper())
    if len(sites) < 100:
        raise RuntimeError(f"bridge panel is unexpectedly small ({len(sites)} loci)")
    return sites


def load_bridge_genotypes(
    path: Path,
) -> tuple[list[str], dict[tuple[str, int], dict[str, str | None]]]:
    """Read the bridge VCF: per-locus sample -> genotype (WGS-oriented)."""
    samples: list[str] = []
    calls: dict[tuple[str, int], dict[str, str | None]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
                continue
            format_fields = fields[8].split(":")
            gt_index = format_fields.index("GT")
            key = (f"chr{fields[0].removeprefix('chr')}", int(fields[1]))
            per_sample: dict[str, str | None] = {}
            for sample, field in zip(samples, fields[9:], strict=True):
                per_sample[sample] = parse_gt_field(field, gt_index)
            calls[key] = per_sample
    return samples, calls


def allele_count_vector(
    keys: list[tuple[str, int]],
    calls: dict[tuple[str, int], dict[str, str | None]],
    sample: str,
) -> np.ndarray:
    """One sample's ALT allele counts over *keys*; NaN where the locus is uncalled."""
    out = np.full(len(keys), np.nan)
    for position, key in enumerate(keys):
        genotype = calls[key].get(sample)
        if genotype is not None:
            out[position] = allele_count(genotype)
    return out
