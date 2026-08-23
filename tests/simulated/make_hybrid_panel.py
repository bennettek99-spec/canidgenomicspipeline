"""Synthetic hybrid-canid bridge panel with known mixture fractions.

Builds, from a fixed seed, the exact on-disk inputs the hybrid diagnostic stages
consume — a bridge VCF, a WGS panel genotype JSON, a two-source reference
genotype JSON, and per-query allele-count tables — so estimators can be checked
against ground truth without any network access or real study data.

The locus layout is deliberately blocked so each test assertion has a mechanism:

===========  =====  ===================================================
Locus block  Count  Role
===========  =====  ===================================================
0-119        120    Coyote vs dog diagnostic (wolf tracks coyote here)
120-179       60    Wolf vs coyote diagnostic (dogs track the background)
180-259       80    Breed/village private, 20 loci per candidate group
260-399      140    Neutral background shared by every group
===========  =====  ===================================================

Queries are drawn from the mixture model the stages fit, ``p_eff = (1 - f_wolf -
f_dog) * p_coyote + f_wolf * p_wolf + f_dog * p_dog_source``, so the recovered
fractions are unbiased estimates of :attr:`HybridPanel.truth`.
"""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

CHROM = "chr1"
N_LOCI = 400
DOG_DIAGNOSTIC = slice(0, 120)
WOLF_DIAGNOSTIC = slice(120, 180)
PRIVATE_START = 180
PRIVATE_PER_GROUP = 20
NEUTRAL_START = 260
# Loci retained in the two-source reference JSON; the rest exercise the
# breed_assign locus filter (which mirrors the real NYC retained-site subset).
RETAINED_LOCI = 300

COYOTE_SAMPLES = [f"Coyote{i:02d}" for i in range(1, 7)]
WOLF_SAMPLES = [f"AlaskanWolf{i:02d}" for i in range(1, 6)]
BREED_SAMPLES = {
    "Beagle": [f"Beagle{i:02d}" for i in range(1, 6)],
    "GermanShepherd": [f"GermanShepherd{i:02d}" for i in range(1, 6)],
    "LabradorRetriever": [f"LabradorRetriever{i:02d}" for i in range(1, 6)],
}
VILLAGE_SAMPLES = [f"VillDog_Peru{i:02d}" for i in range(1, 6)]
# Groups that get a private-allele block, in the order those blocks are laid out.
PRIVATE_GROUPS = ("Beagle", "GermanShepherd", "LabradorRetriever", "VillDog_Peru")

# One bridge column named like a real WGS panel sample. The 2026-07-13 bridge
# VCF mislabels exactly these columns, so stages must be able to drop it.
BUG_COLUMN = "GoldenJackal01"

CALL_DEPTH = 20
LOW_DEPTH = 4
# Every 20th locus is called too shallow to use, so the depth filter and the
# missing-genotype paths are exercised without starving the diagnostic blocks.
LOW_DEPTH_EVERY = 20


@dataclass(frozen=True)
class Query:
    """One simulated query: its true mixture and the dog panel it was drawn from."""

    sample_id: str
    f_wolf: float
    f_dog: float
    dog_source: str  # a PRIVATE_GROUPS name, or "pooled" for a mean-of-dogs source
    region: str  # "eastern" or "western", for the multiway group summary


PEDIGREE_QUERIES = (
    Query("NYFIXTURE_F1", 0.0, 0.50, "Beagle", "eastern"),
    Query("NYFIXTURE_PARENT", 0.0, 0.05, "Beagle", "eastern"),
    Query("NYFIXTURE_OFF1", 0.0, 0.25, "Beagle", "eastern"),
    Query("NYFIXTURE_OFF2", 0.0, 0.28, "Beagle", "eastern"),
)
EASTERN_QUERIES = tuple(Query(f"EC{i:02d}", 0.25, 0.10, "pooled", "eastern") for i in range(1, 11))
WESTERN_QUERIES = tuple(Query(f"WC{i:02d}", 0.0, 0.0, "pooled", "western") for i in range(1, 7))
BUG_QUERY = Query(BUG_COLUMN, 0.0, 0.0, "pooled", "eastern")
ALL_QUERIES = PEDIGREE_QUERIES + EASTERN_QUERIES + WESTERN_QUERIES + (BUG_QUERY,)


@dataclass(frozen=True)
class HybridPanel:
    """Paths to a materialized synthetic panel plus its ground truth."""

    root: Path
    bridge_vcf: Path
    wgs_genotypes: Path
    reference_genotypes: Path
    calls_dir: Path
    dog_fractions: Path
    wgs_samples: list[str]
    reference_samples: list[str]
    coyote_samples: list[str]
    truth: dict[str, Query]

    def query_ids(self, region: str | None = None) -> list[str]:
        return [
            q.sample_id
            for q in self.truth.values()
            if q.sample_id != BUG_COLUMN and (region is None or q.region == region)
        ]

    def two_source_query_ids(self) -> list[str]:
        """Queries with no wolf component, so a coyote-vs-dog fit is well specified."""
        return [s for s in self.query_ids() if self.truth[s].f_wolf == 0.0]

    def wolf_carrying_query_ids(self) -> list[str]:
        """Queries a two-source model cannot describe; their wolf loads onto dog."""
        return [s for s in self.query_ids() if self.truth[s].f_wolf > 0.0]


def group_frequencies(rng: np.random.Generator) -> dict[str, np.ndarray]:
    """True per-group ALT frequencies over the blocked locus layout."""
    background = rng.uniform(0.15, 0.85, N_LOCI)
    dog_groups = [*BREED_SAMPLES, "VillDog_Peru"]
    freqs = {group: background.copy() for group in ("Coyote", "AlaskanWolf", *dog_groups)}

    # Dogs carry the derived allele; wolf sits close to coyote so this block
    # informs the dog axis only.
    freqs["Coyote"][DOG_DIAGNOSTIC] = 0.03
    freqs["AlaskanWolf"][DOG_DIAGNOSTIC] = 0.05
    for group in dog_groups:
        freqs[group][DOG_DIAGNOSTIC] = 0.97

    # Wolf-only derived allele; dogs stay on the shared background.
    freqs["Coyote"][WOLF_DIAGNOSTIC] = 0.03
    freqs["AlaskanWolf"][WOLF_DIAGNOSTIC] = 0.97

    # Private blocks make each candidate panel separable for breed assignment.
    for offset, owner in enumerate(PRIVATE_GROUPS):
        start = PRIVATE_START + offset * PRIVATE_PER_GROUP
        block = slice(start, start + PRIVATE_PER_GROUP)
        for group in dog_groups:
            freqs[group][block] = 0.95 if group == owner else 0.05
    return freqs


def _genotype_string(alt_count: int) -> str:
    """VCF GT for an ALT allele count in {0, 1, 2}."""
    return ("0/0", "0/1", "1/1")[alt_count]


def _panel_membership() -> list[tuple[str, str]]:
    """(sample_id, group) for every WGS panel member, in file order."""
    members = [(s, "Coyote") for s in COYOTE_SAMPLES]
    members += [(s, "AlaskanWolf") for s in WOLF_SAMPLES]
    for breed, samples in BREED_SAMPLES.items():
        members += [(s, breed) for s in samples]
    members += [(s, "VillDog_Peru") for s in VILLAGE_SAMPLES]
    return members


def _query_frequencies(query: Query, freqs: dict[str, np.ndarray]) -> np.ndarray:
    if query.dog_source == "pooled":
        p_dog = np.mean([freqs[g] for g in (*BREED_SAMPLES, "VillDog_Peru")], axis=0)
    else:
        p_dog = freqs[query.dog_source]
    f_coyote = 1.0 - query.f_wolf - query.f_dog
    return f_coyote * freqs["Coyote"] + query.f_wolf * freqs["AlaskanWolf"] + query.f_dog * p_dog


def build_hybrid_panel(out_dir: Path, seed: int = 20260817) -> HybridPanel:
    """Materialize a synthetic bridge panel under *out_dir* and return its paths."""
    out_dir = Path(out_dir)
    calls_dir = out_dir / "calls"
    calls_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    freqs = group_frequencies(rng)
    positions = [(i + 1) * 1000 for i in range(N_LOCI)]
    keys = [f"{CHROM}:{pos}" for pos in positions]

    members = _panel_membership()
    panel_samples = [sample for sample, _ in members]
    panel_counts = np.vstack(
        [rng.binomial(2, freqs[group]) for _, group in members]
    )  # ALT allele counts, (n_panel, n_loci)

    wgs_payload = {
        "samples": panel_samples,
        "loci": {
            key: [_genotype_string(int(c)) for c in panel_counts[:, i]]
            for i, key in enumerate(keys)
        },
    }
    wgs_path = out_dir / "all_sample_genotypes.json"
    wgs_path.write_text(json.dumps(wgs_payload), encoding="utf-8")

    # Two-source reference: coyotes + dogs only. The stage treats every
    # non-coyote entry as dog, so wolves must not appear here.
    reference_samples = [sample for sample, group in members if group not in ("AlaskanWolf",)]
    reference_rows = [panel_samples.index(s) for s in reference_samples]
    reference_payload = {
        key: [_genotype_string(int(c)) for c in panel_counts[reference_rows, i]]
        for i, key in enumerate(keys[:RETAINED_LOCI])
    }
    reference_path = out_dir / "reference_genotypes.json"
    reference_path.write_text(json.dumps(reference_payload), encoding="utf-8")

    # Query genotypes, drawn from the same mixture model the stages fit.
    query_counts: dict[str, Any] = {
        query.sample_id: rng.binomial(2, _query_frequencies(query, freqs)) for query in ALL_QUERIES
    }
    depths = np.full(N_LOCI, CALL_DEPTH)
    depths[::LOW_DEPTH_EVERY] = LOW_DEPTH

    for query in ALL_QUERIES:
        counts = query_counts[query.sample_id]
        path = calls_dir / f"{query.sample_id}_calls.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["chrom", "position", "ref", "alt", "gt", "depth", "ref_count", "alt_count"]
            )
            for i, pos in enumerate(positions):
                depth = int(depths[i])
                alt = round(counts[i] / 2.0 * depth)
                writer.writerow(
                    [
                        CHROM,
                        pos,
                        "A",
                        "G",
                        _genotype_string(int(counts[i])),
                        depth,
                        depth - alt,
                        alt,
                    ]
                )

    bridge_path = out_dir / "bridge.vcf.gz"
    bridge_columns = [q.sample_id for q in ALL_QUERIES]
    with gzip.open(bridge_path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##contig=<ID={CHROM},length={N_LOCI * 1000 + 1000}>\n")
        handle.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(bridge_columns)
            + "\n"
        )
        for i, pos in enumerate(positions):
            calls = "\t".join(_genotype_string(int(query_counts[s][i])) for s in bridge_columns)
            handle.write(f"{CHROM}\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT\t{calls}\n")

    # A plain sample -> fraction map, the simplest breed_assign input shape.
    fractions_path = out_dir / "dog_fractions.json"
    fractions_path.write_text(
        json.dumps({q.sample_id: q.f_dog for q in PEDIGREE_QUERIES}, indent=2),
        encoding="utf-8",
    )

    return HybridPanel(
        root=out_dir,
        bridge_vcf=bridge_path,
        wgs_genotypes=wgs_path,
        reference_genotypes=reference_path,
        calls_dir=calls_dir,
        dog_fractions=fractions_path,
        wgs_samples=panel_samples,
        reference_samples=reference_samples,
        coyote_samples=list(COYOTE_SAMPLES),
        truth={q.sample_id: q for q in ALL_QUERIES},
    )
