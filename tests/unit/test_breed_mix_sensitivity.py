"""Sensitivity of bridge-locus breed assignment to a true 50:50 mix.

The real panel carries a declared Kerry Blue Terrier x Beagle cross; the
synthetic panel lets us construct an F1 of two breeds with known 20-locus
private blocks, so the assertion is mechanistic: a real mix must NOT be
confidently assigned a single breed, while a held-out pure breed still is.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from canidae.analysis.breed_panel import (
    INFERRED_VILLAGE_CODES,
    breed_of,
    classify_group,
)
from canidae.analysis.genotypes import allele_count
from canidae.analysis.reference_mixture import (
    SINGLE_BREED_GAP,
    allele_frequencies,
    score_breed_candidates,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_hybrid_panel import build_hybrid_panel


def _genotype_string(alt_count: int) -> str:
    return ("0/0", "0/1", "1/1")[int(alt_count)]


def _panels(
    keys: list[tuple[str, int]],
    records: dict[tuple[str, int], list[str | None]],
    samples: list[str],
    min_group: int = 3,
) -> tuple[dict[str, list[int]], dict[str, np.ndarray], list[int]]:
    groups: dict[str, list[int]] = defaultdict(list)
    for index, sample in enumerate(samples):
        groups[breed_of(sample)].append(index)
    categories = {group: classify_group(group) for group in groups}
    candidates = {
        group: indices
        for group, indices in groups.items()
        if categories[group] == "breed" and len(indices) >= min_group
    }
    village = [
        index
        for group, indices in groups.items()
        if categories[group] == "village"
        for index in indices
    ]
    if len(village) >= min_group:
        candidates["VillageDog(all regions)"] = village
    for code, region in INFERRED_VILLAGE_CODES.items():
        if len(groups.get(code, [])) >= min_group:
            candidates[f"VillageDog({region}, inferred)"] = groups[code]
    pooled = [
        index
        for group, indices in groups.items()
        if categories[group] != "wild"
        for index in indices
    ]
    panels = {
        group: allele_frequencies(keys, records, indices)
        for group, indices in candidates.items()
    }
    return candidates, panels, pooled


def _synthesize_f1(
    keys: list[tuple[str, int]],
    records: dict[tuple[str, int], list[str | None]],
    parent_a: int,
    parent_b: int,
    rng: np.random.Generator,
) -> int:
    """Append a 50:50 F1 of two panel individuals to ``records`` and return its index."""
    mix_index = len(next(iter(records.values())))
    for key in keys:
        a = allele_count(records[key][parent_a])  # type: ignore[index]
        b = allele_count(records[key][parent_b])  # type: ignore[index]
        alt = rng.binomial(1, a / 2.0) + rng.binomial(1, b / 2.0)
        records[key].append(_genotype_string(int(alt)))
    return mix_index


def test_fifty_fifty_mix_is_rejected_as_a_single_breed(tmp_path: Path) -> None:
    panel = build_hybrid_panel(tmp_path / "panel")
    samples = panel.wgs_samples
    records = {
        (k.rsplit(":", 1)[0], int(k.rsplit(":", 1)[1])): list(v)
        for k, v in __import__("json").loads(
            panel.wgs_genotypes.read_text(encoding="utf-8")
        )["loci"].items()
    }
    keys = sorted(records, key=lambda item: (int(item[0].removeprefix("chr")), item[1]))

    groups: dict[str, list[int]] = defaultdict(list)
    for index, sample in enumerate(samples):
        groups[breed_of(sample)].append(index)
    beagle = groups["Beagle"][0]
    shepherd = groups["GermanShepherd"][0]

    _, panels, _ = _panels(keys, records, samples)
    pooled = allele_frequencies(
        keys, records,
        [i for i, s in enumerate(samples) if classify_group(breed_of(s)) != "wild"],
    )

    mix_index = _synthesize_f1(keys, records, beagle, shepherd, np.random.default_rng(1))
    result = score_breed_candidates(keys, records, mix_index, [*samples, "MIX50"], panels, pooled)

    assert result["single_breed_supported"] is False
    assert result["single_breed_gap"] < SINGLE_BREED_GAP
    # The two source breeds should not individually stand out above any-dog.
    assert result["best_breed"] != "Beagle" or result["single_breed_gap"] < 5.0


def test_held_out_pure_breed_is_confidently_assigned(tmp_path: Path) -> None:
    panel = build_hybrid_panel(tmp_path / "panel")
    samples = panel.wgs_samples
    records = {
        (k.rsplit(":", 1)[0], int(k.rsplit(":", 1)[1])): list(v)
        for k, v in __import__("json").loads(
            panel.wgs_genotypes.read_text(encoding="utf-8")
        )["loci"].items()
    }
    keys = sorted(records, key=lambda item: (int(item[0].removeprefix("chr")), item[1]))

    groups: dict[str, list[int]] = defaultdict(list)
    for index, sample in enumerate(samples):
        groups[breed_of(sample)].append(index)
    beagle = groups["Beagle"]

    _candidates, panels, pooled = _panels(keys, records, samples)
    held_out = beagle[0]
    # Leave the held-out individual out of its own panel and the pooled baseline.
    panels["Beagle"] = allele_frequencies(keys, records, [i for i in beagle if i != held_out])
    any_panel = allele_frequencies(
        keys, records, [i for i in pooled if i != held_out]
    )
    result = score_breed_candidates(keys, records, held_out, samples, panels, any_panel)

    assert result["best_breed"] == "Beagle"
    assert result["single_breed_supported"] is True
    assert result["single_breed_gap"] > SINGLE_BREED_GAP
