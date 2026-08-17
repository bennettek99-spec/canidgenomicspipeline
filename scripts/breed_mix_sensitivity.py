"""Mixed-breed sensitivity check for the bridge-locus breed assignment.

The panel includes declared mixed/unknown dogs (notably the 50:50
``MIX_KerryBlueTerrier_Beagle01`` cross). Scoring them with the same
candidate-breed likelihood model the NYC coydogs use verifies the estimator's
sensitivity: a true mix must NOT be confidently assigned a single breed, while a
(possibly mislabeled) pure-breed sample still can be.

This reads the already-fetched panel JSON (see scripts/breed_assign_nyc_coydog.py)
and writes a CSV + JSON manifest. Pure helpers live in :mod:`canidae.analysis`.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

from canidae.analysis.breed_panel import (
    INFERRED_VILLAGE_CODES,
    breed_of,
    classify_group,
)
from canidae.analysis.reference_mixture import (
    SINGLE_BREED_GAP,
    allele_frequencies,
    score_breed_candidates,
)


def _load_panel(path: Path) -> tuple[list[str], dict[tuple[str, int], list[str | None]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    samples = list(payload["samples"])
    records = {
        (key.rsplit(":", 1)[0], int(key.rsplit(":", 1)[1])): list(values)
        for key, values in payload["loci"].items()
    }
    return samples, records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--wgs-genotypes",
        type=Path,
        default=Path("data/nyc_coydog_breeds/all_sample_genotypes.json"),
    )
    parser.add_argument(
        "--out-dir", type=Path, default=Path("data/nyc_coydog_breeds")
    )
    parser.add_argument("--min-group", type=int, default=3)
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    samples, records = _load_panel(args.wgs_genotypes)
    keys = sorted(records, key=lambda item: (int(item[0].removeprefix("chr")), item[1]))

    groups: dict[str, list[int]] = defaultdict(list)
    for index, sample in enumerate(samples):
        groups[breed_of(sample)].append(index)
    categories = {group: classify_group(group) for group in groups}

    candidates: dict[str, list[int]] = {
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
    if len(village_indices) >= args.min_group:
        candidates["VillageDog(all regions)"] = village_indices
    for code, region in INFERRED_VILLAGE_CODES.items():
        if len(groups.get(code, [])) >= args.min_group:
            candidates[f"VillageDog({region}, inferred)"] = groups[code]
    pooled_dogs = [
        index
        for group, indices in groups.items()
        if categories[group] != "wild"
        for index in indices
    ]
    panels = {
        group: allele_frequencies(keys, records, indices)
        for group, indices in candidates.items()
    }

    mixed_indices = [
        index for index, sample in enumerate(samples)
        if classify_group(breed_of(sample)) == "mixed"
    ]

    rows: list[dict[str, object]] = []
    manifests: dict[str, object] = {}
    for index in mixed_indices:
        # Leave the held-out mixed dog out of the pooled any-dog baseline.
        any_panel = allele_frequencies(
            keys, records, [i for i in pooled_dogs if i != index]
        )
        result = score_breed_candidates(
            keys, records, index, samples, panels, any_panel, top_k=args.top_k
        )
        rows.append({
            "sample_id": result["sample_id"],
            "best_breed": result["best_breed"],
            "single_breed_gap": result["single_breed_gap"],
            "single_breed_supported": result["single_breed_supported"],
            "loci_used": result["loci_used"],
        })
        manifests[str(result["sample_id"])] = result

    csv_path = args.out_dir / "mixed_breed_sensitivity.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "method": {
            "model": "candidate-breed log-likelihood vs pooled any-dog baseline",
            "single_breed_gap_threshold": SINGLE_BREED_GAP,
            "n_candidate_groups": len(candidates),
            "loci": len(keys),
            "note": (
                "a true mixed dog should show single_breed_supported=false; a "
                "confident call on a sample labeled mixed/unknown suggests a label issue"
            ),
        },
        "results": manifests,
    }
    (args.out_dir / "mixed_breed_sensitivity.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    for row in rows:
        print(
            f"{row['sample_id']}: best={row['best_breed']!r} "
            f"gap={row['single_breed_gap']} supported={row['single_breed_supported']}"
        )
    print(f"Wrote {csv_path}")


if __name__ == "__main__":
    main()
