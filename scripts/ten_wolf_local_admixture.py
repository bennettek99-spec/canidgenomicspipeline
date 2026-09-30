"""Recompute shared-panel canid contrasts for requested North American wolf groups.

This script only scores groups with identified genomes in data/wolf_cline/panel.npz.
It never substitutes an unlabeled modern wolf for historical C. l. nubilus.
Run: python3 scripts/ten_wolf_local_admixture.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from canidae.analysis.f_statistics import JackknifeEstimate, d_statistic, f4_ratio, gl_population_frequencies

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "data" / "wolf_cline" / "panel.npz"
MANIFEST = ROOT / "configs" / "examples" / "wolf_cline_samples.csv"
OUT = ROOT / "results" / "north_american_wolf_admixture"

TARGETS = {
    "Great Lakes wolf": ["Wolf18", "Wolf40"],
    "Red wolf": ["Wolf25", "Wolf26"],
    "Eastern/Algonquin wolf": ["AlgonquinWolf13467", "AlgonquinWolf13470"],
    "Mexican wolf": ["Wolf22", "Wolf23"],
    "Yellowstone wolf": ["Wolf28", "Wolf29", "Wolf30"],
    "Alaskan wolf": ["AlaskanWolf"],
}
UNSAMPLED = [
    "Historical Plains wolf (C. l. nubilus; C. variabilis Wied-Neuwied, 1841)",
    "Pacific Coast wolf",
    "Atlantic Coast wolf",
    "West Arctic wolf",
]


def fields(estimate: JackknifeEstimate, prefix: str, *, include_z: bool = True) -> dict[str, float | int]:
    value = float(estimate.estimate)
    se = float(estimate.se)
    result: dict[str, float | int] = {
        prefix: round(value, 6),
        f"{prefix}_se": round(se, 6),
        f"{prefix}_sites": int(estimate.n_sites),
    }
    if include_z:
        result[f"{prefix}_z"] = round(value / se, 3) if se > 0 else 0.0
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        manifest = {row["sample_id"]: row for row in csv.DictReader(handle)}
    with np.load(PANEL) as panel:
        samples = [str(sample) for sample in panel["samples"]]
        if set(samples) != set(manifest):
            raise ValueError("Panel samples do not match the sample manifest")
        missing = [sample for members in TARGETS.values() for sample in members if sample not in samples]
        if missing:
            raise ValueError(f"Target genomes missing from panel: {missing}")
        blocks = panel["chrom"].astype(int)
        pl = panel["pl"]
        rows: list[dict[str, object]] = []
        for name, members in TARGETS.items():
            # Leave the target out of the North American gray-wolf comparator.
            western = [s for s in ("Wolf28", "Wolf29", "Wolf30", "AlaskanWolf") if s not in members]
            populations = {
                "target": members,
                "western": western,
                "eurasian": ["Wolf01", "Wolf02", "Wolf03", "Wolf06", "Wolf24", "Wolf27"],
                "coyote": ["Coyote01", "Coyote02", "Coyote03"],
                "jackal": ["GoldenJackal01"],
            }
            freqs = gl_population_frequencies(pl, samples, populations)
            ratio = f4_ratio(freqs, "eurasian", "jackal", "target", "western", "coyote", blocks)
            coyote_d = d_statistic(freqs, "target", "western", "coyote", "jackal", blocks)
            row: dict[str, object] = {
                "group": name,
                "status": "computed",
                "n_genomes": len(members),
                "genomes": ";".join(members),
                "wolf_reference": ";".join(western),
                "note": "Coyote D tests excess allele sharing relative to the listed western wolves. Wolf f4 ratio is model-dependent, especially near 1. Dog ancestry is not inferred by this contrast.",
            }
            row.update(fields(ratio, "wolf_f4_ratio", include_z=False))
            row.update(fields(coyote_d, "coyote_d"))
            rows.append(row)
        for name in UNSAMPLED:
            rows.append({"group": name, "status": "no_population_matched_genotype", "n_genomes": 0, "genomes": "", "wolf_reference": "", "note": "No authenticated target genome in the local panel; no statistic computed."})
    columns = ["group", "status", "n_genomes", "genomes", "wolf_reference", "wolf_f4_ratio", "wolf_f4_ratio_se", "wolf_f4_ratio_sites", "coyote_d", "coyote_d_se", "coyote_d_z", "coyote_d_sites", "note"]
    with (OUT / "local_contrasts.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    (OUT / "local_run.json").write_text(
        json.dumps({"panel": str(PANEL.relative_to(ROOT)), "samples": len(samples), "variants": len(blocks), "jackknife_blocks": len(np.unique(blocks)), "methods": "genotype-likelihood EM population frequencies; f4-ratio and Patterson D; autosome-block jackknife", "targets": {name: members for name, members in TARGETS.items()}, "unsampled": UNSAMPLED, "manifest": str(MANIFEST.relative_to(ROOT)), "interpretation": "A significant positive coyote D supports excess coyote-lineage allele sharing relative to the listed western comparator. The wolf f4 ratio is conditional on the specified reference model; it is not a direct hybridization fraction. This script does not infer dog ancestry."}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {len(rows)} rows to {OUT / 'local_contrasts.csv'}")


if __name__ == "__main__":
    main()
