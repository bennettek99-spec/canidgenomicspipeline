"""Estimate wolf/dog admixture in eastern coyotes and score the dog component.

This is a 252-locus diagnostic reanalysis, not a publication-grade local
ancestry estimate. Pure helpers live in :mod:`canidae.analysis`.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from canidae.analysis.admixture_ml import (
    BOOTSTRAP_N,
    GRID_STEP,
    MIN_DOG_FRACTION_FOR_BREED,
    MIN_PANEL_SEPARATION,
    SHRINKAGE_PSEUDOCOUNT,
    best_mixture,
    bootstrap_ci,
    calls_mixture_loglik,
    cohort_grid,
    panel_frequencies,
    two_way_fit,
    wgs_allele_count_matrix,
)
from canidae.analysis.breed_panel import (
    INFERRED_VILLAGE_CODES,
    breed_of,
    classify_group,
)
from canidae.analysis.bridge_panel import (
    BRIDGE_WGS_COLUMNS_BUG,
    BRIDGE_WGS_IDS,
    allele_count_vector,
    load_bridge_genotypes,
)

# Western (Arizona) vs eastern (Ontario) RADseq sample numbering from the
# Heppenheimer et al. 2018 sample sheet (see eastern_coyote_wgs_bridge CSV).
WESTERN_IDS = frozenset(f"Coyote{number:03d}" for number in range(1, 13))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bridge-vcf",
        type=Path,
        default=Path("data/eastern_coyote_wgs_bridge/eastern_coyote_wgs_bridge_1000.vcf.gz"),
    )
    parser.add_argument(
        "--wgs-genotypes",
        type=Path,
        default=Path("data/nyc_coydog_breeds/all_sample_genotypes.json"),
    )
    parser.add_argument(
        "--wgs-samples-csv",
        type=Path,
        default=Path("data/eastern_coyote_wgs_bridge/eastern_coyote_wgs_bridge_samples.csv"),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("data/eastern_coyote_ancestry"))
    parser.add_argument("--seed", type=int, default=20260816)
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ---- references (WGS platform only) ----
    wgs = json.loads(args.wgs_genotypes.read_text(encoding="utf-8"))
    wgs_samples: list[str] = wgs["samples"]
    records: dict[tuple[str, int], list[str | None]] = {
        (key.rsplit(":", 1)[0], int(key.rsplit(":", 1)[1])): value
        for key, value in wgs["loci"].items()
    }
    keys = sorted(records, key=lambda item: (int(item[0].removeprefix("chr")), item[1]))
    groups: dict[str, list[int]] = defaultdict(list)
    for index, sample in enumerate(wgs_samples):
        groups[breed_of(sample)].append(index)
    categories = {group: classify_group(group) for group in groups}

    coyote_indices = [
        index
        for group, indices in groups.items()
        if categories[group] == "wild" and "coyote" in group.lower()
        for index in indices
    ]
    gray_wolf_indices = [
        index for group in ("Wolf", "AlaskanWolf") for index in groups.get(group, [])
    ]
    eastern_wolf_indices = [
        index
        for group in ("AlgonquinWolf", "QuebecWolf")
        for index in groups.get(group, [])
    ]
    pooled_dogs = [
        index
        for group, indices in groups.items()
        if categories[group] != "wild"
        for index in indices
    ]
    if not (coyote_indices and gray_wolf_indices and pooled_dogs):
        raise RuntimeError(
            f"reference panels incomplete: coyote={len(coyote_indices)}, "
            f"gray_wolf={len(gray_wolf_indices)}, dogs={len(pooled_dogs)}"
        )
    matrix = wgs_allele_count_matrix(keys, records)
    raw = {}
    for name, indices in (
        ("coyote", coyote_indices),
        ("gray_wolf", gray_wolf_indices),
        ("eastern_wolf", eastern_wolf_indices or []),
        ("dog", pooled_dogs),
    ):
        dosages = matrix[indices]
        counts = np.nansum(dosages, axis=0)
        called = np.sum(~np.isnan(dosages), axis=0)
        freqs = (counts + 1.0) / (2.0 * called + 2.0)
        freqs[called == 0] = np.nan
        raw[name] = freqs
    pooled_prior = np.nanmean(
        np.vstack(
            [raw["coyote"], raw["gray_wolf"], raw["eastern_wolf"], raw["dog"]]
        ),
        axis=0,
    )
    p_coyote = panel_frequencies(matrix, coyote_indices, pooled_prior)
    p_gray = panel_frequencies(matrix, gray_wolf_indices, pooled_prior)
    p_eastern = (
        panel_frequencies(matrix, eastern_wolf_indices, pooled_prior)
        if eastern_wolf_indices
        else None
    )
    p_dog = panel_frequencies(matrix, pooled_dogs, pooled_prior)
    informative = (
        (np.abs(raw["gray_wolf"] - raw["coyote"]) >= MIN_PANEL_SEPARATION)
        | (np.abs(raw["dog"] - raw["coyote"]) >= MIN_PANEL_SEPARATION)
        | (np.abs(raw["eastern_wolf"] - raw["coyote"]) >= MIN_PANEL_SEPARATION)
    ) & np.isfinite(p_coyote) & np.isfinite(p_gray) & np.isfinite(p_dog)
    keep = np.where(informative)[0]
    kept_keys = [keys[i] for i in keep]
    print(
        f"References: {len(coyote_indices)} coyotes, {len(gray_wolf_indices)} gray "
        f"wolves, {len(eastern_wolf_indices)} eastern wolves, {len(pooled_dogs)} dogs; "
        f"{len(keep)}/{len(keys)} loci informative (separation >= "
        f"{MIN_PANEL_SEPARATION})",
        flush=True,
    )
    p_coyote_k, p_gray_k, p_dog_k = p_coyote[keep], p_gray[keep], p_dog[keep]
    p_eastern_k = p_eastern[keep] if p_eastern is not None else None

    print("Method validation (WGS animals, leave-one-out):", flush=True)
    platform_rows: list[dict[str, object]] = []
    validation_targets = [
        ("Coyote01", coyote_indices),
        ("Coyote02", coyote_indices),
        ("AlaskanWolf", gray_wolf_indices),
        ("AlgonquinWolf13467", eastern_wolf_indices),
    ]
    for sample, _panel_indices in validation_targets:
        if sample not in wgs_samples:
            continue
        row_index = wgs_samples.index(sample)
        coyote_loo = [i for i in coyote_indices if i != row_index]
        gray_loo = [i for i in gray_wolf_indices if i != row_index]
        p_c_loo = panel_frequencies(matrix, coyote_loo, pooled_prior)[keep]
        p_g_loo = panel_frequencies(matrix, gray_loo, pooled_prior)[keep]
        dosages = matrix[row_index][keep]
        f_wolf, f_dog, _ = best_mixture(dosages, p_c_loo, p_g_loo, p_dog_k)
        platform_rows.append(
            {"sample_id": sample, "f_wolf": round(f_wolf, 3), "f_dog": round(f_dog, 3)}
        )
        print(f"  {sample}: wolf {f_wolf:.2f} dog {f_dog:.2f}", flush=True)
    coyote_controls = [r for r in platform_rows if str(r["sample_id"]).startswith("Coyote")]
    wolf_control = next(
        (r for r in platform_rows if r["sample_id"] == "AlaskanWolf"), None
    )
    validation_passed = bool(
        wolf_control is not None
        and float(wolf_control["f_wolf"]) >= 0.5
        and all(float(c["f_wolf"]) <= 0.25 for c in coyote_controls)
    )
    if not validation_passed:
        print(
            "WARNING: method validation failed; estimates below are biased "
            "and must not be read as ancestry",
            flush=True,
        )

    bridge_samples, bridge_calls = load_bridge_genotypes(args.bridge_vcf)
    radseq_ids = [
        s
        for s in bridge_samples
        if s not in wgs_samples
        and not (BRIDGE_WGS_COLUMNS_BUG and s in BRIDGE_WGS_IDS)
    ]
    print(f"Queries: {len(radseq_ids)} RADseq coyotes at {len(keep)} loci", flush=True)

    candidates: dict[str, list[int]] = {
        group: indices
        for group, indices in groups.items()
        if categories[group] == "breed" and len(indices) >= 3
    }
    village_indices = [
        index
        for group, indices in groups.items()
        if categories[group] == "village"
        for index in indices
    ]
    if len(village_indices) >= 3:
        candidates["VillageDog(all regions)"] = village_indices
    for code, region in INFERRED_VILLAGE_CODES.items():
        if len(groups.get(code, [])) >= 3:
            candidates[f"VillageDog({region}, inferred)"] = groups[code]

    rows: list[dict[str, object]] = []
    breed_rows: list[dict[str, object]] = []
    summary_top: dict[str, object] = {}
    breed_panels = {
        breed: panel_frequencies(matrix, indices, pooled_prior)[keep]
        for breed, indices in candidates.items()
    }
    for sample in radseq_ids:
        region = "western_arizona" if sample in WESTERN_IDS else "eastern_ontario"
        dosages = allele_count_vector(kept_keys, bridge_calls, sample)
        f_wolf, f_dog, ll = best_mixture(dosages, p_coyote_k, p_gray_k, p_dog_k)
        ci = bootstrap_ci(dosages, p_coyote_k, p_gray_k, p_dog_k, rng)
        f_eastern_wolf = None
        if p_eastern_k is not None:
            f_ew, fd_ew, ll_ew = best_mixture(dosages, p_coyote_k, p_eastern_k, p_dog_k)
            f_eastern_wolf = {
                "f_eastern_wolf": f_ew,
                "f_dog": fd_ew,
                "loglik": round(ll_ew, 1),
            }
        rows.append(
            {
                "sample_id": sample,
                "region": region,
                "loci_called": int(np.sum(~np.isnan(dosages))),
                "f_wolf": round(f_wolf, 3),
                "f_dog": round(f_dog, 3),
                "f_wolf_ci": ci["wolf"],
                "f_dog_ci": ci["dog"],
                "loglik": round(ll, 1),
            }
        )
        summary_top[sample] = rows[-1]
        if f_dog >= MIN_DOG_FRACTION_FOR_BREED:
            entries = []
            for breed, panel in breed_panels.items():
                score = calls_mixture_loglik(
                    dosages, p_coyote_k, p_gray_k, panel, f_wolf, f_dog
                )
                entries.append((score, breed))
            any_dog = calls_mixture_loglik(
                dosages, p_coyote_k, p_gray_k, p_dog_k, f_wolf, f_dog
            )
            entries.sort(reverse=True)
            for rank, (score, breed) in enumerate(entries[: args.top_k], start=1):
                breed_rows.append(
                    {
                        "sample_id": sample,
                        "rank": rank,
                        "breed": breed,
                        "log_likelihood": round(score, 2),
                        "gap_to_any_dog": round(score - any_dog, 2),
                    }
                )
        print(
            f"{sample} ({region}): wolf {f_wolf:.2f} {ci['wolf']}, "
            f"dog {f_dog:.2f} {ci['dog']}"
            + (f", eastern-wolf model {f_eastern_wolf}" if f_eastern_wolf else ""),
            flush=True,
        )

    with (args.out_dir / "ancestry_estimates.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    if breed_rows:
        with (args.out_dir / "breed_scores.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(breed_rows[0]))
            writer.writeheader()
            writer.writerows(breed_rows)

    eastern = [row for row in rows if row["region"] == "eastern_ontario"]
    western = [row for row in rows if row["region"] == "western_arizona"]

    def cohort_fit(sample_ids: list[str], label: str) -> dict[str, object]:
        stacked = np.vstack([allele_count_vector(kept_keys, bridge_calls, s) for s in sample_ids])
        alt_counts = np.nansum(stacked, axis=0)
        chromosome_counts = 2 * np.sum(~np.isnan(stacked), axis=0)
        f_wolf, f_dog, ll = cohort_grid(
            alt_counts, chromosome_counts, p_coyote_k, p_gray_k, p_dog_k
        )
        draws = []
        use = chromosome_counts > 0
        for _ in range(BOOTSTRAP_N):
            pick = rng.integers(0, int(use.sum()), int(use.sum()))
            fw_b, fd_b, _ = cohort_grid(
                alt_counts[use][pick],
                chromosome_counts[use][pick],
                p_coyote_k[use][pick],
                p_gray_k[use][pick],
                p_dog_k[use][pick],
            )
            draws.append((fw_b, fd_b))
        draws_arr = np.array(draws)
        result = {
            "n_samples": len(sample_ids),
            "f_wolf": round(f_wolf, 3),
            "f_dog": round(f_dog, 3),
            "f_wolf_ci": [round(v, 3) for v in np.percentile(draws_arr[:, 0], [2.5, 97.5])],
            "f_dog_ci": [round(v, 3) for v in np.percentile(draws_arr[:, 1], [2.5, 97.5])],
            "loglik": round(ll, 1),
        }
        print(f"Cohort {label}: {result}", flush=True)
        return result

    cohort_eastern = cohort_fit(
        [str(row["sample_id"]) for row in eastern], "eastern (Ontario)"
    )
    cohort_western = cohort_fit(
        [str(row["sample_id"]) for row in western], "western (Arizona, control)"
    )

    def cohort_alt_chrom(sample_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
        stacked = np.vstack(
            [allele_count_vector(kept_keys, bridge_calls, s) for s in sample_ids]
        )
        return np.nansum(stacked, axis=0), 2 * np.sum(~np.isnan(stacked), axis=0)

    alt_e, chrom_e = cohort_alt_chrom([str(row["sample_id"]) for row in eastern])
    f_dog_2w, ll_dog_2w = two_way_fit(alt_e, chrom_e, p_coyote_k, p_dog_k)
    f_wolf_2w, ll_wolf_2w = two_way_fit(alt_e, chrom_e, p_coyote_k, p_gray_k)
    wolf_dog_discrimination = {
        "coyote_to_dog": {"f": round(f_dog_2w, 3), "loglik": round(ll_dog_2w, 1)},
        "coyote_to_wolf": {"f": round(f_wolf_2w, 3), "loglik": round(ll_wolf_2w, 1)},
        "loglik_difference": round(ll_dog_2w - ll_wolf_2w, 2),
        "verdict": (
            "wolf and dog panels are interchangeable for explaining the eastern "
            "signal at this marker count; the wolf/dog split is NOT resolvable "
            "here (see full-marker published studies for that question)"
        ),
    }
    print(
        f"2-way eastern: coyote->dog f={f_dog_2w:.2f} (LL {ll_dog_2w:.1f}); "
        f"coyote->wolf f={f_wolf_2w:.2f} (LL {ll_wolf_2w:.1f})",
        flush=True,
    )

    def stats(values: list[float]) -> dict[str, float]:
        return {
            "mean": round(float(np.mean(values)), 3),
            "min": round(float(np.min(values)), 3),
            "max": round(float(np.max(values)), 3),
        }

    def permutation_test(values_e: list[float], values_w: list[float]) -> float:
        observed = float(np.mean(values_e) - np.mean(values_w))
        combined = np.array(values_e + values_w)
        n_e = len(values_e)
        count = 0
        draws = 10_000
        for _ in range(draws):
            rng.shuffle(combined)
            if float(np.mean(combined[:n_e]) - np.mean(combined[n_e:])) >= observed:
                count += 1
        return round((count + 1) / (draws + 1), 4)

    ew_wolf = [float(r["f_wolf"]) for r in eastern]
    ww_wolf = [float(r["f_wolf"]) for r in western]
    ew_dog = [float(r["f_dog"]) for r in eastern]
    ww_dog = [float(r["f_dog"]) for r in western]

    manifest = {
        "created_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "queries": {
            "n_radseq_coyotes": len(radseq_ids),
            "n_eastern_ontario": len(eastern),
            "n_western_arizona": len(western),
        },
        "references": {
            "coyote": [wgs_samples[i] for i in coyote_indices],
            "gray_wolf_n": len(gray_wolf_indices),
            "eastern_wolf_n": len(eastern_wolf_indices),
            "dogs_n": len(pooled_dogs),
            "loci_total": len(keys),
            "loci_informative": len(keep),
        },
        "platform_validation": platform_rows,
        "validation_passed": validation_passed,
        "data_warning": (
            "The 2026-07-13 bridge VCF's six WGS sample columns are mislabeled "
            "(verified against live source records); this analysis uses only its "
            "36 RADseq coyote columns plus independently fetched WGS genotypes."
        ),
        "group_summary": {
            "eastern_wolf_fraction": stats(ew_wolf),
            "western_wolf_fraction": stats(ww_wolf),
            "eastern_dog_fraction": stats(ew_dog),
            "western_dog_fraction": stats(ww_dog),
            "wolf_contrast_p": permutation_test(ew_wolf, ww_wolf),
            "dog_contrast_p": permutation_test(ew_dog, ww_dog),
        },
        "cohort_estimates": {
            "eastern_ontario": cohort_eastern,
            "western_arizona_control": cohort_western,
            "wolf_vs_dog_discrimination": wolf_dog_discrimination,
        },
        "results": summary_top,
        "method": {
            "model": "3-way maximum-likelihood admixture (coyote/gray wolf/dog)",
            "grid_step": GRID_STEP,
            "bootstrap_n": BOOTSTRAP_N,
            "shrinkage_pseudocount": SHRINKAGE_PSEUDOCOUNT,
            "min_panel_separation": MIN_PANEL_SEPARATION,
            "controls": "Arizona (western) coyotes are negative controls; "
            "bridge-VCF dual-platform animals validate the pipeline",
        },
        "limitations": [
            "252 cross-platform bridge loci; intervals are wide by WGS standards.",
            "Wolf and dog panels share derived alleles, so the wolf/dog split at this "
            "marker count is a coarse decomposition, not a definitive partition.",
            "RADseq and WGS genotypes come from different studies and platforms.",
        ],
    }
    (args.out_dir / "eastern_coyote_ancestry.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest["group_summary"], indent=2), flush=True)
    print(f"Wrote {args.out_dir / 'eastern_coyote_ancestry.json'}", flush=True)


if __name__ == "__main__":
    main()
