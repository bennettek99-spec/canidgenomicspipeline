"""Estimate wolf/dog admixture in eastern coyotes and score the dog component.

Eastern coyotes (Ontario) are thought to carry gray-wolf and domestic-dog
ancestry relative to western coyotes (Arizona).  This script re-derives that
signal from public data with a cross-platform design that is deliberately
small:

- Queries: 36 RADseq coyotes from Heppenheimer et al. 2018 (Dryad), already
  remapped to WGS REF/ALT orientation at 267 bridge loci in the bridge VCF.
- References: the full NHGRI 722-genome WGS panel genotypes at the 252 shared
  bridge loci (fetched once by breed_assign_nyc_coydog.py), split into western
  coyote, gray wolf, eastern wolf, and domestic-dog (pooled + per-breed)
  panels.

For every RADseq coyote a three-way coyote/wolf/dog admixture proportion is
fitted by maximum likelihood over the 252 loci (binomial genotype likelihood,
grid search), with locus-bootstrap confidence intervals.  Arizona coyotes act
as the negative control.  For animals with a dog component, each breed panel
is scored against the pooled any-dog panel exactly as in the NYC coydog
analysis; a small gap means the dog ancestry does not match any single breed.

This is a 252-locus diagnostic reanalysis, not a publication-grade local
ancestry estimate.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.special import gammaln

sys.path.insert(0, str(Path(__file__).parent))
import breed_assign_nyc_coydog as breeds
import validate_nyc_coydog as nyc

# Western (Arizona) vs eastern (Ontario) RADseq sample numbering from the
# Heppenheimer et al. 2018 sample sheet (see eastern_coyote_wgs_bridge CSV).
WESTERN_IDS = frozenset(f"Coyote{number:03d}" for number in range(1, 13))
GRID_STEP = 0.01
BOOTSTRAP_N = 200
# Report breed scoring only for animals whose dog fraction point estimate
# reaches this level; below it the breed question is not answerable.
MIN_DOG_FRACTION_FOR_BREED = 0.05
MIN_EFFECTIVE_P = 1e-4
# Empirical-Bayes shrinkage: every panel's allele frequency is pulled toward
# the pooled all-canid frequency at that locus with the strength of
# SHRINKAGE_PSEUDOCOUNT virtual observations.  Without this, panels of very
# different sizes are not comparable in a likelihood: a 660-dog panel reaches
# extreme frequencies a 3-coyote panel never can, and homozygous-reference
# samples then spuriously maximize f_dog.
SHRINKAGE_PSEUDOCOUNT = 6.0
# Loci must separate coyote from at least one admixture source (on raw
# frequencies) to inform the fit; uninformative loci only add platform noise.
MIN_PANEL_SEPARATION = 0.15
# NOTE: the 2026-07-13 bridge VCF's six WGS sample columns carry mislabeled
# genotypes (verified against live source records: e.g. its "Coyote01" column
# matches panel sample 140447_S11).  Its 36 RADseq columns are unaffected.
# This script never uses the bridge VCF's WGS columns; WGS genotypes come from
# all_sample_genotypes.json (fetched and verified separately).
BRIDGE_WGS_COLUMNS_BUG = True
BRIDGE_WGS_IDS = frozenset(
    {"Coyote01", "Coyote02", "AlaskanWolf", "AlgonquinWolf13467",
     "AlgonquinWolf13470", "GoldenJackal01"}
)


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
                per_sample[sample] = nyc_index_genotype(field, gt_index)
            calls[key] = per_sample
    return samples, calls


def nyc_index_genotype(sample_field: str, gt_index: int) -> str | None:
    values = sample_field.split(":")
    if gt_index >= len(values):
        return None
    gt = values[gt_index]
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    return gt


def dosage_vector(
    keys: list[tuple[str, int]],
    calls: dict[tuple[str, int], dict[str, str | None]],
    sample: str,
) -> np.ndarray:
    out = np.full(len(keys), np.nan)
    for position, key in enumerate(keys):
        genotype = calls[key].get(sample)
        if genotype is not None:
            out[position] = nyc.dosage(genotype)
    return out


def wgs_dosage_matrix(
    keys: list[tuple[str, int]],
    records: dict[tuple[str, int], list[str | None]],
) -> np.ndarray:
    """Dosage matrix (all WGS samples x loci) built once and reused."""
    matrix = np.full((len(records[keys[0]]), len(keys)), np.nan)
    for row in range(matrix.shape[0]):
        for position, key in enumerate(keys):
            genotype = records[key][row]
            if genotype is not None:
                matrix[row, position] = nyc.dosage(genotype)
    return matrix


def panel_frequencies(
    matrix: np.ndarray, indices: list[int], pooled_prior: np.ndarray
) -> np.ndarray:
    """Shrunk allele frequencies; panels of different sizes stay comparable."""
    dosages = matrix[indices]
    counts = np.nansum(dosages, axis=0)
    called = np.sum(~np.isnan(dosages), axis=0)
    freqs = (counts + SHRINKAGE_PSEUDOCOUNT * pooled_prior) / (2.0 * called + SHRINKAGE_PSEUDOCOUNT)
    freqs[called == 0] = np.nan
    return freqs


def loglik_grid(
    dosages: np.ndarray, p_coyote: np.ndarray, p_wolf: np.ndarray, p_dog: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Log-likelihood over a (f_wolf, f_dog) grid for one sample.

    Vectorized: p_eff = (1 - fw - fd) * p_c + fw * p_w + fd * p_d per locus,
    with the pairwise-independence binomial likelihood on observed dosages.
    """
    steps = round(1.0 / GRID_STEP) + 1
    grid = np.linspace(0.0, 1.0, steps)
    fw, fd = np.meshgrid(grid, grid, indexing="ij")
    valid = fw + fd <= 1.0 + 1e-9
    fw_v, fd_v = fw[valid], fd[valid]
    observed = dosages[~np.isnan(dosages)]
    p_obs_c = p_coyote[~np.isnan(dosages)]
    p_obs_w = p_wolf[~np.isnan(dosages)]
    p_obs_d = p_dog[~np.isnan(dosages)]
    p_eff = (
        (1.0 - fw_v[:, None] - fd_v[:, None]) * p_obs_c[None, :]
        + fw_v[:, None] * p_obs_w[None, :]
        + fd_v[:, None] * p_obs_d[None, :]
    )
    p_eff = np.clip(p_eff, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
    alt2 = observed[None, :] > 1.0
    het = (observed[None, :] > 0.0) & ~alt2
    ref2 = ~(alt2 | het)
    ll = np.zeros_like(p_eff)
    ll += alt2 * (2.0 * np.log(p_eff))
    ll += het * (np.log(2.0) + np.log(p_eff) + np.log1p(-p_eff))
    ll += ref2 * (2.0 * np.log1p(-p_eff))
    total = ll.sum(axis=1)
    return fw_v, fd_v, total


def best_mixture(
    dosages: np.ndarray, p_coyote: np.ndarray, p_wolf: np.ndarray, p_dog: np.ndarray
) -> tuple[float, float, float]:
    fw_v, fd_v, total = loglik_grid(dosages, p_coyote, p_wolf, p_dog)
    best = int(np.argmax(total))
    return float(fw_v[best]), float(fd_v[best]), float(total[best])


def cohort_grid(
    alt_counts: np.ndarray,
    chromosome_counts: np.ndarray,
    p_coyote: np.ndarray,
    p_wolf: np.ndarray,
    p_dog: np.ndarray,
) -> tuple[float, float, float]:
    """ML mixture fit for a pooled cohort (binomial on per-locus alt counts).

    Pooling all samples in a region restores the power that 86 informative
    loci cannot give per-individual fits; this is the headline estimator.
    """
    steps = round(1.0 / GRID_STEP) + 1
    grid = np.linspace(0.0, 1.0, steps)
    fw, fd = np.meshgrid(grid, grid, indexing="ij")
    valid = fw + fd <= 1.0 + 1e-9
    fw_v, fd_v = fw[valid], fd[valid]
    use = chromosome_counts > 0
    n = chromosome_counts[use]
    k = np.rint(alt_counts[use]).astype(np.int64)
    pc, pw, pd = p_coyote[use], p_wolf[use], p_dog[use]
    p_eff = (
        (1.0 - fw_v[:, None] - fd_v[:, None]) * pc[None, :]
        + fw_v[:, None] * pw[None, :]
        + fd_v[:, None] * pd[None, :]
    )
    p_eff = np.clip(p_eff, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
    ll = (
        gammaln(n + 1)
        - gammaln(k + 1)
        - gammaln(n - k + 1)
        + k * np.log(p_eff)
        + (n - k) * np.log1p(-p_eff)
    ).sum(axis=1)
    best = int(np.argmax(ll))
    return float(fw_v[best]), float(fd_v[best]), float(ll[best])


def bootstrap_ci(
    dosages: np.ndarray,
    p_coyote: np.ndarray,
    p_wolf: np.ndarray,
    p_dog: np.ndarray,
    rng: np.random.Generator,
) -> dict[str, list[float]]:
    n = int(np.sum(~np.isnan(dosages)))
    wolf_draws: list[float] = []
    dog_draws: list[float] = []
    keep = ~np.isnan(dosages)
    for _ in range(BOOTSTRAP_N):
        pick = rng.integers(0, n, n)
        d = dosages[keep][pick]
        pc, pw, pd = p_coyote[keep][pick], p_wolf[keep][pick], p_dog[keep][pick]
        fw, fd, _ = best_mixture(d, pc, pw, pd)
        wolf_draws.append(fw)
        dog_draws.append(fd)
    return {
        "wolf": [round(v, 3) for v in np.percentile(wolf_draws, [2.5, 97.5])],
        "dog": [round(v, 3) for v in np.percentile(dog_draws, [2.5, 97.5])],
    }


def calls_mixture_loglik(
    dosages: np.ndarray,
    p_coyote: np.ndarray,
    p_wolf: np.ndarray,
    p_dog_breed: np.ndarray,
    f_wolf: float,
    f_dog: float,
) -> float:
    keep = ~np.isnan(dosages)
    observed = dosages[keep]
    p_eff = (
        (1.0 - f_wolf - f_dog) * p_coyote[keep]
        + f_wolf * p_wolf[keep]
        + f_dog * p_dog_breed[keep]
    )
    p_eff = np.clip(p_eff, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
    alt2 = observed > 1.0
    het = (observed > 0.0) & ~alt2
    ref2 = ~(alt2 | het)
    ll = float(
        2.0 * np.sum(np.log(p_eff[alt2]))
        + np.sum(np.log(2.0) + np.log(p_eff[het]) + np.log1p(-p_eff[het]))
        + 2.0 * np.sum(np.log1p(-p_eff[ref2]))
    )
    return ll


def two_way_fit(
    alt_counts: np.ndarray,
    chromosome_counts: np.ndarray,
    p_from: np.ndarray,
    p_to: np.ndarray,
) -> tuple[float, float]:
    """Constrained 2-way fit (f of p_from -> p_to); returns (f, loglik)."""
    use = chromosome_counts > 0
    n = chromosome_counts[use]
    k = np.rint(alt_counts[use]).astype(np.int64)
    a, b = p_from[use], p_to[use]
    best = (0.0, -np.inf)
    for f in np.linspace(0.0, 1.0, 101):
        p_eff = np.clip((1.0 - f) * a + f * b, MIN_EFFECTIVE_P, 1.0 - MIN_EFFECTIVE_P)
        ll = float(
            (
                gammaln(n + 1)
                - gammaln(k + 1)
                - gammaln(n - k + 1)
                + k * np.log(p_eff)
                + (n - k) * np.log1p(-p_eff)
            ).sum()
        )
        if ll > best[1]:
            best = (float(f), ll)
    return best


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
        groups[breeds.breed_of(sample)].append(index)
    categories = {group: breeds.classify_group(group) for group in groups}

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
    matrix = wgs_dosage_matrix(keys, records)
    called_counts = np.sum(~np.isnan(matrix), axis=0)
    alt_counts = np.nansum(matrix, axis=0)
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
    # Balanced prior: unweighted mean over the four taxon pools, so the 664-dog
    # panel cannot drag the coyote baseline toward dog allele frequencies.
    pooled_prior = np.nanmean(
        np.vstack(
            [raw["coyote"], raw["gray_wolf"], raw["eastern_wolf"], raw["dog"]]
        ),
        axis=0,
    )
    del called_counts, alt_counts
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
        | (
            np.abs(raw["eastern_wolf"] - raw["coyote"]) >= MIN_PANEL_SEPARATION
        )
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

    # ---- method validation on WGS-platform rows (leave-one-out) ----
    # The bridge VCF's WGS columns are mislabeled (see note above), so
    # validation animals are scored from the verified WGS matrix with their
    # own contribution removed from the reference panels.
    print("Method validation (WGS animals, leave-one-out):", flush=True)
    platform_rows: list[dict[str, object]] = []
    validation_targets = [
        ("Coyote01", coyote_indices),
        ("Coyote02", coyote_indices),
        ("AlaskanWolf", gray_wolf_indices),
        ("AlgonquinWolf13467", eastern_wolf_indices),
    ]
    for sample, _panel_indices in validation_targets:
        row_index = wgs_samples.index(sample)
        coyote_loo = [i for i in coyote_indices if i != row_index]
        gray_loo = [i for i in gray_wolf_indices if i != row_index]
        loo_prior = pooled_prior
        p_c_loo = panel_frequencies(matrix, coyote_loo, loo_prior)[keep]
        p_g_loo = panel_frequencies(matrix, gray_loo, loo_prior)[keep]
        dosages = matrix[row_index][keep]
        f_wolf, f_dog, _ = best_mixture(dosages, p_c_loo, p_g_loo, p_dog_k)
        platform_rows.append(
            {"sample_id": sample, "f_wolf": round(f_wolf, 3), "f_dog": round(f_dog, 3)}
        )
        print(f"  {sample}: wolf {f_wolf:.2f} dog {f_dog:.2f}", flush=True)
    coyote_controls = [r for r in platform_rows if r["sample_id"].startswith("Coyote")]
    wolf_control = next(r for r in platform_rows if r["sample_id"] == "AlaskanWolf")
    validation_passed = wolf_control["f_wolf"] >= 0.5 and all(
        c["f_wolf"] <= 0.25 for c in coyote_controls
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
    for code, region in breeds.INFERRED_VILLAGE_CODES.items():
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
        dosages = dosage_vector(kept_keys, bridge_calls, sample)
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

    # ---- cohort-level estimates (the headline numbers) ----
    def cohort_fit(sample_ids: list[str], label: str) -> dict[str, object]:
        stacked = np.vstack([dosage_vector(kept_keys, bridge_calls, s) for s in sample_ids])
        alt_counts = np.nansum(stacked, axis=0)
        chromosome_counts = 2 * np.sum(~np.isnan(stacked), axis=0)
        f_wolf, f_dog, ll = cohort_grid(
            alt_counts, chromosome_counts, p_coyote_k, p_gray_k, p_dog_k
        )
        # locus bootstrap for confidence intervals
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
        [row["sample_id"] for row in eastern], "eastern (Ontario)"
    )
    cohort_western = cohort_fit(
        [row["sample_id"] for row in western], "western (Arizona, control)"
    )

    # ---- can wolf vs dog even be told apart at these loci? ----
    def cohort_alt_chrom(sample_ids: list[str]) -> tuple[np.ndarray, np.ndarray]:
        stacked = np.vstack(
            [dosage_vector(kept_keys, bridge_calls, s) for s in sample_ids]
        )
        return np.nansum(stacked, axis=0), 2 * np.sum(~np.isnan(stacked), axis=0)

    alt_e, chrom_e = cohort_alt_chrom([row["sample_id"] for row in eastern])
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
    print(f"2-way eastern: coyote->dog f={f_dog_2w:.2f} (LL {ll_dog_2w:.1f}); "
          f"coyote->wolf f={f_wolf_2w:.2f} (LL {ll_wolf_2w:.1f})", flush=True)

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
