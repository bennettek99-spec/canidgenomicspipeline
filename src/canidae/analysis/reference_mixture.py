"""Two-source mixture scoring and breed-panel likelihoods on bridge loci."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import numpy as np

from canidae.analysis.genotypes import allele_count, dosage

MIN_EFFECTIVE_P = 1e-4
# Log-likelihood units; below this the best single breed is not considered
# distinguishable from the pooled any-dog panel.
SINGLE_BREED_GAP = 5.0


def allele_frequencies(
    keys: Sequence[tuple[str, int]],
    records: Mapping[tuple[str, int], Sequence[str | None]],
    indices: Sequence[int],
) -> np.ndarray:
    """Laplace-smoothed ALT frequencies per locus for one sample group."""
    freqs = np.full(len(keys), np.nan)
    for position, key in enumerate(keys):
        counts = [
            allele_count(records[key][index])
            for index in indices
            if records[key][index] is not None
        ]
        if counts:
            freqs[position] = (sum(counts) + 1.0) / (2 * len(counts) + 2.0)
    return freqs


def binom2_logpmf(alt_count: float, p: float) -> float:
    """Binomial(2, p) log-pmf of an ALT allele count in {0, 1, 2}."""
    p = min(max(p, MIN_EFFECTIVE_P), 1.0 - MIN_EFFECTIVE_P)
    if alt_count > 1.0:
        return 2.0 * math.log(p)
    if alt_count > 0.0:
        return math.log(2.0) + math.log(p) + math.log1p(-p)
    return 2.0 * math.log1p(-p)


def log_likelihood(
    keys: Sequence[tuple[str, int]],
    records: Mapping[tuple[str, int], Sequence[str | None]],
    sample_index: int,
    panel: np.ndarray,
) -> tuple[float, int]:
    total = 0.0
    used = 0
    for position, key in enumerate(keys):
        p = panel[position]
        if not math.isfinite(p):
            continue
        genotype = records[key][sample_index]
        if genotype is None:
            continue
        total += binom2_logpmf(allele_count(genotype), float(p))
        used += 1
    return total, used


def calls_log_likelihood(
    calls: Mapping[tuple[str, int], tuple[str | None, int, int, int]],
    keys: Sequence[tuple[str, int]],
    panel: np.ndarray,
    coyote_panel: np.ndarray,
    dog_fraction: float,
) -> tuple[float, int]:
    """Mixture likelihood: p_eff = f*p_breed + (1-f)*p_coyote per locus."""
    total = 0.0
    used = 0
    for position, key in enumerate(keys):
        p_breed = panel[position]
        p_coyote = coyote_panel[position]
        if not (math.isfinite(p_breed) and math.isfinite(p_coyote)):
            continue
        genotype, _, _, _ = calls[key]
        if genotype is None:
            continue
        p_eff = dog_fraction * float(p_breed) + (1.0 - dog_fraction) * float(p_coyote)
        total += binom2_logpmf(allele_count(genotype), p_eff)
        used += 1
    return total, used


def leave_one_out(
    keys: Sequence[tuple[str, int]],
    records: Mapping[tuple[str, int], Sequence[str | None]],
    samples: Sequence[str],
    groups: Mapping[str, Sequence[int]],
    min_group: int,
) -> tuple[list[dict[str, object]], dict[str, float | int | None]]:
    eligible = {
        group: list(indices)
        for group, indices in groups.items()
        if len(indices) >= min_group
    }
    panels = {
        group: allele_frequencies(keys, records, indices) for group, indices in eligible.items()
    }
    rows: list[dict[str, object]] = []
    correct = 0
    total = 0
    for group, members in sorted(eligible.items()):
        for sample_index in members:
            own_minus = allele_frequencies(
                keys, records, [i for i in members if i != sample_index]
            )
            scored: list[tuple[float, str]] = []
            for candidate, panel in panels.items():
                effective = own_minus if candidate == group else panel
                ll, used = log_likelihood(keys, records, sample_index, effective)
                if used < 0.9 * len(keys):
                    continue
                scored.append((ll, candidate))
            scored.sort(reverse=True)
            if not scored:
                continue
            sample_id = samples[sample_index]
            prediction = scored[0][1]
            is_correct = prediction == group
            correct += is_correct
            total += 1
            rows.append(
                {
                    "sample_id": sample_id,
                    "true_breed": group,
                    "predicted_breed": prediction,
                    "correct": is_correct,
                    "runner_up": scored[1][1] if len(scored) > 1 else "",
                    "log_likelihood_gap": (
                        round(scored[0][0] - scored[1][0], 2) if len(scored) > 1 else None
                    ),
                }
            )
    summary: dict[str, float | int | None] = {
        "n_tested": total,
        "top1_correct": correct,
        "top1_accuracy": correct / total if total else None,
    }
    return rows, summary


def infer_dog_fraction(
    calls: Mapping[tuple[str, int], tuple[str | None, int, int, int]],
    reference: Mapping[tuple[str, int], Sequence[str]],
    *,
    wgs_samples: Sequence[str],
    coyote_samples: frozenset[str],
    min_difference: float = 0.40,
    min_called: int = 20,
) -> tuple[float | None, int, float | None]:
    """Least-squares dog fraction of a query against coyote vs dog reference means.

    Works on the [0, 1] dosage scale throughout: the query genotype and both
    reference group means share that scale, so the fitted fraction is unaffected
    by it. ``min_difference`` is a dosage-scale separation threshold.
    """
    numerator = denominator = 0.0
    called = 0
    for key, genotypes in reference.items():
        coyote = np.mean(
            [
                dosage(gt)
                for gt, sample in zip(genotypes, wgs_samples, strict=True)
                if sample in coyote_samples
            ]
        )
        dog = np.mean(
            [
                dosage(gt)
                for gt, sample in zip(genotypes, wgs_samples, strict=True)
                if sample not in coyote_samples
            ]
        )
        difference = dog - coyote
        observed = calls[key][0]
        if observed is None or abs(difference) < min_difference:
            continue
        called += 1
        numerator += (dosage(observed) - coyote) * difference
        denominator += difference * difference
    if called < min_called or denominator == 0:
        return None, called, None
    estimate = max(0.0, min(1.0, numerator / denominator))
    # Per-locus residual RMSE is a useful stability diagnostic for this small panel.
    residuals = []
    for key, genotypes in reference.items():
        observed = calls[key][0]
        if observed is None:
            continue
        coyote = np.mean(
            [
                dosage(gt)
                for gt, sample in zip(genotypes, wgs_samples, strict=True)
                if sample in coyote_samples
            ]
        )
        dog = np.mean(
            [
                dosage(gt)
                for gt, sample in zip(genotypes, wgs_samples, strict=True)
                if sample not in coyote_samples
            ]
        )
        if abs(dog - coyote) >= min_difference:
            residuals.append((dosage(observed) - (coyote + estimate * (dog - coyote))) ** 2)
    return estimate, called, math.sqrt(float(np.mean(residuals))) if residuals else None
