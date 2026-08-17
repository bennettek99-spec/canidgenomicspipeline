"""Small genotype parsing helpers shared across hybrid-canid analyses.

Two scales are in play and they are not interchangeable. ``allele_count`` is the
ALT allele count in {0, 1, 2} — the convention the rest of CANIS uses (see
``local_ancestry.hmm`` and ``popgen.plink``) and the one every binomial
likelihood and Laplace-smoothed frequency estimator here expects. ``dosage``
rescales that to [0, 1] and is only for estimators that compare a genotype
against group *mean* dosages on the same scale.
"""

from __future__ import annotations


def allele_count(gt: str) -> int:
    """ALT allele count in {0, 1, 2} for a biallelic diploid genotype string."""
    return sum(int(allele) for allele in gt.replace("|", "/").split("/"))


def dosage(gt: str) -> float:
    """ALT allele dosage in [0, 1] for a biallelic diploid genotype string."""
    return allele_count(gt) / 2.0


def parse_gt_field(sample_field: str, gt_index: int) -> str | None:
    """Extract a biallelic GT from a VCF sample field, or None if missing/invalid."""
    values = sample_field.split(":")
    if gt_index >= len(values):
        return None
    gt = values[gt_index]
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    return gt
