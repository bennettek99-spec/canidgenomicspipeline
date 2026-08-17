"""Pure analysis helpers for hybrid-canid diagnostics.

Library code extracted from the NYC coydog / eastern coyote / breed-assignment
scripts so stages can reuse it without depending on ``scripts/``.
"""

from __future__ import annotations

from canidae.analysis.breed_panel import (
    DOMESTIC_EXCEPTIONS,
    INFERRED_VILLAGE_CODES,
    UNRESOLVED_CODES,
    VILLAGE_KEYWORDS,
    VILLAGE_PREFIXES,
    WILD_KEYWORDS,
    breed_of,
    classify_group,
)
from canidae.analysis.genotypes import allele_count, dosage, parse_gt_field
from canidae.analysis.reference_mixture import (
    MIN_EFFECTIVE_P,
    SINGLE_BREED_GAP,
    allele_frequencies,
    binom2_logpmf,
    calls_log_likelihood,
    infer_dog_fraction,
    leave_one_out,
    log_likelihood,
)

__all__ = [
    "DOMESTIC_EXCEPTIONS",
    "INFERRED_VILLAGE_CODES",
    "MIN_EFFECTIVE_P",
    "SINGLE_BREED_GAP",
    "UNRESOLVED_CODES",
    "VILLAGE_KEYWORDS",
    "VILLAGE_PREFIXES",
    "WILD_KEYWORDS",
    "allele_count",
    "allele_frequencies",
    "binom2_logpmf",
    "breed_of",
    "calls_log_likelihood",
    "classify_group",
    "dosage",
    "infer_dog_fraction",
    "leave_one_out",
    "log_likelihood",
    "parse_gt_field",
]
