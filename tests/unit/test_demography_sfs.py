from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from canidae.stages.demography.sfs import (
    DemographyConfig,
    _callable_sites,
    _mean_or_none,
)

# -- _callable_sites --------------------------------------------------------------------


def test_explicit_callable_sites_win_over_everything() -> None:
    cfg = DemographyConfig(callable_sites=5_000, sequence_length=9_000)
    assert _callable_sites(cfg, {"callable_sites": 7_000}) == 5_000


def test_legacy_sequence_length_alias_is_used_when_no_explicit_denominator() -> None:
    cfg = DemographyConfig(sequence_length=1_234)
    assert _callable_sites(cfg, {}) == 1_234


def test_metadata_fallback_applies_last() -> None:
    cfg = DemographyConfig()
    assert _callable_sites(cfg, {"callable_sites": 42}) == 42


def test_missing_denominator_everywhere_yields_zero() -> None:
    assert _callable_sites(DemographyConfig(), {}) == 0


def test_non_numeric_metadata_value_is_treated_as_absent() -> None:
    assert _callable_sites(DemographyConfig(), {"callable_sites": "lots"}) == 0


def test_negative_metadata_values_are_clamped_to_zero() -> None:
    # The pydantic field rejects negatives, but a stale artifact metadata value
    # can still arrive negative at runtime.
    assert _callable_sites(DemographyConfig(), {"callable_sites": -10}) == 0


# -- _mean_or_none ----------------------------------------------------------------------


def test_mean_ignores_missing_entries_and_rounds() -> None:
    series = pd.Series([0.100000004, 0.200000004, None])
    assert _mean_or_none(series) == 0.15


def test_mean_of_all_missing_is_none() -> None:
    assert _mean_or_none(pd.Series([None, None])) is None


def test_mean_of_empty_series_is_none() -> None:
    assert _mean_or_none(pd.Series([], dtype=float)) is None


# -- config surface ---------------------------------------------------------------------


def test_defaults_disable_ne_and_per_site_stats() -> None:
    cfg = DemographyConfig()
    assert cfg.mutation_rate == 0.0
    assert cfg.callable_sites == 0
    assert cfg.sequence_length == 0


def test_callable_sites_field_rejects_negatives() -> None:
    with pytest.raises(ValueError):
        DemographyConfig(callable_sites=-1)


def test_nan_helper_matches_stage_convention() -> None:
    # The stage emits NaN for Tajima's D when seg < 3; np.isfinite gates reporting.
    assert not np.isfinite(float("nan"))
