"""The paired contrast retains site and chromosome covariance."""

from __future__ import annotations

import numpy as np

from canidae.analysis.f_statistics import paired_f4_ratio_difference


def test_paired_ratio_difference_matches_common_site_ratio() -> None:
    frequencies = {
        "A": np.array([0.8, 0.7, 0.9, 0.8, 0.7, 0.9]),
        "O": np.array([0.1] * 6),
        "B": np.array([0.8] * 6),
        "C": np.array([0.1] * 6),
        "X": np.array([0.7, 0.6, 0.8, 0.7, 0.6, 0.8]),
        "Y": np.array([0.4, 0.3, 0.5, 0.4, np.nan, 0.5]),
    }
    blocks = np.array([1, 1, 2, 2, 3, 3])
    result = paired_f4_ratio_difference(frequencies, "A", "O", "X", "Y", "B", "C", blocks)
    keep = np.isfinite(frequencies["Y"])
    numerator = ((frequencies["A"] - frequencies["O"]) * (frequencies["X"] - frequencies["Y"]))[
        keep
    ]
    denominator = ((frequencies["A"] - frequencies["O"]) * (frequencies["B"] - frequencies["C"]))[
        keep
    ]
    assert np.isclose(result.estimate, numerator.sum() / denominator.sum())
    assert result.n_sites == 5
    assert result.n_blocks == 3
    assert np.isfinite(result.se)
