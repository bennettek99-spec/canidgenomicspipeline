"""Unit tests for the hybrid-diagnostic figure functions.

These plot functions run headless (Agg) so they must produce a valid PNG from
the CSV shapes the hybrid stages actually emit, including the missing-estimate
cases (``NY04`` has no dog fraction; breed gaps are signed and can be negative).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from canidae.stages.reporting import figures


def _png(path: Path) -> bytes:
    assert path.exists(), f"figure not written: {path}"
    data = path.read_bytes()
    assert data.startswith(b"\x89PNG"), "output is not a PNG"
    return data


def test_dog_fraction_barplot_writes_a_png(tmp_path: Path) -> None:
    csv = tmp_path / "reference_mixture.csv"
    pd.DataFrame(
        [
            {"sample_id": "T211", "dog_fraction": 0.235, "diagnostic_sites_called": 21},
            {"sample_id": "NY04", "dog_fraction": None, "diagnostic_sites_called": 7},
            {"sample_id": "NY01", "dog_fraction": 0.402, "diagnostic_sites_called": 32},
        ]
    ).to_csv(csv, index=False)
    out = figures.dog_fraction_barplot(csv, tmp_path / "mixture.png")
    assert _png(out)


def test_dog_fraction_barplot_handles_all_missing(tmp_path: Path) -> None:
    csv = tmp_path / "reference_mixture.csv"
    pd.DataFrame(
        [
            {"sample_id": "NY04", "dog_fraction": None, "diagnostic_sites_called": 7},
        ]
    ).to_csv(csv, index=False)
    out = figures.dog_fraction_barplot(csv, tmp_path / "mixture.png")
    assert _png(out)


def test_breed_gap_barplot_writes_a_png(tmp_path: Path) -> None:
    csv = tmp_path / "breed_assign.csv"
    pd.DataFrame(
        [
            {"sample_id": "NY01", "single_breed_gap": -0.53, "single_breed_supported": False},
            {"sample_id": "T211", "single_breed_gap": 0.33, "single_breed_supported": False},
        ]
    ).to_csv(csv, index=False)
    out = figures.breed_gap_barplot(csv, tmp_path / "breed_gap.png", threshold=5.0)
    assert _png(out)


def test_breed_ranking_barplot_writes_a_png(tmp_path: Path) -> None:
    scores = tmp_path / "breed_scores.csv"
    pd.DataFrame(
        [
            {
                "sample_id": "NYFIXTURE_F1",
                "rank": 1,
                "breed": "Beagle",
                "log_likelihood": -10.0,
                "gap_to_next": 12.0,
                "loci_used": 285,
            },
            {
                "sample_id": "NYFIXTURE_F1",
                "rank": 2,
                "breed": "LabradorRetriever",
                "log_likelihood": -22.0,
                "gap_to_next": None,
                "loci_used": 285,
            },
            {
                "sample_id": "NYFIXTURE_OFF1",
                "rank": 1,
                "breed": "Beagle",
                "log_likelihood": -30.0,
                "gap_to_next": 0.5,
                "loci_used": 285,
            },
        ]
    ).to_csv(scores, index=False)
    primary = tmp_path / "breed_assign.csv"
    pd.DataFrame(
        [
            {"sample_id": "NYFIXTURE_F1", "any_dog_log_likelihood": -22.5},
            {"sample_id": "NYFIXTURE_OFF1", "any_dog_log_likelihood": -30.0},
        ]
    ).to_csv(primary, index=False)
    out = figures.breed_ranking_barplot(scores, primary, tmp_path / "breed_ranking.png")
    assert _png(out)


def test_breed_ranking_barplot_handles_empty_scores(tmp_path: Path) -> None:
    scores = tmp_path / "breed_scores.csv"
    pd.DataFrame(columns=["sample_id", "rank", "breed", "log_likelihood"]).to_csv(
        scores, index=False
    )
    primary = tmp_path / "breed_assign.csv"
    pd.DataFrame(columns=["sample_id", "any_dog_log_likelihood"]).to_csv(primary, index=False)
    out = figures.breed_ranking_barplot(scores, primary, tmp_path / "breed_ranking.png")
    assert _png(out)


def test_admixture_scatter_writes_a_png(tmp_path: Path) -> None:
    csv = tmp_path / "multiway_admixture.csv"
    pd.DataFrame(
        [
            {
                "sample_id": "Coyote001",
                "region": "western",
                "f_wolf": 0.0,
                "f_dog": 0.0,
                "f_wolf_ci_lo": 0.0,
                "f_wolf_ci_hi": 0.0,
                "f_dog_ci_lo": 0.0,
                "f_dog_ci_hi": 0.0,
            },
            {
                "sample_id": "Coyote168",
                "region": "eastern",
                "f_wolf": 0.11,
                "f_dog": 0.0,
                "f_wolf_ci_lo": 0.0,
                "f_wolf_ci_hi": 0.2,
                "f_dog_ci_lo": 0.0,
                "f_dog_ci_hi": 0.05,
            },
        ]
    ).to_csv(csv, index=False)
    out = figures.admixture_scatter(csv, tmp_path / "multiway_scatter.png")
    assert _png(out)
