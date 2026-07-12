"""Smoke tests for the local browser interface configuration layer."""

from __future__ import annotations

from pathlib import Path

import pytest

from canidae.ui import UiRequestError, build_ui_config, estimate_laptop_run


def test_public_ui_config_starts_with_integrated_reduced_panel(tmp_path: Path) -> None:
    cfg = build_ui_config({
        "dataset_mode": "redwolf_public",
        "panel_preset": "10k",
        "sample_sheet": "samples.csv",
        "selected_samples": "Wolf25, Wolf26",
        "analyses": ["tree", "introgression"],
    }, tmp_path)

    assert cfg.pipeline[0] == "reduced_panel"
    assert {"qc", "load_genotypes", "distance", "nj_tree", "dstats", "report"} <= set(cfg.pipeline)
    assert cfg.stage_config("reduced_panel")["selected_samples"] == ["Wolf25", "Wolf26"]


def test_local_ui_config_requires_input_paths(tmp_path: Path) -> None:
    with pytest.raises(UiRequestError, match="VCF/BCF"):
        build_ui_config({"dataset_mode": "local"}, tmp_path)


def test_laptop_estimate_mentions_exact_preflight() -> None:
    estimate = estimate_laptop_run("25k")
    assert estimate["download"] == "about 2-8 GB"
    assert "exact indexed-range estimate" in estimate["note"]
