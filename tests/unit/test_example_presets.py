"""Every shipped example preset and metadata sheet must stay loadable.

A public preset is a promise: the YAML parses, every stage it names is
registered, every stage config block validates, and any sample sheet it ships
alongside satisfies the sheet schema. These are cheap checks, but they are the
ones that break silently when a stage config gains a field or a stage is
renamed.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from canidae.core.config import GlobalConfig
from canidae.pipeline import instantiate_stages
from canidae.stages.acquisition.sample_sheet import read_sample_sheet

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "configs" / "examples"
PRESETS = sorted(EXAMPLES.glob("*.yaml"))
SHEETS = sorted(EXAMPLES.glob("*.csv"))


def test_the_repository_ships_example_presets() -> None:
    assert PRESETS, "expected recipes under configs/examples/"
    assert SHEETS, "expected sample sheets under configs/examples/"


@pytest.mark.parametrize("preset", PRESETS, ids=lambda p: p.stem)
def test_preset_instantiates_every_stage_it_names(preset: Path, tmp_path: Path) -> None:
    cfg = GlobalConfig.load(preset, overrides={"paths.root": str(tmp_path)})
    assert cfg.pipeline, f"{preset.name} declares an empty pipeline"
    stages = instantiate_stages(cfg)
    assert [stage.name for stage in stages] == list(cfg.pipeline)


@pytest.mark.parametrize("preset", PRESETS, ids=lambda p: p.stem)
def test_preset_configures_only_stages_it_runs(preset: Path) -> None:
    """A stages: block for a stage absent from pipeline: is silently dead config."""
    raw = yaml.safe_load(preset.read_text(encoding="utf-8")) or {}
    configured = set(raw.get("stages") or {})
    declared = set(raw.get("pipeline") or [])
    assert configured <= declared, (
        f"{preset.name} configures stages it never runs: {sorted(configured - declared)}"
    )


@pytest.mark.parametrize("sheet", SHEETS, ids=lambda p: p.stem)
def test_sample_sheet_validates(sheet: Path) -> None:
    df = read_sample_sheet(sheet)
    assert not df.empty
    assert df["sample_id"].is_unique
    assert (df["taxon"] != "").all()
    assert (df["population"] != "").all()


def test_nyc_sheet_documents_exactly_the_preset_query_samples() -> None:
    """The NYC metadata sheet carries the accessions behind the preset's queries."""
    sheet = read_sample_sheet(EXAMPLES / "nyc_coydog_samples.csv")
    raw = yaml.safe_load((EXAMPLES / "nyc_coydog_validation.yaml").read_text(encoding="utf-8"))
    queries = set(raw["stages"]["reference_mixture"]["query_samples"])
    assert set(sheet["sample_id"]) == queries
    assert set(raw["stages"]["breed_assign"]["query_samples"]) == queries
    # Each query is traceable back to the archive run it came from. These must
    # be the schema's own column names, or read_sample_sheet drops the values.
    assert (sheet["accession"] != "").all()
    assert (sheet["source_study"] != "").all()
