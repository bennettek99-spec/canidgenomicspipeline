from __future__ import annotations

from pathlib import Path

import pytest

from canidae.core.errors import StageInputError
from canidae.core.model import CanidTaxon
from canidae.stages.acquisition.sample_sheet import read_sample_sheet, to_individuals

_GOOD = (
    "sample_id,taxon,population,region,latitude,longitude\n"
    "W1,gray_wolf,eurasia,siberia,61.0,105.0\n"
    "C1,coyote,north_america,,39.0,-98.0\n"
)


def _write(tmp_path: Path, text: str) -> Path:
    p = tmp_path / "samples.csv"
    p.write_text(text, encoding="utf-8")
    return p


def test_valid_sheet_parses(tmp_path: Path) -> None:
    df = read_sample_sheet(_write(tmp_path, _GOOD))
    assert list(df["sample_id"]) == ["W1", "C1"]
    inds = to_individuals(df)
    assert inds[0].taxon.label is CanidTaxon.GRAY_WOLF
    assert inds[0].locality.latitude == 61.0
    assert inds[1].locality.longitude == -98.0


def test_missing_required_column(tmp_path: Path) -> None:
    with pytest.raises(StageInputError, match="required column"):
        read_sample_sheet(_write(tmp_path, "sample_id,taxon\nW1,gray_wolf\n"))


def test_duplicate_sample_id(tmp_path: Path) -> None:
    text = "sample_id,taxon,population\nW1,gray_wolf,a\nW1,coyote,b\n"
    with pytest.raises(StageInputError, match="duplicate"):
        read_sample_sheet(_write(tmp_path, text))


def test_unknown_taxon_rejected(tmp_path: Path) -> None:
    text = "sample_id,taxon,population\nX1,unicorn,a\n"
    with pytest.raises(StageInputError, match="taxon"):
        read_sample_sheet(_write(tmp_path, text))


def test_out_of_range_coordinate(tmp_path: Path) -> None:
    text = "sample_id,taxon,population,latitude\nW1,gray_wolf,a,999\n"
    with pytest.raises(StageInputError, match="latitude"):
        read_sample_sheet(_write(tmp_path, text))


def test_missing_file(tmp_path: Path) -> None:
    with pytest.raises(StageInputError, match="not found"):
        read_sample_sheet(tmp_path / "nope.csv")


def test_sample_sheet_stage_emits_artifact(tmp_context, tmp_path: Path) -> None:
    from canidae.core.model import ArtifactKind
    from canidae.stages.acquisition.ingest import SampleSheetConfig, SampleSheetStage

    csv = _write(tmp_path, _GOOD)
    stage = SampleSheetStage(SampleSheetConfig(path=csv))
    result = stage.run(tmp_context)
    for art in result.artifacts:
        tmp_context.datastore.register(art)
    assert tmp_context.datastore.has(ArtifactKind.SAMPLE_SHEET, "sample_sheet")
    assert result.metrics["n_samples"] == 2
