"""Shared test fixtures and lightweight dummy stages.

The dummy stages exercise the Stage/executor contract without needing any external
bioinformatics tool, so the foundation is fully testable in CI on any platform.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.executor import build_context
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.provenance import ProvenanceWriter
from canidae.core.runtime import LocalRunner
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult


def _assert_close(observed: Any, expected: Any, path: str) -> None:
    if isinstance(expected, dict):
        assert isinstance(observed, dict), f"{path}: expected a mapping"
        assert observed.keys() == expected.keys(), f"{path}: keys differ"
        for key, value in expected.items():
            _assert_close(observed[key], value, f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(observed, (list, tuple)), f"{path}: expected a sequence"
        assert len(observed) == len(expected), f"{path}: length differs"
        for index, (got, want) in enumerate(zip(observed, expected, strict=True)):
            _assert_close(got, want, f"{path}[{index}]")
    elif isinstance(expected, float):
        assert observed == pytest.approx(expected, rel=1e-6, abs=2e-6), (
            f"{path}: {observed} != {expected}"
        )
    else:
        assert observed == expected, f"{path}: {observed!r} != {expected!r}"


@pytest.fixture()
def assert_golden_close() -> Callable[[Any, Any], None]:
    """Compare a result to a golden snapshot: exact structure, floats to 1e-6.

    Goldens store floats rounded to six decimals, and the last digit can move between
    BLAS builds or Python versions; a relative 1e-6 tolerance absorbs that and nothing
    more.
    """

    def check(observed: Any, expected: Any) -> None:
        _assert_close(observed, expected, "$")

    return check


@pytest.fixture()
def tmp_context(tmp_path: Path) -> RunContext:
    store = DataStore(tmp_path / "store")
    prov = ProvenanceWriter(tmp_path / "run", config_digest="test", seed=1)
    return build_context(GlobalConfig(), store, LocalRunner(), prov, run_dir=tmp_path / "run")


class ProducerStage(Stage):
    name = "producer"
    config_model = StageConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.QC_TABLE, "metrics")]

    def run(self, ctx: RunContext) -> StageResult:
        out = ctx.datastore.path_for(self.name, "metrics.csv")
        out.write_text("sample,coverage\nA,30\nB,25\n", encoding="utf-8")
        art = Artifact(ArtifactKind.QC_TABLE, "metrics", out, fmt=FileFormat.CSV)
        return StageResult(artifacts=[art], metrics={"rows": 2})


class ConsumerStage(Stage):
    name = "consumer"
    config_model = StageConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.QC_TABLE, "metrics")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "summary")]

    def run(self, ctx: RunContext) -> StageResult:
        inputs = self.gather_inputs(ctx)
        assert inputs and inputs[0].path.exists()
        out = ctx.datastore.path_for(self.name, "summary.txt")
        out.write_text("ok\n", encoding="utf-8")
        art = Artifact(ArtifactKind.ANALYSIS_RESULT, "summary", out, fmt=FileFormat.OTHER)
        return StageResult(artifacts=[art])


@pytest.fixture()
def producer() -> ProducerStage:
    return ProducerStage()


@pytest.fixture()
def consumer() -> ConsumerStage:
    return ConsumerStage()
