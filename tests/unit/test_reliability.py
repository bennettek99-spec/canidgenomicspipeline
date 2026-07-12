"""Focused tests for safe resume, atomic stage promotion, and recovery provenance."""

from __future__ import annotations

from pathlib import Path

from pydantic import Field

from canidae.core.executor import NativeExecutor
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.resources import ResourceManager
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.core.stage_cache import StageCache


class _ValueConfig(StageConfig):
    value: int = Field(default=1, ge=0)


class _CountingStage(Stage):
    name = "cache_counter"
    config_model = _ValueConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "cache_counter")]

    def run(self, ctx: RunContext) -> StageResult:
        out = ctx.datastore.path_for(self.name, "counter.txt")
        out.write_text(str(self.config.value), encoding="utf-8")
        return StageResult([
            Artifact(ArtifactKind.ANALYSIS_RESULT, "cache_counter", out, FileFormat.OTHER)
        ])


class _InputStage(Stage):
    name = "input_cache"
    config_model = StageConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.QC_TABLE, "raw_metrics")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "input_summary")]

    def run(self, ctx: RunContext) -> StageResult:
        raw = ctx.datastore.get(ArtifactKind.QC_TABLE, "raw_metrics").path
        out = ctx.datastore.path_for(self.name, "summary.txt")
        out.write_text(raw.read_text(encoding="utf-8"), encoding="utf-8")
        return StageResult([
            Artifact(ArtifactKind.ANALYSIS_RESULT, "input_summary", out, FileFormat.OTHER)
        ])


class _BrokenStage(Stage):
    name = "broken_atomic"
    config_model = StageConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.REPORT, "broken")]

    def run(self, ctx: RunContext) -> StageResult:
        ctx.datastore.path_for(self.name, "partial.html").write_text("partial", encoding="utf-8")
        raise RuntimeError("intentional failure")


def test_safe_resume_skips_only_unchanged_stage(tmp_context: RunContext) -> None:
    executor = NativeExecutor(resume=True)
    first = executor.run([_CountingStage(_ValueConfig(value=1))], tmp_context)
    second = executor.run([_CountingStage(_ValueConfig(value=1))], tmp_context)
    third = executor.run([_CountingStage(_ValueConfig(value=2))], tmp_context)

    assert first.executed == ["cache_counter"]
    assert second.skipped == ["cache_counter"]
    assert third.executed == ["cache_counter"]  # stage configuration changed
    out = tmp_context.datastore.get(ArtifactKind.ANALYSIS_RESULT, "cache_counter").path
    assert out.read_text(encoding="utf-8") == "2"


def test_safe_resume_invalidates_when_an_input_changes(tmp_context: RunContext) -> None:
    raw = tmp_context.datastore.path_for("seed", "raw.csv")
    raw.write_text("before\n", encoding="utf-8")
    tmp_context.datastore.register(
        Artifact(ArtifactKind.QC_TABLE, "raw_metrics", raw, FileFormat.CSV)
    )
    executor = NativeExecutor(resume=True)
    stage = _InputStage()
    assert executor.run([stage], tmp_context).executed == ["input_cache"]
    assert executor.run([stage], tmp_context).skipped == ["input_cache"]

    raw.write_text("after\n", encoding="utf-8")
    assert executor.run([stage], tmp_context).executed == ["input_cache"]
    summary = tmp_context.datastore.get(ArtifactKind.ANALYSIS_RESULT, "input_summary").path
    assert summary.read_text(encoding="utf-8") == "after\n"


def test_safe_resume_fingerprint_includes_stage_code(tmp_context: RunContext, monkeypatch) -> None:
    stage = _CountingStage(_ValueConfig(value=3))
    NativeExecutor(resume=True).run([stage], tmp_context)
    cache = StageCache(tmp_context.datastore.root.parent / "cache" / "stage-cache")
    original = cache.fingerprint(stage, tmp_context, [])
    assert cache.is_current(stage, tmp_context, original)

    import canidae.core.stage_cache as stage_cache_module

    monkeypatch.setattr(stage_cache_module, "_package_fingerprint", lambda _path: "changed-code")
    changed = cache.fingerprint(stage, tmp_context, [])
    assert not cache.is_current(stage, tmp_context, changed)


def test_failed_stage_never_promotes_partial_outputs_or_index_entries(
    tmp_context: RunContext,
) -> None:
    report = NativeExecutor().run([_BrokenStage()], tmp_context)

    assert "broken_atomic" in report.failed
    assert not tmp_context.datastore.has(ArtifactKind.REPORT, "broken")
    assert not (tmp_context.datastore.root / "broken_atomic").exists()
    manifest = (Path(tmp_context.run_dir) / "manifest.json").read_text(encoding="utf-8")
    assert '"status": "failed"' in manifest
    assert "intentional failure" in manifest


def test_cancel_marker_stops_at_a_safe_boundary_with_recovery_state(
    tmp_context: RunContext,
) -> None:
    manager = ResourceManager(
        tmp_context.config.resource_manager, run_dir=Path(tmp_context.run_dir)
    )
    manager.request_cancel()

    report = NativeExecutor(resource_manager=manager).run([_CountingStage()], tmp_context)

    assert report.cancelled == ["cache_counter"]
    assert not tmp_context.datastore.has(ArtifactKind.ANALYSIS_RESULT, "cache_counter")
    manifest = (Path(tmp_context.run_dir) / "manifest.json").read_text(encoding="utf-8")
    assert "cancel_marker_present" in manifest
    assert "Remove" in manifest
