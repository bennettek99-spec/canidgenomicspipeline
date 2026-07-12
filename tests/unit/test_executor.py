from __future__ import annotations

from canidae.core.errors import WorkflowError
from canidae.core.executor import NativeExecutor, WorkflowGraph
from canidae.core.model import ArtifactKind
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult


def test_graph_orders_producer_before_consumer(producer, consumer) -> None:
    levels = WorkflowGraph.from_stages([consumer, producer]).topological_levels()
    order = [n for level in levels for n in level]
    assert order.index("producer") < order.index("consumer")


def test_external_inputs_detected(consumer) -> None:
    graph = WorkflowGraph.from_stages([consumer])
    assert ("qc_table", "metrics") in graph.external_inputs


def test_duplicate_producer_rejected(producer) -> None:
    class Dup(Stage):
        name = "dup"
        config_model = StageConfig

        def required_inputs(self):
            return []

        def produced_outputs(self):
            return [ArtifactSpec(ArtifactKind.QC_TABLE, "metrics")]

        def run(self, ctx):
            return StageResult()

    try:
        WorkflowGraph.from_stages([producer, Dup()])
    except WorkflowError:
        return
    raise AssertionError("expected WorkflowError for duplicate producer")


def test_cycle_detected() -> None:
    class A(Stage):
        name = "a"
        config_model = StageConfig

        def required_inputs(self):
            return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "b_out")]

        def produced_outputs(self):
            return [ArtifactSpec(ArtifactKind.QC_TABLE, "a_out")]

        def run(self, ctx):
            return StageResult()

    class B(Stage):
        name = "b"
        config_model = StageConfig

        def required_inputs(self):
            return [ArtifactSpec(ArtifactKind.QC_TABLE, "a_out")]

        def produced_outputs(self):
            return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "b_out")]

        def run(self, ctx):
            return StageResult()

    graph = WorkflowGraph.from_stages([A(), B()])
    try:
        graph.topological_levels()
    except WorkflowError:
        return
    raise AssertionError("expected WorkflowError for cycle")


def test_end_to_end_execution(tmp_context: RunContext, producer, consumer) -> None:
    executor = NativeExecutor(max_workers=2)
    report = executor.run([consumer, producer], tmp_context)
    assert report.ok
    assert set(report.executed) == {"producer", "consumer"}
    assert tmp_context.datastore.has(ArtifactKind.ANALYSIS_RESULT, "summary")
    # provenance manifest was written
    assert (tmp_context.run_dir / "manifest.json").exists()


def test_resume_skips_completed(tmp_context: RunContext, producer, consumer) -> None:
    executor = NativeExecutor(max_workers=2, resume=True)
    executor.run([producer, consumer], tmp_context)
    second = executor.run([producer, consumer], tmp_context)
    assert set(second.skipped) == {"producer", "consumer"}
    assert not second.executed


def test_dry_run_does_not_execute(tmp_context: RunContext, producer, consumer) -> None:
    executor = NativeExecutor()
    report = executor.run([producer, consumer], tmp_context, dry_run=True)
    assert report.executed == ["producer", "consumer"]
    # nothing actually produced
    assert not tmp_context.datastore.has(ArtifactKind.QC_TABLE, "metrics")


def test_failing_stage_is_reported(tmp_context: RunContext) -> None:
    class Boom(Stage):
        name = "boom"
        config_model = StageConfig

        def required_inputs(self):
            return []

        def produced_outputs(self):
            return [ArtifactSpec(ArtifactKind.REPORT, "r")]

        def run(self, ctx):
            raise RuntimeError("kaboom")

    report = NativeExecutor().run([Boom()], tmp_context)
    assert not report.ok
    assert "boom" in report.failed
