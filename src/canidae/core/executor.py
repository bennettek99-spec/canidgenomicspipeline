"""Native Python DAG executor.

The executor assembles a dependency graph from stages' *declared* inputs and outputs
(matched by ``(kind, role)``), topologically sorts it into levels of independent stages,
and runs each level — optionally in parallel — with resume/skip and dry-run support.

The heavy work inside a stage is typically a subprocess (an aligner, a caller), so
orchestration parallelism uses threads: it is I/O-bound and avoids pickling the datastore.
Outputs are registered back into the datastore *serially* between levels, so downstream
levels always see a consistent view. This is the only executor for now; a SLURM/cloud
backend can replace it behind the same :meth:`Executor.run` signature.
"""

from __future__ import annotations

from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path

from canidae.core.errors import DataStoreError, ResourceLimitError, StageInputError, WorkflowError
from canidae.core.logging import get_logger
from canidae.core.model import Artifact, Cohort
from canidae.core.provenance import ProvenanceRecord
from canidae.core.resources import ResourceManager
from canidae.core.stage import RunContext, Stage, StageResult
from canidae.core.stage_cache import StageCache
from canidae.core.staging import StagedDataStore, rewrite_staged_paths

_log = get_logger("executor")


@dataclass(slots=True)
class WorkflowGraph:
    """A resolved DAG of stages."""

    stages: dict[str, Stage]
    edges: dict[str, set[str]]          # producer -> {consumers}
    external_inputs: set[tuple[str, str]]  # (kind, role) not produced by any stage

    @classmethod
    def from_stages(cls, stages: list[Stage]) -> WorkflowGraph:
        by_name: dict[str, Stage] = {}
        producer: dict[tuple, str] = {}
        for stage in stages:
            if stage.name in by_name:
                raise WorkflowError(f"duplicate stage name in pipeline: '{stage.name}'")
            by_name[stage.name] = stage
            for spec in stage.produced_outputs():
                if spec.key in producer:
                    raise WorkflowError(
                        f"two stages produce {spec.kind.value}:{spec.role}: "
                        f"'{producer[spec.key]}' and '{stage.name}'"
                    )
                producer[spec.key] = stage.name

        edges: dict[str, set[str]] = {name: set() for name in by_name}
        external: set[tuple[str, str]] = set()
        for stage in stages:
            for spec in stage.required_inputs():
                src = producer.get(spec.key)
                if src is not None:
                    edges[src].add(stage.name)
                elif not spec.optional:
                    external.add((spec.kind.value, spec.role))
        return cls(stages=by_name, edges=edges, external_inputs=external)

    def topological_levels(self) -> list[list[str]]:
        """Kahn's algorithm, grouped into levels of mutually-independent stages."""
        indeg: dict[str, int] = defaultdict(int)
        for consumers in self.edges.values():
            for c in consumers:
                indeg[c] += 1
        ready = deque(sorted(n for n in self.stages if indeg[n] == 0))
        levels: list[list[str]] = []
        seen = 0
        while ready:
            level = sorted(ready)
            ready.clear()
            levels.append(level)
            for name in level:
                seen += 1
                for consumer in sorted(self.edges[name]):
                    indeg[consumer] -= 1
                    if indeg[consumer] == 0:
                        ready.append(consumer)
        if seen != len(self.stages):
            raise WorkflowError("workflow graph contains a cycle")
        return levels


@dataclass(slots=True)
class ExecutionReport:
    """Summary of a run."""

    executed: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    failed: dict[str, str] = field(default_factory=dict)
    paused: list[str] = field(default_factory=list)
    cancelled: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.failed


@dataclass(slots=True)
class _StageExecution:
    stage: Stage
    result: StageResult
    record: ProvenanceRecord
    datastore: StagedDataStore


class NativeExecutor:
    """Runs a DAG locally with fingerprinted resume and transactional outputs."""

    def __init__(
        self,
        max_workers: int = 4,
        *,
        fail_fast: bool = True,
        resume: bool = True,
        resource_manager: ResourceManager | None = None,
        stage_cache: StageCache | None = None,
    ) -> None:
        self.max_workers = max_workers
        self.fail_fast = fail_fast
        self.resume = resume
        self.resource_manager = resource_manager
        self.stage_cache = stage_cache

    def plan(self, stages: list[Stage]) -> list[list[str]]:
        return WorkflowGraph.from_stages(stages).topological_levels()

    def run(
        self,
        stages: list[Stage],
        ctx: RunContext,
        *,
        dry_run: bool = False,
    ) -> ExecutionReport:
        graph = WorkflowGraph.from_stages(stages)
        levels = graph.topological_levels()
        report = ExecutionReport()
        manager = self._resource_manager_for(ctx)
        cache = self._cache_for(ctx)

        _log.info("workflow: %d stages in %d levels", len(stages), len(levels))
        if graph.external_inputs:
            need = ", ".join(f"{k}:{r}" for k, r in sorted(graph.external_inputs))
            _log.info("expects pre-seeded external inputs: %s", need)

        if dry_run:
            for i, level in enumerate(levels):
                _log.info("level %d: %s", i, ", ".join(level))
            report.executed = [n for level in levels for n in level]
            return report

        stop = False
        for level_index, level in enumerate(levels):
            if manager.cancel_requested():
                self._record_cancel(graph, levels, level_index, report, ctx, manager)
                break
            if manager.pause_requested():
                self._record_pause(graph, levels, level_index, report, ctx, manager)
                break

            stages_to_run: list[Stage] = []
            fingerprints: dict[str, dict] = {}
            for name in level:
                stage = graph.stages[name]
                if not stage.config.enabled:
                    self._record_skip(
                        stage,
                        ctx,
                        report,
                        reason="disabled_by_config",
                        recovery=f"Set stages.{stage.name}.enabled=true to include this stage.",
                    )
                    continue
                if self.resume:
                    try:
                        stage.validate_inputs(ctx)
                        inputs = stage.gather_inputs(ctx)
                        fingerprint = cache.fingerprint(stage, ctx, inputs)
                    except Exception as exc:
                        self._record_failure(stage, ctx, report, exc)
                        if self.fail_fast:
                            stop = True
                            break
                        continue
                    if cache.is_current(stage, ctx, fingerprint):
                        _log.info("stage '%s' fingerprint unchanged; safe-resume skip", name)
                        record = ctx.provenance.start(stage.name)
                        ctx.provenance.record_inputs(record, inputs)
                        ctx.provenance.finish(
                            record,
                            status="skipped",
                            skip_reason="cache_hit",
                            recovery=(
                                f"Delete {cache.path_for(stage.name)} or set executor.resume=false "
                                "to force a recalculation."
                            ),
                        )
                        report.skipped.append(stage.name)
                        continue
                    fingerprints[stage.name] = fingerprint
                stages_to_run.append(stage)

            if stop:
                break
            if not stages_to_run:
                continue

            batches, limit_stop = self._batches_with_limits(
                stages_to_run, ctx, report, manager
            )
            if limit_stop:
                if self.fail_fast:
                    break
                if not batches:
                    continue

            for batch_index, batch in enumerate(batches):
                if manager.cancel_requested():
                    remaining_names = [
                        stage.name for remaining in batches[batch_index:] for stage in remaining
                    ]
                    self._record_cancel(
                        graph,
                        levels,
                        level_index,
                        report,
                        ctx,
                        manager,
                        current_remaining=remaining_names,
                    )
                    stop = True
                    break
                if manager.pause_requested():
                    remaining_names = [
                        stage.name for remaining in batches[batch_index:] for stage in remaining
                    ]
                    self._record_pause(
                        graph,
                        levels,
                        level_index,
                        report,
                        ctx,
                        manager,
                        current_remaining=remaining_names,
                    )
                    stop = True
                    break

                executions = self._run_batch(batch, ctx, report, manager)
                for execution in executions:
                    stage = execution.stage
                    result = execution.result
                    record = execution.record
                    if result.skipped:
                        execution.datastore.cleanup()
                        ctx.provenance.finish(
                            record,
                            metrics=result.metrics,
                            status="skipped",
                            skip_reason="stage_reported_skip",
                            recovery=(
                                "Review the stage metrics, then rerun after changing its "
                                "inputs if needed."
                            ),
                        )
                        report.skipped.append(stage.name)
                        continue
                    try:
                        stamped = self._register_outputs(
                            stage, result, record, ctx, execution.datastore
                        )
                    except Exception as exc:
                        execution.datastore.cleanup()
                        self._record_failure(stage, ctx, report, exc, record=record)
                        continue

                    result.metrics = rewrite_staged_paths(
                        result.metrics,
                        execution.datastore.stage_root,
                        ctx.datastore.root / stage.name,
                    )

                    result.metrics.setdefault(
                        "cache_fingerprint",
                        fingerprints.get(stage.name, {}).get("fingerprint"),
                    )
                    try:
                        fingerprint = fingerprints.get(stage.name)
                        if fingerprint is None:
                            fingerprint = cache.fingerprint(stage, ctx, stage.gather_inputs(ctx))
                        cache_path = cache.store(stage, fingerprint, stamped)
                        result.metrics["cache_metadata"] = str(cache_path)
                    # Output correctness takes priority over cache convenience.
                    except Exception as exc:
                        _log.warning(
                            "could not persist cache metadata for stage '%s': %s", stage.name, exc
                        )
                        result.metrics["cache_metadata_error"] = f"{type(exc).__name__}: {exc}"

                    ctx.provenance.finish(
                        record,
                        outputs=stamped,
                        metrics=result.metrics,
                        status="succeeded",
                    )
                    report.executed.append(stage.name)

                if self.fail_fast and report.failed:
                    stop = True
                    break
            if stop:
                break

        return report

    # -- scheduling and lifecycle helpers -------------------------------------------

    def _resource_manager_for(self, ctx: RunContext) -> ResourceManager:
        if self.resource_manager is not None:
            return self.resource_manager
        if ctx.resource_manager is not None:
            return ctx.resource_manager
        run_dir = Path(ctx.run_dir) if ctx.run_dir is not None else ctx.datastore.root
        return ResourceManager(ctx.config.resource_manager, run_dir=run_dir)

    def _cache_for(self, ctx: RunContext) -> StageCache:
        if self.stage_cache is not None:
            return self.stage_cache
        existing = ctx.scratch.get("_stage_cache")
        if isinstance(existing, StageCache):
            return existing
        # Cache metadata must outlive a timestamped run directory.  Large artifacts remain
        # in the datastore; this persistent folder contains only compact fingerprints.
        cache_root = Path(ctx.config.paths.cache_root)
        if not cache_root.is_absolute():
            cache_root = ctx.datastore.root.parent / "cache"
        root = cache_root / "stage-cache"
        cache = StageCache(root)
        ctx.scratch["_stage_cache"] = cache
        return cache

    def _batches_with_limits(
        self,
        stages: list[Stage],
        ctx: RunContext,
        report: ExecutionReport,
        manager: ResourceManager,
    ) -> tuple[list[list[Stage]], bool]:
        remaining = list(stages)
        saw_limit = False
        while remaining:
            try:
                batches = manager.batches(
                    remaining, ctx, requested_workers=self.max_workers
                )
                return batches, saw_limit
            except ResourceLimitError as exc:
                saw_limit = True
                stage = next(stage for stage in remaining if stage.name == exc.stage)
                self._record_failure(stage, ctx, report, exc)
                remaining.remove(stage)
                if self.fail_fast:
                    return [], True
        return [], saw_limit

    def _run_batch(
        self,
        stages: list[Stage],
        ctx: RunContext,
        report: ExecutionReport,
        manager: ResourceManager,
    ) -> list[_StageExecution]:
        results: list[_StageExecution] = []
        workers = max(1, min(manager.effective_workers(self.max_workers), len(stages)))
        records = {stage.name: ctx.provenance.start(stage.name) for stage in stages}

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(self._execute_one, stage, ctx, records[stage.name]): stage
                for stage in stages
            }
            for future in futures:
                stage = futures[future]
                try:
                    result, staged = future.result()
                    results.append(_StageExecution(stage, result, records[stage.name], staged))
                except Exception as exc:
                    _log.exception("stage '%s' failed", stage.name)
                    self._record_failure(stage, ctx, report, exc, record=records[stage.name])
        return results

    def _execute_one(
        self, stage: Stage, ctx: RunContext, record: ProvenanceRecord
    ) -> tuple[StageResult, StagedDataStore]:
        staged = StagedDataStore(ctx.datastore, stage.name)
        try:
            stage.validate_inputs(ctx)
            inputs = stage.gather_inputs(ctx)
            ctx.provenance.record_inputs(record, inputs)
            stage_ctx = replace(
                ctx,
                datastore=staged,
                scratch={**ctx.scratch, "_record": record},
            )
            result = stage.run(stage_ctx)
            result.provenance = record
            return result, staged
        except Exception:
            staged.cleanup()
            raise

    def _register_outputs(
        self,
        stage: Stage,
        result: StageResult,
        record: ProvenanceRecord,
        ctx: RunContext,
        staged: StagedDataStore,
    ) -> list[Artifact]:
        handles = [self._stamp(artifact, stage.name, record.id) for artifact in result.artifacts]
        promoted = staged.promote(handles)
        stamped = ctx.datastore.register_many(promoted)
        missing = [
            spec for spec in stage.produced_outputs()
            if not ctx.datastore.has(spec.kind, spec.role)
        ]
        if missing:
            pretty = ", ".join(f"{spec.kind.value}:{spec.role}" for spec in missing)
            raise DataStoreError(f"stage '{stage.name}' did not return declared outputs: {pretty}")
        result.artifacts = stamped
        return stamped

    def _record_skip(
        self,
        stage: Stage,
        ctx: RunContext,
        report: ExecutionReport,
        *,
        reason: str,
        recovery: str,
    ) -> None:
        record = ctx.provenance.start(stage.name)
        ctx.provenance.finish(record, status="skipped", skip_reason=reason, recovery=recovery)
        report.skipped.append(stage.name)

    def _record_failure(
        self,
        stage: Stage,
        ctx: RunContext,
        report: ExecutionReport,
        exc: Exception,
        *,
        record: ProvenanceRecord | None = None,
    ) -> None:
        error = f"{type(exc).__name__}: {exc}"
        report.failed[stage.name] = error
        current = record or ctx.provenance.start(stage.name)
        ctx.provenance.finish(
            current,
            status="failed",
            error=error,
            recovery=getattr(
                exc,
                "recovery",
                "Inspect the error, correct the input or configuration, then rerun safely.",
            ),
        )

    def _record_pause(
        self,
        graph: WorkflowGraph,
        levels: list[list[str]],
        level_index: int,
        report: ExecutionReport,
        ctx: RunContext,
        manager: ResourceManager,
        *,
        current_remaining: list[str] | None = None,
    ) -> None:
        already = (
            set(report.executed) | set(report.skipped) | set(report.failed) | set(report.paused)
        )
        names = current_remaining if current_remaining is not None else list(levels[level_index])
        names += [name for level in levels[level_index + 1:] for name in level]
        for name in names:
            if name in already:
                continue
            stage = graph.stages[name]
            if not stage.config.enabled:
                continue
            self._record_skip(
                stage,
                ctx,
                report,
                reason="pause_marker_present",
                recovery=(
                    f"Remove {manager.pause_marker} and rerun with the same run id to resume "
                    "from the last completed cached stage."
                ),
            )
            report.paused.append(name)

    def _record_cancel(
        self,
        graph: WorkflowGraph,
        levels: list[list[str]],
        level_index: int,
        report: ExecutionReport,
        ctx: RunContext,
        manager: ResourceManager,
        *,
        current_remaining: list[str] | None = None,
    ) -> None:
        """Record a cooperative stop at a batch boundary with resume instructions."""
        already = (
            set(report.executed) | set(report.skipped) | set(report.failed)
            | set(report.paused) | set(report.cancelled)
        )
        names = current_remaining if current_remaining is not None else list(levels[level_index])
        names += [name for level in levels[level_index + 1:] for name in level]
        for name in names:
            if name in already:
                continue
            stage = graph.stages[name]
            if not stage.config.enabled:
                continue
            self._record_skip(
                stage,
                ctx,
                report,
                reason="cancel_marker_present",
                recovery=(
                    f"Remove {manager.cancel_marker} and rerun with the same configuration; "
                    "unchanged completed stages will be recovered from the safe cache."
                ),
            )
            report.cancelled.append(name)

    @staticmethod
    def _stamp(artifact: Artifact, stage_name: str, prov_id: str | None) -> Artifact:
        return Artifact(
            kind=artifact.kind,
            role=artifact.role,
            path=artifact.path,
            fmt=artifact.fmt,
            checksum=artifact.checksum,
            produced_by=artifact.produced_by or stage_name,
            provenance_id=artifact.provenance_id or prov_id,
            schema_version=artifact.schema_version,
            metadata=artifact.metadata,
        )


def build_context(
    config,
    datastore,
    runner,
    provenance,
    *,
    cohort: Cohort | None = None,
    run_dir: Path | None = None,
    resource_manager: ResourceManager | None = None,
) -> RunContext:
    """Assemble the shared :class:`RunContext` for a run."""
    return RunContext(
        config=config, datastore=datastore, runner=runner, provenance=provenance,
        cohort=cohort, run_dir=run_dir, resource_manager=resource_manager,
    )


# Re-exported for callers that catch input errors distinctly.
__all__ = [
    "ExecutionReport",
    "NativeExecutor",
    "StageInputError",
    "WorkflowGraph",
    "build_context",
]
