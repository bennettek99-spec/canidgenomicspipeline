"""Pipeline assembly: turn a validated config into runnable, wired stages.

This is the thin layer that binds the foundation together:

1. resolve the ordered ``config.pipeline`` names against the stage :data:`STAGES` registry,
2. validate each stage's config block against that stage's own model,
3. build the shared run context (datastore, runner, provenance, logging), and
4. hand the stage list to the executor.

The artifact datastore lives under ``data_root/store`` so it persists across invocations
(enabling resume); each invocation gets its own timestamped run directory for logs and the
provenance manifest.
"""

from __future__ import annotations

import random
from datetime import datetime
from pathlib import Path

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.executor import ExecutionReport, NativeExecutor, build_context
from canidae.core.logging import configure_logging, get_logger
from canidae.core.model import Cohort
from canidae.core.provenance import ProvenanceWriter
from canidae.core.registry import STAGES, load_entry_point_plugins
from canidae.core.resources import ResourceManager
from canidae.core.runtime import make_runner
from canidae.core.stage import Stage
from canidae.core.stage_cache import StageCache

_log = get_logger("pipeline")


def set_global_seed(seed: int) -> None:
    random.seed(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:  # pragma: no cover
        pass


def instantiate_stages(config: GlobalConfig) -> list[Stage]:
    """Build stage instances for the configured pipeline, validating each config block."""
    from canidae.stages import load_builtin_stages  # local import: optional heavy deps

    load_builtin_stages()
    load_entry_point_plugins()
    stages: list[Stage] = []
    for name in config.pipeline:
        stage_cls = STAGES.get(name)
        stage_config = config.parse_stage_config(name, stage_cls.config_model)
        stages.append(stage_cls(stage_config))
    return stages


def make_run_dir(config: GlobalConfig, run_id: str | None = None) -> Path:
    run_id = run_id or datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = config.paths.run_root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def run_pipeline(
    config: GlobalConfig,
    *,
    cohort: Cohort | None = None,
    dry_run: bool = False,
    run_id: str | None = None,
) -> ExecutionReport:
    """Execute the configured pipeline end-to-end."""
    run_dir = make_run_dir(config, run_id)
    configure_logging(config.logging, run_dir=run_dir)
    set_global_seed(config.seed)

    _log.info("project '%s' | config digest %s", config.project_name,
              config.digest()[:16])

    datastore = DataStore(config.paths.data_root / "store")
    runner = make_runner(config.containers)
    provenance = ProvenanceWriter(
        run_dir, config_digest=config.digest(), seed=config.seed,
        repo_dir=config.paths.root,
    )
    resource_manager = ResourceManager(config.resource_manager, run_dir=run_dir)
    try:
        disk_checks = resource_manager.preflight(
            [config.paths.data_root, config.paths.run_root, config.paths.cache_root]
        )
        applied_threads = resource_manager.apply_thread_limits()
        stale_transactions = resource_manager.cleanup_stale_transactions(datastore.root)
        provenance.add_run_metadata(
            resource_preflight=[
                {"path": str(check.path), "free_mb": check.free_mb,
                 "required_free_mb": check.required_free_mb}
                for check in disk_checks
            ],
            thread_limits=applied_threads,
            stale_transactions_removed=stale_transactions,
            pause_marker=str(resource_manager.pause_marker),
            cancel_marker=str(resource_manager.cancel_marker),
        )
    except Exception as exc:
        provenance.record_run_event(
            "resource_preflight", status="failed", error=f"{type(exc).__name__}: {exc}",
            recovery=getattr(
                exc, "recovery", "Resolve the local resource preflight error and rerun."
            ),
        )
        provenance.flush()
        raise
    ctx = build_context(
        config, datastore, runner, provenance, cohort=cohort, run_dir=run_dir,
        resource_manager=resource_manager,
    )

    stages = instantiate_stages(config)
    executor = NativeExecutor(
        max_workers=config.executor.max_workers,
        fail_fast=config.executor.fail_fast,
        resume=config.executor.resume,
        resource_manager=resource_manager,
        stage_cache=StageCache(config.paths.cache_root / "stage-cache"),
    )
    report = executor.run(stages, ctx, dry_run=dry_run)
    provenance.flush()

    _log.info("done: %d executed, %d skipped, %d failed",
              len(report.executed), len(report.skipped), len(report.failed))
    return report
