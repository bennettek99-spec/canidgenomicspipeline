"""Laptop-oriented resource guards for the native executor."""

from __future__ import annotations

import os
import shutil
import time
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from canidae.core.config import ResourceManagerConfig
from canidae.core.errors import ResourceLimitError, ResourcePreflightError
from canidae.core.runtime import ResourceSpec

if TYPE_CHECKING:
    from canidae.core.stage import RunContext, Stage


@dataclass(frozen=True, slots=True)
class DiskSpaceCheck:
    path: Path
    free_mb: int
    required_free_mb: int


class ResourceManager:
    """Apply conservative local execution limits and recovery-friendly controls."""

    def __init__(self, config: ResourceManagerConfig, *, run_dir: Path) -> None:
        self.config = config
        self.run_dir = Path(run_dir)

    @property
    def pause_marker(self) -> Path:
        configured = self.config.pause_marker
        if configured is None:
            return self.run_dir / ".canidae.pause"
        return configured if configured.is_absolute() else self.run_dir / configured

    def pause_requested(self) -> bool:
        return self.pause_marker.exists()

    @property
    def cancel_marker(self) -> Path:
        configured = self.config.cancel_marker
        if configured is None:
            return self.run_dir / ".canidae.cancel"
        return configured if configured.is_absolute() else self.run_dir / configured

    def cancel_requested(self) -> bool:
        return self.cancel_marker.exists()

    def request_pause(self) -> Path:
        """Ask an active executor to pause at its next safe stage boundary."""
        self.pause_marker.parent.mkdir(parents=True, exist_ok=True)
        self.pause_marker.touch()
        return self.pause_marker

    def resume(self) -> None:
        """Clear a pause request; rerun with the same run id to continue from cache."""
        self.pause_marker.unlink(missing_ok=True)

    def request_cancel(self) -> Path:
        """Ask an active executor to stop at its next safe stage boundary."""
        self.cancel_marker.parent.mkdir(parents=True, exist_ok=True)
        self.cancel_marker.touch()
        return self.cancel_marker

    def clear_cancel(self) -> None:
        self.cancel_marker.unlink(missing_ok=True)

    def effective_workers(self, requested: int) -> int:
        limit = self.config.max_workers
        return max(1, min(requested, limit)) if limit else max(1, requested)

    def stage_resources(self, stage: Stage, ctx: RunContext) -> ResourceSpec:
        requested = stage.resources() or ResourceSpec.from_defaults(ctx.config.resources)
        if self.config.temperature_friendly:
            return ResourceSpec(
                cpus=min(requested.cpus, self.config.max_threads_per_stage),
                mem_mb=requested.mem_mb,
                gpus=requested.gpus,
                time_min=requested.time_min,
            )
        return requested

    def batches(
        self, stages: list[Stage], ctx: RunContext, *, requested_workers: int
    ) -> list[list[Stage]]:
        """Partition independent stages so declared memory never exceeds the local cap."""
        worker_limit = self.effective_workers(requested_workers)
        memory_limit = self.config.max_memory_mb
        batches: list[list[Stage]] = []
        current: list[Stage] = []
        current_memory = 0
        for stage in stages:
            request = self.stage_resources(stage, ctx)
            if memory_limit and request.mem_mb > memory_limit:
                raise ResourceLimitError(
                    stage.name,
                    requested_mb=request.mem_mb,
                    limit_mb=memory_limit,
                )
            exceeds_workers = len(current) >= worker_limit
            exceeds_memory = bool(
                memory_limit and current and current_memory + request.mem_mb > memory_limit
            )
            if exceeds_workers or exceeds_memory:
                batches.append(current)
                current = []
                current_memory = 0
            current.append(stage)
            current_memory += request.mem_mb
        if current:
            batches.append(current)
        return batches

    def preflight(self, paths: Iterable[Path]) -> list[DiskSpaceCheck]:
        """Check the filesystem(s) that will host run data before work begins."""
        checks: list[DiskSpaceCheck] = []
        seen: set[tuple[str, int]] = set()
        for raw_path in paths:
            path = _existing_parent(Path(raw_path))
            usage = shutil.disk_usage(path)
            # ``st_dev`` avoids duplicate checks when multiple run paths share one volume.
            identity = (str(path.drive).lower(), path.stat().st_dev)
            if identity in seen:
                continue
            seen.add(identity)
            check = DiskSpaceCheck(
                path=path,
                free_mb=usage.free // (1024 * 1024),
                required_free_mb=self.config.min_free_disk_mb,
            )
            checks.append(check)
            if check.free_mb < check.required_free_mb:
                raise ResourcePreflightError(
                    f"only {check.free_mb} MiB free at {path}; the configured safety floor is "
                    f"{check.required_free_mb} MiB",
                    recovery=(
                        "Free disk space, lower resource_manager.min_free_disk_mb only if "
                        "you have independently checked the run's projected size, then rerun."
                    ),
                )
        return checks

    def apply_thread_limits(self) -> dict[str, str]:
        """Set conservative numerical-library defaults without overriding user choices."""
        if not self.config.temperature_friendly:
            return {}
        value = str(self.config.max_threads_per_stage)
        applied: dict[str, str] = {}
        for name in (
            "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"
        ):
            if name not in os.environ:
                os.environ[name] = value
                applied[name] = value
        return applied

    def cleanup_stale_transactions(self, datastore_root: Path) -> int:
        """Remove abandoned transaction directories old enough not to be an active run."""
        if not self.config.cleanup_intermediates:
            return 0
        staging_root = Path(datastore_root) / ".staging"
        if not staging_root.exists():
            return 0
        oldest = time.time() - (self.config.stale_transaction_hours * 3600)
        removed = 0
        for child in staging_root.iterdir():
            try:
                if child.is_dir() and child.stat().st_mtime < oldest:
                    shutil.rmtree(child, ignore_errors=True)
                    removed += 1
            except OSError:
                continue
        return removed


def _existing_parent(path: Path) -> Path:
    candidate = path.resolve()
    while not candidate.exists() and candidate.parent != candidate:
        candidate = candidate.parent
    return candidate
