"""Configuration system for CANIS.

Every run is fully described by a layered configuration:

    packaged defaults  <  project config  <  extra YAML files  <  dotted CLI overrides

All layers are deep-merged, then validated by Pydantic. The resulting :class:`GlobalConfig`
is immutable and hashable, so its digest can be recorded in provenance and used for
cache invalidation. No analysis parameter is ever hard-coded in a stage; stages read a
validated slice of this config.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Literal, TypeVar

import yaml
from pydantic import BaseModel, ConfigDict, Field

from canidae.core.errors import ConfigError

_T = TypeVar("_T", bound=BaseModel)

_PACKAGE_DEFAULTS = Path(__file__).resolve().parents[3] / "configs" / "defaults.yaml"


class _Frozen(BaseModel):
    """Base for immutable, extra-forbidding config models."""

    model_config = ConfigDict(frozen=True, extra="forbid")


class PathsConfig(_Frozen):
    """Filesystem layout. Relative paths are resolved against ``root``."""

    root: Path = Path(".")
    data_root: Path = Path("data")
    run_root: Path = Path("runs")
    reference_root: Path = Path("data/references")
    cache_root: Path = Path("data/cache")

    def resolved(self) -> PathsConfig:
        base = self.root.resolve()

        def r(p: Path) -> Path:
            return p if p.is_absolute() else base / p

        return PathsConfig(
            root=base,
            data_root=r(self.data_root),
            run_root=r(self.run_root),
            reference_root=r(self.reference_root),
            cache_root=r(self.cache_root),
        )


class ResourceDefaults(_Frozen):
    """Default resource request applied to a stage that does not specify its own."""

    cpus: int = Field(default=1, ge=1)
    mem_mb: int = Field(default=4096, ge=256)
    gpus: int = Field(default=0, ge=0)
    time_min: int = Field(default=60, ge=1)


class ExecutorConfig(_Frozen):
    backend: Literal["native"] = "native"
    max_workers: int = Field(default=4, ge=1)
    fail_fast: bool = True
    resume: bool = True  # skip only when outputs and their safe-resume fingerprint agree


class ResourceManagerConfig(_Frozen):
    """Local-machine guardrails layered on top of per-stage resource requests.

    ``0`` disables a numerical cap where that is useful on shared workstations.  The
    laptop profile supplies concrete limits; packaged defaults remain portable.
    """

    max_memory_mb: int = Field(default=0, ge=0)
    max_workers: int = Field(default=0, ge=0)
    min_free_disk_mb: int = Field(default=2048, ge=0)
    pause_marker: Path | None = None
    cancel_marker: Path | None = None
    cleanup_intermediates: bool = True
    stale_transaction_hours: int = Field(default=24, ge=0)
    temperature_friendly: bool = True
    max_threads_per_stage: int = Field(default=4, ge=1)
    process_chromosomes_sequentially: bool = True


class ContainerConfig(_Frozen):
    engine: Literal["none", "apptainer", "docker"] = "none"
    image_dir: Path = Path("containers")
    bind_paths: tuple[str, ...] = ()


class LoggingConfig(_Frozen):
    level: Literal["DEBUG", "INFO", "WARNING", "ERROR"] = "INFO"
    json_format: bool = False
    rich_console: bool = True


class GlobalConfig(_Frozen):
    """The fully-resolved, validated configuration for a run."""

    project_name: str = "canis"
    seed: int = 1234
    paths: PathsConfig = PathsConfig()
    executor: ExecutorConfig = ExecutorConfig()
    resources: ResourceDefaults = ResourceDefaults()
    resource_manager: ResourceManagerConfig = ResourceManagerConfig()
    containers: ContainerConfig = ContainerConfig()
    logging: LoggingConfig = LoggingConfig()

    # Pipeline recipe: ordered stage names to run (resolved against the registry).
    pipeline: tuple[str, ...] = ()

    # Per-stage config blocks, validated lazily by each stage against its own model.
    stages: dict[str, dict[str, Any]] = Field(default_factory=dict)

    # Optional pinned tool versions / binary path overrides (tool name -> spec).
    tools: dict[str, str] = Field(default_factory=dict)

    # -- construction ------------------------------------------------------------------

    @classmethod
    def load(
        cls,
        *files: str | Path,
        overrides: Mapping[str, Any] | None = None,
        include_defaults: bool = True,
    ) -> GlobalConfig:
        """Load and deep-merge config layers, then validate.

        ``files`` are merged left-to-right on top of the packaged defaults. ``overrides``
        is a mapping of dotted keys (``executor.max_workers``) to values, applied last.
        """
        merged: dict[str, Any] = {}
        layers: list[Path] = []
        if include_defaults and _PACKAGE_DEFAULTS.exists():
            layers.append(_PACKAGE_DEFAULTS)
        layers.extend(Path(f) for f in files)

        for layer in layers:
            if not layer.exists():
                raise ConfigError(f"config file not found: {layer}")
            try:
                data = yaml.safe_load(layer.read_text(encoding="utf-8")) or {}
            except yaml.YAMLError as exc:  # pragma: no cover - passthrough
                raise ConfigError(f"invalid YAML in {layer}: {exc}") from exc
            if not isinstance(data, Mapping):
                raise ConfigError(f"top-level YAML in {layer} must be a mapping")
            _deep_merge(merged, dict(data))

        for dotted, value in (overrides or {}).items():
            _assign_dotted(merged, dotted, value)

        try:
            cfg = cls.model_validate(merged)
        except Exception as exc:  # pydantic ValidationError
            raise ConfigError(f"configuration failed validation: {exc}") from exc
        return cfg.with_resolved_paths()

    def with_resolved_paths(self) -> GlobalConfig:
        return self.model_copy(update={"paths": self.paths.resolved()})

    # -- stage access ------------------------------------------------------------------

    def stage_config(self, stage_name: str) -> dict[str, Any]:
        """Return the raw config block for a stage (empty dict if absent)."""
        return dict(self.stages.get(stage_name, {}))

    def parse_stage_config(self, stage_name: str, model: type[_T]) -> _T:
        """Validate a stage's config block against its own Pydantic model."""
        try:
            return model.model_validate(self.stage_config(stage_name))
        except Exception as exc:
            raise ConfigError(f"config for stage '{stage_name}' failed validation: {exc}") from exc

    # -- provenance --------------------------------------------------------------------

    def digest(self) -> str:
        """A stable SHA-256 over the canonicalized config (used for cache invalidation)."""
        payload = json.dumps(self.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    def to_yaml(self) -> str:
        return yaml.safe_dump(self.model_dump(mode="json"), sort_keys=True)


# --------------------------------------------------------------------------------------
# Merge helpers
# --------------------------------------------------------------------------------------


def _deep_merge(base: dict[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge ``overlay`` into ``base`` in place. Mappings merge; scalars and
    lists replace."""
    for key, value in overlay.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, Mapping):
            _deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def _assign_dotted(target: dict[str, Any], dotted: str, value: Any) -> None:
    """Assign ``value`` into ``target`` at a dotted path, creating intermediate dicts."""
    parts = dotted.split(".")
    node = target
    for part in parts[:-1]:
        nxt = node.get(part)
        if not isinstance(nxt, dict):
            nxt = {}
            node[part] = nxt
        node = nxt
    node[parts[-1]] = _coerce_scalar(value)


def _coerce_scalar(value: Any) -> Any:
    """Best-effort coercion of CLI string overrides to int/float/bool/None."""
    if not isinstance(value, str):
        return value
    lowered = value.strip().lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    if lowered in {"null", "none"}:
        return None
    for caster in (int, float):
        try:
            return caster(value)
        except ValueError:
            continue
    return value
