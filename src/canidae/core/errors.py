"""Typed exception hierarchy for CANIS.

Every error the platform raises deliberately derives from :class:`CanidaeError`, so callers
(the CLI, notebooks, test harnesses) can distinguish *expected* pipeline failures from
genuine bugs with a single ``except CanidaeError``.
"""

from __future__ import annotations


class CanidaeError(Exception):
    """Base class for all CANIS errors."""


class ConfigError(CanidaeError):
    """Configuration is missing, malformed, or fails validation."""


class RegistryError(CanidaeError):
    """A stage, dataset, reference, or analysis was requested but not registered
    (or was registered twice under the same name)."""


class DataStoreError(CanidaeError):
    """An artifact could not be registered, located, or retrieved."""


class IntegrityError(CanidaeError):
    """A checksum, size, or content-integrity check failed."""


class MissingToolError(CanidaeError):
    """A required external tool is absent or its version does not satisfy the constraint."""


class ExternalToolError(CanidaeError):
    """An external command exited non-zero or failed to produce its declared outputs."""

    def __init__(self, message: str, *, returncode: int | None = None,
                 stderr: str | None = None) -> None:
        super().__init__(message)
        self.returncode = returncode
        self.stderr = stderr


class StageInputError(CanidaeError):
    """A stage's declared inputs are unsatisfied or fail validation."""


class WorkflowError(CanidaeError):
    """The workflow graph is invalid (cycle, unsatisfiable dependency, duplicate producer)."""


class ResourcePreflightError(CanidaeError):
    """The local machine cannot safely start a run with the configured resource floor."""

    def __init__(self, message: str, *, recovery: str) -> None:
        super().__init__(message)
        self.recovery = recovery


class ResourceLimitError(CanidaeError):
    """A stage's declared memory requirement cannot fit inside the configured limit."""

    def __init__(self, stage: str, *, requested_mb: int, limit_mb: int) -> None:
        super().__init__(
            f"stage '{stage}' requests {requested_mb} MiB but the local memory limit is "
            f"{limit_mb} MiB"
        )
        self.stage = stage
        self.requested_mb = requested_mb
        self.limit_mb = limit_mb
        self.recovery = (
            "Increase resource_manager.max_memory_mb, lower the stage resource request, "
            "or use a smaller reduced-panel preset before resuming."
        )
