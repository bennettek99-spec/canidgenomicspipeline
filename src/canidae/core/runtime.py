"""Runtime layer: a uniform way to invoke external bioinformatics tools.

Stages never call :func:`subprocess.run` directly. Instead they describe a tool with a
:class:`ToolSpec` (binary name + how to read its version + a version constraint) and a
resource request with :class:`ResourceSpec`, then hand a command to a :class:`ToolRunner`.

Runners come in flavours behind one interface:

* :class:`LocalRunner`     — run on the local host.
* :class:`ContainerRunner` — run inside Apptainer/Docker for reproducibility.

(A cluster/SLURM runner and a dry-run runner slot in behind the same ``run`` signature
without touching any stage.) The runner records a :class:`ToolInvocation` on the active
provenance record, so every external command is captured automatically.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

from canidae.core.config import ContainerConfig, ResourceDefaults
from canidae.core.errors import ExternalToolError, MissingToolError
from canidae.core.logging import get_logger
from canidae.core.provenance import ProvenanceRecord, ToolInvocation

_log = get_logger("runtime")


@dataclass(frozen=True, slots=True)
class ResourceSpec:
    """A resource request for a single command."""

    cpus: int = 1
    mem_mb: int = 4096
    gpus: int = 0
    time_min: int = 60

    @classmethod
    def from_defaults(cls, defaults: ResourceDefaults) -> ResourceSpec:
        return cls(defaults.cpus, defaults.mem_mb, defaults.gpus, defaults.time_min)


@dataclass(frozen=True, slots=True)
class ToolSpec:
    """Declares an external tool and how to check its version."""

    name: str
    binary: str | None = None  # defaults to ``name``
    version_args: tuple[str, ...] = ("--version",)
    version_regex: str = r"(\d+\.\d+(?:\.\d+)?)"
    min_version: str | None = None
    container_image: str | None = None  # image name/URI when containerized

    @property
    def executable(self) -> str:
        return self.binary or self.name


@dataclass(slots=True)
class CommandResult:
    """The outcome of running a command."""

    argv: list[str]
    returncode: int
    stdout: str
    stderr: str
    duration_s: float

    @property
    def ok(self) -> bool:
        return self.returncode == 0


@runtime_checkable
class ToolRunner(Protocol):
    """Interface every runner satisfies."""

    def which(self, spec: ToolSpec) -> str | None: ...

    def version(self, spec: ToolSpec) -> str | None: ...

    def ensure(self, spec: ToolSpec) -> None: ...

    def run(
        self,
        spec: ToolSpec,
        args: list[str],
        *,
        resources: ResourceSpec | None = None,
        expect_outputs: list[Path] | None = None,
        record: ProvenanceRecord | None = None,
        cwd: Path | None = None,
        check: bool = True,
        dry_run: bool = False,
    ) -> CommandResult: ...


class _BaseRunner:
    """Shared version-checking, output-validation, and provenance-recording logic."""

    def which(self, spec: ToolSpec) -> str | None:
        return shutil.which(spec.executable)

    def version(self, spec: ToolSpec) -> str | None:
        path = self.which(spec)
        if path is None:
            return None
        try:
            out = subprocess.run(
                [path, *spec.version_args], capture_output=True, text=True, timeout=15
            )
        except (OSError, subprocess.SubprocessError):
            return None
        text = f"{out.stdout}\n{out.stderr}"
        match = re.search(spec.version_regex, text)
        return match.group(1) if match else None

    def ensure(self, spec: ToolSpec) -> None:
        if self.which(spec) is None:
            raise MissingToolError(
                f"required tool '{spec.name}' (binary '{spec.executable}') not found on PATH"
            )
        if spec.min_version:
            found = self.version(spec)
            if found and _version_tuple(found) < _version_tuple(spec.min_version):
                raise MissingToolError(
                    f"tool '{spec.name}' version {found} < required {spec.min_version}"
                )

    def _finalize(
        self,
        spec: ToolSpec,
        argv: list[str],
        proc: subprocess.CompletedProcess[str],
        started: float,
        expect_outputs: list[Path] | None,
        record: ProvenanceRecord | None,
        check: bool,
    ) -> CommandResult:
        duration = time.monotonic() - started
        result = CommandResult(argv, proc.returncode, proc.stdout, proc.stderr, duration)
        if record is not None:
            record.add_tool(
                ToolInvocation(
                    tool=spec.name,
                    version=self.version(spec),
                    argv=argv,
                    returncode=proc.returncode,
                    duration_s=round(duration, 3),
                )
            )
        if check and not result.ok:
            raise ExternalToolError(
                f"tool '{spec.name}' exited {proc.returncode}",
                returncode=proc.returncode,
                stderr=proc.stderr[-4000:],
            )
        if check and expect_outputs:
            missing = [p for p in expect_outputs if not Path(p).exists()]
            if missing:
                raise ExternalToolError(
                    f"tool '{spec.name}' did not produce expected outputs: {missing}"
                )
        return result


class LocalRunner(_BaseRunner):
    """Runs commands directly on the local host."""

    def run(
        self,
        spec: ToolSpec,
        args: list[str],
        *,
        resources: ResourceSpec | None = None,
        expect_outputs: list[Path] | None = None,
        record: ProvenanceRecord | None = None,
        cwd: Path | None = None,
        check: bool = True,
        dry_run: bool = False,
    ) -> CommandResult:
        argv = [spec.executable, *args]
        if dry_run:
            _log.info("[dry-run] %s", " ".join(argv))
            return CommandResult(argv, 0, "", "", 0.0)
        self.ensure(spec)
        _log.debug("exec: %s", " ".join(argv))
        started = time.monotonic()
        proc = subprocess.run(argv, capture_output=True, text=True, cwd=str(cwd) if cwd else None)
        return self._finalize(spec, argv, proc, started, expect_outputs, record, check)


class ContainerRunner(_BaseRunner):
    """Runs commands inside an Apptainer/Docker container for reproducibility."""

    def __init__(self, config: ContainerConfig) -> None:
        self.config = config

    def _wrap(self, spec: ToolSpec, argv: list[str]) -> list[str]:
        image = spec.container_image or ""
        if not image:
            raise MissingToolError(
                f"container engine '{self.config.engine}' selected but tool "
                f"'{spec.name}' has no container_image"
            )
        binds: list[str] = []
        for b in self.config.bind_paths:
            binds += ["--bind", b] if self.config.engine == "apptainer" else ["-v", b]
        if self.config.engine == "apptainer":
            return ["apptainer", "exec", *binds, image, *argv]
        return ["docker", "run", "--rm", *binds, image, *argv]

    def run(
        self,
        spec: ToolSpec,
        args: list[str],
        *,
        resources: ResourceSpec | None = None,
        expect_outputs: list[Path] | None = None,
        record: ProvenanceRecord | None = None,
        cwd: Path | None = None,
        check: bool = True,
        dry_run: bool = False,
    ) -> CommandResult:
        inner = [spec.executable, *args]
        argv = self._wrap(spec, inner)
        if dry_run:
            _log.info("[dry-run] %s", " ".join(argv))
            return CommandResult(argv, 0, "", "", 0.0)
        started = time.monotonic()
        proc = subprocess.run(argv, capture_output=True, text=True, cwd=str(cwd) if cwd else None)
        return self._finalize(spec, argv, proc, started, expect_outputs, record, check)


def make_runner(container: ContainerConfig) -> ToolRunner:
    """Factory: pick a runner based on the container configuration."""
    if container.engine == "none":
        return LocalRunner()
    return ContainerRunner(container)


def _version_tuple(v: str) -> tuple[int, ...]:
    return tuple(int(x) for x in re.findall(r"\d+", v))
