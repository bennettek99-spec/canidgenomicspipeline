"""Provenance capture for CANIS.

Reproducibility is a first-class output. Every stage emits a :class:`ProvenanceRecord`
capturing exactly how its artifacts were produced: the tool versions and command lines,
input/output content hashes, the config digest, the git commit, RNG seed, host, and
timing. A :class:`ProvenanceWriter` accumulates records for a run and writes both a
machine-readable ``manifest.json`` and a human-readable ``manifest.md``.

A run is intended to be reconstructible from its manifest alone.
"""

from __future__ import annotations

import getpass
import json
import platform
import socket
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from canidae.core.atomic import atomic_write_text
from canidae.core.model import Artifact
from canidae.version import __version__


@dataclass(slots=True)
class ToolInvocation:
    """A single external command executed by a stage."""

    tool: str
    version: str | None
    argv: list[str]
    returncode: int | None = None
    duration_s: float | None = None


@dataclass(slots=True)
class ProvenanceRecord:
    """Everything needed to explain (and re-run) one stage."""

    id: str
    stage: str
    started_at: str
    status: str = "running"  # running | succeeded | failed | skipped
    finished_at: str | None = None
    duration_s: float | None = None
    config_digest: str = ""
    seed: int | None = None
    git_commit: str | None = None
    inputs: list[dict[str, Any]] = field(default_factory=list)
    outputs: list[dict[str, Any]] = field(default_factory=list)
    tools: list[ToolInvocation] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    environment: dict[str, str] = field(default_factory=dict)
    skip_reason: str | None = None
    error: str | None = None
    recovery: str | None = None

    def add_tool(self, invocation: ToolInvocation) -> None:
        self.tools.append(invocation)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _now_iso() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _git_commit(cwd: Path | None = None) -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(cwd) if cwd else None,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return out.stdout.strip() or None if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def _environment_snapshot() -> dict[str, str]:
    return {
        "canidae_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "host": socket.gethostname(),
        "user": _safe_user(),
    }


def _safe_user() -> str:
    try:
        return getpass.getuser()
    except Exception:  # pragma: no cover
        return "unknown"


class ProvenanceWriter:
    """Accumulates provenance records for a run and persists the manifest."""

    def __init__(
        self,
        run_dir: Path,
        *,
        config_digest: str = "",
        seed: int | None = None,
        repo_dir: Path | None = None,
    ) -> None:
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.config_digest = config_digest
        self.seed = seed
        self.git_commit = _git_commit(repo_dir)
        self.records: list[ProvenanceRecord] = []
        self.run_metadata: dict[str, Any] = {}
        self.run_events: list[dict[str, Any]] = []
        self._counter = 0
        self._start_times: dict[str, float] = {}

    def start(self, stage: str) -> ProvenanceRecord:
        """Open a provenance record for a stage. Call :meth:`finish` when done."""
        self._counter += 1
        record = ProvenanceRecord(
            id=f"{stage}-{self._counter:03d}",
            stage=stage,
            started_at=_now_iso(),
            config_digest=self.config_digest,
            seed=self.seed,
            git_commit=self.git_commit,
            environment=_environment_snapshot(),
        )
        self._start_times[record.id] = time.monotonic()
        # Persist the explicit running state immediately. If a process is interrupted,
        # the manifest tells the user which stage needs recovery instead of omitting it.
        self.records.append(record)
        self.flush()
        return record

    def record_inputs(self, record: ProvenanceRecord, artifacts: list[Artifact]) -> None:
        record.inputs.extend(_artifact_stub(a) for a in artifacts)

    def finish(
        self,
        record: ProvenanceRecord,
        *,
        outputs: list[Artifact] | None = None,
        metrics: dict[str, Any] | None = None,
        status: str = "succeeded",
        skip_reason: str | None = None,
        error: str | None = None,
        recovery: str | None = None,
    ) -> ProvenanceRecord:
        if status not in {"succeeded", "failed", "skipped"}:
            raise ValueError(f"invalid stage provenance status: {status}")
        t0 = self._start_times.pop(record.id, None)
        record.status = status
        record.finished_at = _now_iso()
        record.duration_s = round(time.monotonic() - t0, 3) if t0 else None
        record.skip_reason = skip_reason
        record.error = error
        record.recovery = recovery
        if outputs:
            record.outputs.extend(_artifact_stub(a) for a in outputs)
        if metrics:
            record.metrics.update(metrics)
        if record not in self.records:
            self.records.append(record)
        self.flush()  # persist incrementally so a crash still leaves a partial manifest
        return record

    def add_run_metadata(self, **metadata: Any) -> None:
        """Persist small run-wide context such as disk preflight and thread limits."""
        self.run_metadata.update(metadata)
        self.flush()

    def record_run_event(
        self,
        event: str,
        *,
        status: str,
        error: str | None = None,
        recovery: str | None = None,
        details: dict[str, Any] | None = None,
    ) -> None:
        self.run_events.append({
            "event": event,
            "status": status,
            "at": _now_iso(),
            "error": error,
            "recovery": recovery,
            "details": details or {},
        })
        self.flush()

    def flush(self) -> Path:
        manifest = {
            "run_dir": str(self.run_dir),
            "generated_at": _now_iso(),
            "config_digest": self.config_digest,
            "git_commit": self.git_commit,
            "seed": self.seed,
            "run_metadata": self.run_metadata,
            "run_events": self.run_events,
            "records": [r.to_dict() for r in self.records],
        }
        path = self.run_dir / "manifest.json"
        atomic_write_text(path, json.dumps(manifest, indent=2, default=str))
        atomic_write_text(self.run_dir / "manifest.md", self._render_markdown())
        return path

    def _render_markdown(self) -> str:
        lines = [
            "# CANIS run manifest",
            "",
            f"- generated: {_now_iso()}",
            f"- config digest: `{self.config_digest[:16]}`",
            f"- git commit: `{self.git_commit or 'n/a'}`",
            f"- seed: {self.seed}",
            "",
            "## Stages",
            "",
        ]
        for r in self.records:
            lines.append(f"### {r.id} — {r.stage}")
            lines.append(f"- status: {r.status}")
            lines.append(f"- duration: {r.duration_s}s")
            if r.skip_reason:
                lines.append(f"- skip reason: {r.skip_reason}")
            if r.error:
                lines.append(f"- error: {r.error}")
            if r.recovery:
                lines.append(f"- recovery: {r.recovery}")
            if r.tools:
                lines.append("- tools:")
                lines += [
                    f"    - `{t.tool}` ({t.version or 'n/a'}) rc={t.returncode}"
                    for t in r.tools
                ]
            if r.outputs:
                lines.append("- outputs:")
                lines += [f"    - {o['kind']}:{o['role']} → `{o['path']}`" for o in r.outputs]
            lines.append("")
        if self.run_events:
            lines.extend(["## Run events", ""])
            for event in self.run_events:
                lines.append(f"- {event['event']}: {event['status']}")
                if event.get("error"):
                    lines.append(f"  - error: {event['error']}")
                if event.get("recovery"):
                    lines.append(f"  - recovery: {event['recovery']}")
            lines.append("")
        return "\n".join(lines)


def _artifact_stub(a: Artifact) -> dict[str, Any]:
    return {
        "kind": a.kind.value,
        "role": a.role,
        "path": str(a.path),
        "fmt": a.fmt.value,
        "checksum": a.checksum,
        "produced_by": a.produced_by,
    }
