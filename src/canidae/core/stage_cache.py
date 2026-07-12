"""Metadata-only safe-resume cache for pipeline stages."""

from __future__ import annotations

import hashlib
import inspect
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

from canidae.core.atomic import atomic_write_text
from canidae.core.hashing import hash_file

if TYPE_CHECKING:
    from canidae.core.model import Artifact
    from canidae.core.stage import RunContext, Stage


class StageCache:
    """Persist only compact stage fingerprints and output handles.

    The cache never copies a VCF, BAM, BCF, or genotype array.  It records fingerprints of
    declared inputs, relevant configured paths, the complete resolved configuration digest,
    and the source file defining the stage class.  A cache hit is accepted only when the
    current registered output files still match the recorded fingerprints.
    """

    schema_version = 1

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def fingerprint(self, stage: Stage, ctx: RunContext, inputs: list[Artifact]) -> dict[str, Any]:
        code = _code_fingerprint(stage)
        payload: dict[str, Any] = {
            "schema_version": self.schema_version,
            "stage": stage.name,
            "config_digest": ctx.config.digest(),
            "stage_config": stage.config.model_dump(mode="json"),
            "code": code,
            "inputs": [_artifact_fingerprint(artifact) for artifact in inputs],
            "configured_paths": _configured_path_fingerprints(stage, ctx),
        }
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
        payload["fingerprint"] = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
        return payload

    def is_current(self, stage: Stage, ctx: RunContext, fingerprint: dict[str, Any]) -> bool:
        if not stage.produced_outputs():
            return False
        cached = self._load(stage.name)
        if cached is None or cached.get("fingerprint") != fingerprint.get("fingerprint"):
            return False

        cached_outputs = {
            (entry.get("kind"), entry.get("role")): entry
            for entry in cached.get("outputs", [])
        }
        for spec in stage.produced_outputs():
            if not ctx.datastore.has(spec.kind, spec.role):
                return False
            artifact = ctx.datastore.get(spec.kind, spec.role)
            recorded = cached_outputs.get((spec.kind.value, spec.role))
            if recorded is None:
                return False
            try:
                current = _artifact_fingerprint(artifact)
            except OSError:
                return False
            if current.get("fingerprint") != recorded.get("fingerprint"):
                return False
        return True

    def store(
        self,
        stage: Stage,
        fingerprint: dict[str, Any],
        outputs: list[Artifact],
    ) -> Path:
        record = {
            **fingerprint,
            "stored_at": datetime.now(UTC).isoformat(timespec="seconds"),
            "outputs": [_artifact_fingerprint(artifact) for artifact in outputs],
        }
        path = self._path_for(stage.name)
        atomic_write_text(path, json.dumps(record, indent=2, sort_keys=True, default=str))
        return path

    def invalidate(self, stage_name: str) -> None:
        self._path_for(stage_name).unlink(missing_ok=True)

    def path_for(self, stage_name: str) -> Path:
        """Return the compact metadata record path for user-facing recovery guidance."""
        return self._path_for(stage_name)

    def _path_for(self, stage_name: str) -> Path:
        safe = "".join(c if c.isalnum() or c in {"-", "_", "."} else "_" for c in stage_name)
        return self.root / f"{safe}.json"

    def _load(self, stage_name: str) -> dict[str, Any] | None:
        path = self._path_for(stage_name)
        if not path.exists():
            return None
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        return data if isinstance(data, dict) else None


def _artifact_fingerprint(artifact: Artifact) -> dict[str, Any]:
    return {
        "kind": artifact.kind.value,
        "role": artifact.role,
        "path": str(artifact.path),
        "fingerprint": hash_file(artifact.path),
    }


def _code_fingerprint(stage: Stage) -> dict[str, str | None]:
    source = inspect.getsourcefile(type(stage))
    if source is None:
        return {"module": type(stage).__module__, "path": None, "fingerprint": None}
    path = Path(source)
    return {
        "module": type(stage).__module__,
        "path": str(path),
        "fingerprint": hash_file(path, mode="full") if path.exists() else None,
        # A stage commonly delegates numerical work to sibling modules (for example,
        # DStats -> fstats).  Fingerprint the compact source tree as well so a helper-code
        # change cannot incorrectly reuse an old stage result.  This is metadata only.
        "package_fingerprint": _package_fingerprint(path),
    }


def _package_fingerprint(source: Path) -> str | None:
    package = next((parent for parent in (source.parent, *source.parents)
                    if parent.name == "canidae"), None)
    if package is None or not package.exists():
        return None
    digest = hashlib.sha256()
    for path in sorted(package.rglob("*.py")):
        digest.update(str(path.relative_to(package)).encode("utf-8"))
        digest.update(hash_file(path, mode="full").encode("ascii"))
    return digest.hexdigest()


def _configured_path_fingerprints(stage: Stage, ctx: RunContext) -> list[dict[str, str | None]]:
    paths = sorted({
        _resolve_configured_path(value, ctx.config.paths.root)
        for value in _walk_paths(stage.config.model_dump(mode="python"))
    }, key=lambda path: str(path))
    entries: list[dict[str, str | None]] = []
    for path in paths:
        entries.append({
            "path": str(path),
            "fingerprint": hash_file(path) if path.exists() else None,
        })
    return entries


def _walk_paths(value: Any) -> list[Path]:
    if isinstance(value, Path):
        return [value]
    if isinstance(value, dict):
        paths: list[Path] = []
        for item in value.values():
            paths.extend(_walk_paths(item))
        return paths
    if isinstance(value, (list, tuple, set)):
        paths = []
        for item in value:
            paths.extend(_walk_paths(item))
        return paths
    return []


def _resolve_configured_path(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else root / path
