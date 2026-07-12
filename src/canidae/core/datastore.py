"""The DataStore: the shared, typed artifact registry modules exchange through.

Because modules never call each other, the DataStore is the single mechanism by which a
downstream stage discovers upstream outputs. It:

* registers :class:`~canidae.core.model.Artifact` handles keyed by ``(kind, role)``,
* persists an on-disk ``artifacts.json`` index so runs are resumable,
* verifies integrity via content hashes, and
* offers a deterministic path layout for stages to write their outputs into.

It stores *handles and metadata*, never the genomic data itself (which lives on disk in
CRAM/BCF/Zarr under the run/data directories).
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from canidae.core.atomic import atomic_write_text
from canidae.core.errors import DataStoreError, IntegrityError
from canidae.core.hashing import hash_file
from canidae.core.logging import get_logger
from canidae.core.model import Artifact, ArtifactKind, FileFormat

_log = get_logger("datastore")


class DataStore:
    """A content-addressable workspace for a run's artifacts."""

    def __init__(self, root: Path) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self._index_path = self.root / "artifacts.json"
        self._artifacts: dict[tuple[ArtifactKind, str], Artifact] = {}
        self._recover_interrupted_promotions()
        if self._index_path.exists():
            self._load_index()

    # -- path layout -------------------------------------------------------------------

    def stage_dir(self, stage_name: str) -> Path:
        """A dedicated output directory for a stage; created on demand."""
        path = self.root / stage_name
        path.mkdir(parents=True, exist_ok=True)
        return path

    def path_for(self, stage_name: str, filename: str) -> Path:
        return self.stage_dir(stage_name) / filename

    # -- registration ------------------------------------------------------------------

    def register(self, artifact: Artifact, *, compute_checksum: bool = True) -> Artifact:
        """Register (or replace) an artifact, optionally computing its checksum.

        Registering the same ``(kind, role)`` again replaces the handle — the most recent
        producer wins, which is what re-running a stage should do.
        """
        stored = self.register_many([artifact], compute_checksum=compute_checksum)[0]
        _log.debug("registered artifact %s:%s -> %s", stored.kind.value, stored.role,
                   stored.path)
        return stored

    def register_many(
        self,
        artifacts: list[Artifact],
        *,
        compute_checksum: bool = True,
    ) -> list[Artifact]:
        """Validate and register a stage's outputs in one atomic index update.

        The underlying files are already promoted by the executor.  This method keeps a
        multi-output stage from exposing an index with only some of its output handles.
        """
        stored: list[Artifact] = []
        for artifact in artifacts:
            if not artifact.path.exists():
                raise DataStoreError(
                    f"cannot register missing artifact {artifact.kind.value}:{artifact.role} "
                    f"at {artifact.path}"
                )
            current = artifact
            if compute_checksum and artifact.checksum is None:
                current = Artifact(
                    kind=artifact.kind,
                    role=artifact.role,
                    path=artifact.path,
                    fmt=artifact.fmt,
                    checksum=hash_file(artifact.path),
                    produced_by=artifact.produced_by,
                    provenance_id=artifact.provenance_id,
                    schema_version=artifact.schema_version,
                    metadata=artifact.metadata,
                )
            stored.append(current)

        previous = self._artifacts
        updated = dict(previous)
        for artifact in stored:
            updated[artifact.key] = artifact
        self._artifacts = updated
        try:
            self._save_index()
        except Exception:
            self._artifacts = previous
            raise
        return stored

    def add(
        self,
        kind: ArtifactKind,
        role: str,
        path: Path,
        *,
        fmt: FileFormat = FileFormat.OTHER,
        produced_by: str | None = None,
        provenance_id: str | None = None,
        metadata: dict | None = None,
    ) -> Artifact:
        """Convenience constructor + register."""
        artifact = Artifact(
            kind=kind, role=role, path=Path(path), fmt=fmt, produced_by=produced_by,
            provenance_id=provenance_id, metadata=metadata or {},
        )
        return self.register(artifact)

    # -- lookup ------------------------------------------------------------------------

    def get(self, kind: ArtifactKind, role: str) -> Artifact:
        try:
            return self._artifacts[(kind, role)]
        except KeyError:
            raise DataStoreError(
                f"no artifact registered for {kind.value}:{role}"
            ) from None

    def find(self, kind: ArtifactKind, role: str | None = None) -> list[Artifact]:
        """Return all artifacts of a kind, optionally filtered by role."""
        return [
            a
            for (k, r), a in self._artifacts.items()
            if k == kind and (role is None or r == role)
        ]

    def has(self, kind: ArtifactKind, role: str) -> bool:
        return (kind, role) in self._artifacts

    def all(self) -> list[Artifact]:
        return list(self._artifacts.values())

    # -- integrity ---------------------------------------------------------------------

    def verify(self, artifact: Artifact) -> None:
        """Raise :class:`IntegrityError` if the artifact is missing or its hash changed."""
        if not artifact.path.exists():
            raise IntegrityError(f"artifact vanished: {artifact.path}")
        if artifact.checksum is None:
            return
        # Only full sha256 checksums are re-verifiable; quick fingerprints are advisory.
        if artifact.checksum.startswith("sha256:"):
            current = hash_file(artifact.path, mode="full")
            if current != artifact.checksum:
                raise IntegrityError(
                    f"checksum mismatch for {artifact.path}: "
                    f"expected {artifact.checksum}, got {current}"
                )

    # -- persistence -------------------------------------------------------------------

    def _save_index(self) -> None:
        payload = [
            {
                "kind": a.kind.value,
                "role": a.role,
                "path": str(a.path),
                "fmt": a.fmt.value,
                "checksum": a.checksum,
                "produced_by": a.produced_by,
                "provenance_id": a.provenance_id,
                "schema_version": a.schema_version,
                "metadata": a.metadata,
            }
            for a in self._artifacts.values()
        ]
        atomic_write_text(self._index_path, json.dumps(payload, indent=2, default=str))

    def _recover_interrupted_promotions(self) -> None:
        """Restore a previous stage directory if a process died mid-promotion.

        ``atomic_replace_directory`` normally rolls back immediately.  This inexpensive
        startup check covers a hard process kill between the two renames.
        """
        marker = ".previous-"
        for backup in self.root.glob(".*.previous-*"):
            name = backup.name
            if not name.startswith(".") or marker not in name:
                continue
            stage_name = name[1:].split(marker, 1)[0]
            if not stage_name or Path(stage_name).name != stage_name:
                continue
            target = self.root / stage_name
            try:
                if target.exists():
                    if backup.is_dir():
                        shutil.rmtree(backup)
                    else:
                        backup.unlink()
                else:
                    os.replace(backup, target)
                    _log.warning("recovered interrupted stage-output promotion: %s", target)
            except OSError as exc:
                _log.warning("could not recover staged output backup %s: %s", backup, exc)

    def _load_index(self) -> None:
        try:
            payload = json.loads(self._index_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise DataStoreError(f"corrupt artifact index at {self._index_path}: {exc}") \
                from exc
        for row in payload:
            artifact = Artifact(
                kind=ArtifactKind(row["kind"]),
                role=row["role"],
                path=Path(row["path"]),
                fmt=FileFormat(row.get("fmt", "other")),
                checksum=row.get("checksum"),
                produced_by=row.get("produced_by"),
                provenance_id=row.get("provenance_id"),
                schema_version=row.get("schema_version", 1),
                metadata=row.get("metadata", {}),
            )
            self._artifacts[artifact.key] = artifact
