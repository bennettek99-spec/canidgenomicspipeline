"""Per-stage transactional datastore view.

Stages historically write directly to :class:`DataStore` paths.  ``StagedDataStore`` keeps
that public interface but redirects a running stage's own output directory into a temporary
same-filesystem workspace.  The executor validates returned artifacts, promotes the complete
directory, and only then updates the shared artifact index.
"""

from __future__ import annotations

import shutil
import uuid
from dataclasses import replace
from pathlib import Path
from typing import Any

from canidae.core.atomic import atomic_replace_directory
from canidae.core.datastore import DataStore
from canidae.core.errors import DataStoreError
from canidae.core.hashing import hash_file
from canidae.core.model import Artifact, ArtifactKind, FileFormat


class StagedDataStore:
    """A write-isolated view of a datastore for one stage execution.

    Reads delegate to the shared datastore.  Registrations are retained locally until the
    executor has promoted the staged output directory, preventing half-written artifacts from
    becoming visible to downstream stages.
    """

    def __init__(self, base: DataStore, stage_name: str) -> None:
        self.base = base
        self.stage_name = stage_name
        self.root = base.root
        self._workspace = base.root / ".staging" / f"{stage_name}-{uuid.uuid4().hex}"
        self._stage_root = self._workspace / stage_name
        self._local_artifacts: dict[tuple[ArtifactKind, str], Artifact] = {}

    @property
    def workspace(self) -> Path:
        return self._workspace

    @property
    def stage_root(self) -> Path:
        return self._stage_root

    # -- redirected writes -----------------------------------------------------------

    def stage_dir(self, stage_name: str) -> Path:
        if stage_name != self.stage_name:
            return self.base.stage_dir(stage_name)
        self._stage_root.mkdir(parents=True, exist_ok=True)
        return self._stage_root

    def path_for(self, stage_name: str, filename: str) -> Path:
        return self.stage_dir(stage_name) / filename

    def register(self, artifact: Artifact, *, compute_checksum: bool = True) -> Artifact:
        """Record a local handle without mutating the shared index.

        Checksums are deliberately computed only after promotion.  The executor's
        ``DataStore.register_many`` call then validates all final paths in one metadata
        transaction.
        """
        if not artifact.path.exists():
            raise DataStoreError(
                f"cannot register missing artifact {artifact.kind.value}:{artifact.role} "
                f"at {artifact.path}"
            )
        self._local_artifacts[artifact.key] = artifact
        return artifact

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
        return self.register(
            Artifact(
                kind=kind,
                role=role,
                path=Path(path),
                fmt=fmt,
                produced_by=produced_by,
                provenance_id=provenance_id,
                metadata=metadata or {},
            )
        )

    # -- reads -----------------------------------------------------------------------

    def get(self, kind: ArtifactKind, role: str) -> Artifact:
        return self._local_artifacts.get((kind, role)) or self.base.get(kind, role)

    def has(self, kind: ArtifactKind, role: str) -> bool:
        return (kind, role) in self._local_artifacts or self.base.has(kind, role)

    def find(self, kind: ArtifactKind, role: str | None = None) -> list[Artifact]:
        merged = {artifact.key: artifact for artifact in self.base.find(kind, role)}
        merged.update(
            {
                key: artifact
                for key, artifact in self._local_artifacts.items()
                if key[0] == kind and (role is None or key[1] == role)
            }
        )
        return list(merged.values())

    def all(self) -> list[Artifact]:
        merged = {artifact.key: artifact for artifact in self.base.all()}
        merged.update(self._local_artifacts)
        return list(merged.values())

    def verify(self, artifact: Artifact) -> None:
        if artifact.key in self._local_artifacts:
            if not artifact.path.exists():
                raise DataStoreError(f"artifact vanished before promotion: {artifact.path}")
            return
        self.base.verify(artifact)

    # -- transaction finalization ----------------------------------------------------

    def promote(self, artifacts: list[Artifact]) -> list[Artifact]:
        """Validate and atomically promote this stage's local output directory.

        Artifacts outside the staging directory (for example, a referenced public VCF)
        remain in place but are still subsequently validated by the base datastore.
        """
        for artifact in artifacts:
            if not artifact.path.exists():
                raise DataStoreError(
                    f"stage '{self.stage_name}' returned a missing output "
                    f"{artifact.kind.value}:{artifact.role} at {artifact.path}"
                )
            # A read verifies accessibility and produces no duplicate data.  It also catches
            # a directory whose member files disappeared before promotion.
            hash_file(artifact.path)

        staged_exists = self._stage_root.exists()
        final_root = self.base.root / self.stage_name
        if staged_exists:
            atomic_replace_directory(self._stage_root, final_root)

        promoted: list[Artifact] = []
        for artifact in artifacts:
            try:
                relative = artifact.path.relative_to(self._stage_root)
            except ValueError:
                promoted.append(artifact)
            else:
                promoted.append(
                    replace(
                        artifact,
                        path=final_root / relative,
                        metadata=rewrite_staged_paths(
                            artifact.metadata, self._stage_root, final_root
                        ),
                    )
                )
        self.cleanup()
        return promoted

    def cleanup(self) -> None:
        """Remove unpromoted transaction files; shared final outputs are untouched."""
        if self._workspace.exists():
            shutil.rmtree(self._workspace, ignore_errors=True)

    def __getattr__(self, name: str) -> Any:
        """Delegate uncommon read-only datastore APIs without widening this wrapper."""
        return getattr(self.base, name)


def rewrite_staged_paths(value: Any, staged_root: Path, final_root: Path) -> Any:
    """Recursively rewrite temporary transaction paths after atomic promotion."""
    staged = str(staged_root)
    final = str(final_root)
    if isinstance(value, str):
        return value.replace(staged, final)
    if isinstance(value, dict):
        return {
            key: rewrite_staged_paths(item, staged_root, final_root) for key, item in value.items()
        }
    if isinstance(value, list):
        return [rewrite_staged_paths(item, staged_root, final_root) for item in value]
    if isinstance(value, tuple):
        return tuple(rewrite_staged_paths(item, staged_root, final_root) for item in value)
    return value
