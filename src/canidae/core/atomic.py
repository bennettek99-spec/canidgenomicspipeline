"""Small same-filesystem atomic-write primitives used by core persistence layers.

The pipeline never needs to copy a large genomic artifact merely to make metadata safe.
These helpers keep metadata writes and stage-directory promotion on the same filesystem so
readers see either the old complete value or the new complete value.
"""

from __future__ import annotations

import os
import shutil
import uuid
from pathlib import Path


def atomic_write_text(path: Path, text: str, *, encoding: str = "utf-8") -> None:
    """Durably replace a small text file without exposing a partial write."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("w", encoding=encoding, newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink(missing_ok=True)


def atomic_replace_directory(source: Path, target: Path) -> None:
    """Promote ``source`` to ``target`` with rollback of an existing target.

    ``source`` and ``target`` must live on the same filesystem.  If a previous result is
    present, it is first renamed to a private backup, then removed only after the new
    directory has been promoted.  A failure between those operations restores the backup.
    """
    source = Path(source)
    target = Path(target)
    if not source.is_dir():
        raise ValueError(f"staged output directory does not exist: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)

    backup: Path | None = None
    if target.exists():
        backup = target.with_name(f".{target.name}.previous-{uuid.uuid4().hex}")

    try:
        if backup is not None:
            os.replace(target, backup)
        os.replace(source, target)
    except Exception:
        if backup is not None and backup.exists() and not target.exists():
            os.replace(backup, target)
        raise
    else:
        if backup is not None and backup.exists():
            shutil.rmtree(backup)
