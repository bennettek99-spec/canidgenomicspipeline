"""Content hashing helpers used by the datastore and provenance layers.

Hashing large genomic files fully can be expensive. We therefore support two modes:

* ``full``  — stream the whole file through SHA-256 (default for small/medium files).
* ``quick`` — hash ``(size, mtime_ns, head_bytes, tail_bytes)`` for very large files, a
  cheap fingerprint sufficient for cache invalidation (not for integrity guarantees).

Directories (e.g. a ``.zarr`` store) are hashed as a manifest of relative paths + per-file
quick fingerprints.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Literal

HashMode = Literal["full", "quick"]

_CHUNK = 1 << 20  # 1 MiB
_QUICK_EDGE = 1 << 16  # 64 KiB from head and tail
_QUICK_THRESHOLD = 1 << 28  # 256 MiB: above this, default to quick fingerprints


def hash_file(path: Path, mode: HashMode | None = None) -> str:
    """Return a hex digest for ``path``. If ``mode`` is None, choose based on size."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.is_dir():
        return hash_dir(path)
    size = path.stat().st_size
    chosen = mode or ("quick" if size >= _QUICK_THRESHOLD else "full")
    return _quick_fingerprint(path, size) if chosen == "quick" else _full_sha256(path)


def _full_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return f"sha256:{h.hexdigest()}"


def _quick_fingerprint(path: Path, size: int) -> str:
    st = path.stat()
    h = hashlib.sha256()
    h.update(str(size).encode())
    h.update(str(st.st_mtime_ns).encode())
    with path.open("rb") as fh:
        h.update(fh.read(_QUICK_EDGE))
        if size > _QUICK_EDGE:
            fh.seek(max(0, size - _QUICK_EDGE))
            h.update(fh.read(_QUICK_EDGE))
    return f"quick:{h.hexdigest()}"


def hash_dir(path: Path) -> str:
    """Hash a directory as a stable manifest of relative paths + quick fingerprints."""
    h = hashlib.sha256()
    for child in sorted(p for p in path.rglob("*") if p.is_file()):
        rel = child.relative_to(path).as_posix()
        h.update(rel.encode())
        h.update(_quick_fingerprint(child, child.stat().st_size).encode())
    return f"dir:{h.hexdigest()}"


def hash_bytes(data: bytes) -> str:
    return f"sha256:{hashlib.sha256(data).hexdigest()}"
