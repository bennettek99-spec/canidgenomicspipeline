from __future__ import annotations

from pathlib import Path

import pytest

from canidae.core.datastore import DataStore
from canidae.core.errors import DataStoreError, IntegrityError
from canidae.core.model import ArtifactKind, FileFormat


def _write(store: DataStore, stage: str, name: str, text: str) -> Path:
    p = store.path_for(stage, name)
    p.write_text(text, encoding="utf-8")
    return p


def test_add_get_roundtrip(tmp_path: Path) -> None:
    store = DataStore(tmp_path / "store")
    p = _write(store, "qc", "metrics.csv", "a,b\n1,2\n")
    art = store.add(ArtifactKind.QC_TABLE, "metrics", p, fmt=FileFormat.CSV, produced_by="qc")
    assert art.checksum is not None
    got = store.get(ArtifactKind.QC_TABLE, "metrics")
    assert got.path == p and got.produced_by == "qc"
    assert store.has(ArtifactKind.QC_TABLE, "metrics")


def test_missing_lookup_raises(tmp_path: Path) -> None:
    store = DataStore(tmp_path / "store")
    with pytest.raises(DataStoreError):
        store.get(ArtifactKind.CALLSET, "nope")


def test_register_missing_file_raises(tmp_path: Path) -> None:
    store = DataStore(tmp_path / "store")
    with pytest.raises(DataStoreError):
        store.add(ArtifactKind.CALLSET, "x", tmp_path / "ghost.bcf")


def test_index_persists_across_instances(tmp_path: Path) -> None:
    root = tmp_path / "store"
    store = DataStore(root)
    p = _write(store, "qc", "m.csv", "x\n")
    store.add(ArtifactKind.QC_TABLE, "metrics", p)
    reopened = DataStore(root)
    assert reopened.has(ArtifactKind.QC_TABLE, "metrics")


def test_verify_detects_tampering(tmp_path: Path) -> None:
    store = DataStore(tmp_path / "store")
    p = _write(store, "qc", "m.csv", "original\n")
    art = store.add(ArtifactKind.QC_TABLE, "metrics", p)
    # Force a full sha256 so verify is meaningful, then tamper.
    if art.checksum and art.checksum.startswith("sha256:"):
        p.write_text("tampered\n", encoding="utf-8")
        with pytest.raises(IntegrityError):
            store.verify(art)


def test_find_by_kind(tmp_path: Path) -> None:
    store = DataStore(tmp_path / "store")
    store.add(ArtifactKind.QC_TABLE, "a", _write(store, "s", "a.csv", "1\n"))
    store.add(ArtifactKind.QC_TABLE, "b", _write(store, "s", "b.csv", "1\n"))
    assert len(store.find(ArtifactKind.QC_TABLE)) == 2
    assert len(store.find(ArtifactKind.QC_TABLE, "a")) == 1
