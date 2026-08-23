from __future__ import annotations

from pathlib import Path

import pytest

from canidae.core.config import GlobalConfig
from canidae.core.errors import ConfigError


def test_defaults_load() -> None:
    cfg = GlobalConfig.load()
    assert cfg.project_name == "canis"
    assert cfg.executor.backend == "native"
    assert cfg.executor.max_workers == 4


def test_layered_override(tmp_path: Path) -> None:
    project = tmp_path / "project.yaml"
    project.write_text("project_name: wolves\nexecutor:\n  max_workers: 16\n", encoding="utf-8")
    cfg = GlobalConfig.load(project)
    assert cfg.project_name == "wolves"
    assert cfg.executor.max_workers == 16
    # unspecified values fall back to packaged defaults
    assert cfg.executor.fail_fast is True


def test_dotted_cli_override_coercion() -> None:
    cfg = GlobalConfig.load(overrides={"executor.max_workers": "8", "seed": "42"})
    assert cfg.executor.max_workers == 8
    assert cfg.seed == 42


def test_digest_is_stable_and_sensitive() -> None:
    a = GlobalConfig.load(overrides={"seed": "1"})
    b = GlobalConfig.load(overrides={"seed": "1"})
    c = GlobalConfig.load(overrides={"seed": "2"})
    assert a.digest() == b.digest()
    assert a.digest() != c.digest()


def test_paths_resolved_absolute() -> None:
    cfg = GlobalConfig.load()
    assert cfg.paths.data_root.is_absolute()


def test_extra_key_rejected(tmp_path: Path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("executor:\n  nonsense: 1\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        GlobalConfig.load(bad)


def test_missing_file_raises() -> None:
    with pytest.raises(ConfigError):
        GlobalConfig.load("does-not-exist.yaml")
