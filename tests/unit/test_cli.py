from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from canidae.cli import app
from canidae.version import __version__

runner = CliRunner()


def test_version_command() -> None:
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert __version__ in result.stdout


def test_config_show() -> None:
    result = runner.invoke(app, ["config", "--set", "seed=99"])
    assert result.exit_code == 0
    assert "seed: 99" in result.stdout


def test_stages_empty_ok() -> None:
    result = runner.invoke(app, ["stages"])
    assert result.exit_code == 0


def test_run_dry_run_with_empty_pipeline() -> None:
    result = runner.invoke(app, ["run", "--dry-run"])
    assert result.exit_code == 0
    assert "nothing to run" in result.stdout


def test_run_dry_run_pipeline(tmp_path: Path) -> None:
    # A project config that references the built-in dummy 'noop' stage would go here once
    # real stages exist; for now assert the empty-pipeline guard path.
    cfg = tmp_path / "p.yaml"
    cfg.write_text("pipeline: []\n", encoding="utf-8")
    result = runner.invoke(app, ["run", "-c", str(cfg), "--dry-run"])
    assert result.exit_code == 0
