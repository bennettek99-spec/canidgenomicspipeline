from __future__ import annotations

import sys

import pytest

from canidae.core.errors import ExternalToolError, MissingToolError
from canidae.core.provenance import ProvenanceWriter
from canidae.core.runtime import (
    LocalRunner,
    ResourceSpec,
    ToolSpec,
    _version_tuple,
)


def test_version_tuple_ordering() -> None:
    assert _version_tuple("1.9.2") < _version_tuple("1.10.0")
    assert _version_tuple("2.0") > _version_tuple("1.99.99")


def test_missing_tool_raises() -> None:
    runner = LocalRunner()
    spec = ToolSpec(name="definitely-not-a-real-binary-xyz")
    assert runner.which(spec) is None
    with pytest.raises(MissingToolError):
        runner.ensure(spec)


def test_dry_run_skips_execution() -> None:
    runner = LocalRunner()
    spec = ToolSpec(name="definitely-not-a-real-binary-xyz")
    result = runner.run(spec, ["--help"], dry_run=True)
    assert result.ok and result.returncode == 0


def test_real_command_records_provenance(tmp_path) -> None:
    # Use the current Python interpreter as a guaranteed-present executable.
    runner = LocalRunner()
    spec = ToolSpec(name="python", binary=sys.executable, version_args=("--version",))
    prov = ProvenanceWriter(tmp_path, config_digest="d", seed=1)
    record = prov.start("demo")
    result = runner.run(
        spec, ["-c", "print('hello')"], record=record, resources=ResourceSpec(cpus=1)
    )
    assert result.ok and "hello" in result.stdout
    prov.finish(record)
    assert record.tools and record.tools[0].tool == "python"


def test_nonzero_exit_raises(tmp_path) -> None:
    runner = LocalRunner()
    spec = ToolSpec(name="python", binary=sys.executable)
    with pytest.raises(ExternalToolError):
        runner.run(spec, ["-c", "import sys; sys.exit(3)"])
