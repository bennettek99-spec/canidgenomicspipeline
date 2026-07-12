"""Logging setup for CANIS.

Provides a single :func:`configure_logging` entry point that installs a rich console
handler and, when a run directory is given, a per-run file handler. A small
:func:`get_logger` wrapper keeps logger names namespaced under ``canidae``.

We keep this dependency-light: rich is used when available and configured, otherwise we
fall back to the standard :class:`logging.StreamHandler`.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from canidae.core.config import LoggingConfig

_LOGGER_ROOT = "canidae"
_CONFIGURED = False


class _JsonFormatter(logging.Formatter):
    """Minimal structured formatter: one JSON object per line."""

    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, Any] = {
            "ts": self.formatTime(record, "%Y-%m-%dT%H:%M:%S"),
            "level": record.levelname,
            "logger": record.name,
            "msg": record.getMessage(),
        }
        if record.exc_info:
            payload["exc"] = self.formatException(record.exc_info)
        for key, value in getattr(record, "extra_fields", {}).items():
            payload[key] = value
        return json.dumps(payload, default=str)


def configure_logging(
    config: LoggingConfig | None = None,
    *,
    run_dir: Path | None = None,
) -> logging.Logger:
    """Configure the ``canidae`` logger hierarchy. Idempotent per process."""
    global _CONFIGURED
    config = config or LoggingConfig()
    root = logging.getLogger(_LOGGER_ROOT)
    root.setLevel(config.level)
    root.handlers.clear()
    root.propagate = False

    console_handler = _make_console_handler(config)
    root.addHandler(console_handler)

    if run_dir is not None:
        Path(run_dir).mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(Path(run_dir) / "run.log", encoding="utf-8")
        file_handler.setFormatter(_JsonFormatter())
        root.addHandler(file_handler)

    _CONFIGURED = True
    return root


def _make_console_handler(config: LoggingConfig) -> logging.Handler:
    if config.rich_console and not config.json_format:
        try:
            from rich.logging import RichHandler

            handler: logging.Handler = RichHandler(
                rich_tracebacks=True, show_path=False, markup=False
            )
            handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%X]"))
            return handler
        except Exception:  # pragma: no cover - rich always present in practice
            pass
    handler = logging.StreamHandler()
    handler.setFormatter(
        _JsonFormatter()
        if config.json_format
        else logging.Formatter("%(asctime)s %(levelname)-7s %(name)s | %(message)s")
    )
    return handler


def get_logger(name: str | None = None) -> logging.Logger:
    """Return a namespaced logger. ``get_logger("stage.qc")`` -> ``canidae.stage.qc``."""
    if not _CONFIGURED:
        configure_logging()
    return logging.getLogger(_LOGGER_ROOT if not name else f"{_LOGGER_ROOT}.{name}")
