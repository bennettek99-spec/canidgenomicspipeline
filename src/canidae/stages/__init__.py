"""Pipeline stages.

Each subpackage here implements one pipeline module as one or more
:class:`~canidae.core.stage.Stage` subclasses, registered against
:data:`canidae.core.registry.STAGES`. Stages are the *only* place external tools are
invoked and the *only* place analysis logic lives; they communicate exclusively through
typed artifacts in the datastore.

Implemented (Phase 1): acquisition · qc · popgen · reporting
Planned: processing · phylogenetics · introgression · comparative · geographic
Future: selection · demography · ancient · sv
"""

from __future__ import annotations

import importlib

from canidae.core.logging import get_logger

_log = get_logger("stages")

# Built-in stage packages, imported on demand. Importing a package runs its stages'
# registration decorators. Kept lazy so ``import canidae`` stays cheap and so a
# foundation-only install (without the analysis extras) still works.
_BUILTIN_PACKAGES = (
    "canidae.stages.acquisition",
    "canidae.stages.processing",
    "canidae.stages.qc",
    "canidae.stages.popgen",
    "canidae.stages.phylogenetics",
    "canidae.stages.introgression",
    "canidae.stages.local_ancestry",
    "canidae.stages.selection",
    "canidae.stages.demography",
    "canidae.stages.ancient",
    "canidae.stages.geographic",
    "canidae.stages.sv",
    "canidae.stages.reporting",
)


def load_builtin_stages() -> int:
    """Import all built-in stage packages so they self-register. Returns the count that
    loaded successfully; packages whose optional dependencies are missing are skipped
    with a warning rather than raising."""
    loaded = 0
    for pkg in _BUILTIN_PACKAGES:
        try:
            importlib.import_module(pkg)
            loaded += 1
        except ImportError as exc:  # optional analysis deps absent
            _log.warning("stage package '%s' unavailable: %s", pkg, exc)
    return loaded
