"""CANIS — Canid Evolutionary Genomics Pipeline.

A configuration-driven, reproducible platform for comparative evolutionary genomics of
canids. The top-level package intentionally re-exports only lightweight symbols so that
``import canidae`` stays cheap; heavy stage machinery is imported from its subpackages on
demand.
"""

from __future__ import annotations

from canidae.version import __version__

__all__ = ["__version__"]
