"""Hybrid-canid diagnostic stages (bridge-panel mixture / breed assignment).

These stages answer coydog / eastern-coyote style questions on small cross-platform
bridge panels. Results are exploratory diagnostics, not whole-genome ancestry.
"""

from __future__ import annotations

from canidae.stages.hybrid import breed_assign, multiway_admixture, reference_mixture

__all__ = ["breed_assign", "multiway_admixture", "reference_mixture"]
