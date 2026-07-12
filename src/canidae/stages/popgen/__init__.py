"""Population-genomics module.

Analyses over a materialized genotype matrix: PCA, F_ST, diversity (Phase 1), plus
distance, clustering, admixture, and runs of homozygosity (Phase 2).
"""

from __future__ import annotations

from canidae.stages.popgen import (
    admixture,
    cluster,
    distance,
    diversity,
    fst,
    load,
    pca,
    roh,
)

__all__ = ["admixture", "cluster", "distance", "diversity", "fst", "load", "pca", "roh"]
