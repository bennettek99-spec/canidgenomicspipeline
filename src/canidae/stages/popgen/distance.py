"""Pairwise genetic-distance matrix between individuals.

Computes an allele-difference distance: the mean per-site absolute difference in ALT-allele
count, normalized by two to [0, 1] (0 = identical, 1 = maximally different). This feeds clustering,
neighbor-joining trees (Phase 3), and isolation-by-distance (geographic module).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import allele_difference_matrix, load_genotypes


class DistanceConfig(StageConfig):
    metric: str = "allele_difference"


@STAGES.register("distance")
class DistanceStage(Stage):
    name = "distance"
    config_model = DistanceConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "distance")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: DistanceConfig = self.config  # type: ignore[assignment]
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        ac = geno.calls.count_alleles()
        n_alt = geno.calls.to_n_alt()[ac.is_segregating()]  # (n_seg_sites, n_samples)
        if n_alt.shape[0] < 1:
            raise StageInputError("no segregating sites for distance computation")

        dmatrix = allele_difference_matrix(n_alt)  # mean ALT-count difference / 2 -> [0, 1]
        matrix = pd.DataFrame(dmatrix, index=geno.samples, columns=geno.samples)

        out = ctx.datastore.path_for(self.name, "distance_matrix.csv")
        matrix.to_csv(out)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "distance", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "distance", "metric": cfg.metric,
                       "n_sites": int(n_alt.shape[0])},
        )
        iu = np.triu_indices(dmatrix.shape[0], k=1)
        return StageResult(
            artifacts=[art],
            metrics={"n_samples": geno.n_samples,
                     "mean_distance": round(float(dmatrix[iu].mean()), 6)},
        )
