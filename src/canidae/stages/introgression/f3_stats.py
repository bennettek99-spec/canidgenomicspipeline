"""Outgroup-f3 stage: a shared-drift matrix among populations.

Outgroup f3(O; A, B) measures the shared genetic drift of A and B since their divergence
from an outgroup O — a robust summary of relatedness widely used to build affinity
matrices. Higher values mean more shared history.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.introgression import fstats
from canidae.stages.introgression.dstats import _auto_outgroup
from canidae.stages.popgen.store import (
    load_genotypes,
    load_sample_labels,
    population_indices,
)


class F3Config(StageConfig):
    outgroup: str = ""
    allow_auto_outgroup: bool = False
    n_blocks: int = 20


@STAGES.register("f3")
class F3Stage(Stage):
    name = "f3"
    config_model = F3Config

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "f3")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: F3Config = self.config  # type: ignore[assignment]
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        groups = population_indices(geno, labels)
        freqs = fstats.allele_frequencies(geno, groups)

        if not cfg.outgroup and not cfg.allow_auto_outgroup:
            raise StageInputError(
                "outgroup-f3 requires an explicit biological outgroup; set stages.f3.outgroup"
            )
        outgroup = cfg.outgroup or _auto_outgroup(freqs)
        ingroup = [p for p in groups if p != outgroup]
        if len(ingroup) < 2:
            raise StageInputError("outgroup-f3 needs >= 2 non-outgroup populations")

        matrix = pd.DataFrame(np.nan, index=ingroup, columns=ingroup)
        for i, a in enumerate(ingroup):
            for b in ingroup[i:]:
                val = fstats.outgroup_f3(freqs[outgroup], freqs[a], freqs[b],
                                         n_blocks=cfg.n_blocks).estimate
                matrix.loc[a, b] = matrix.loc[b, a] = round(val, 6)

        out = ctx.datastore.path_for(self.name, "outgroup_f3.csv")
        matrix.to_csv(out)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "f3", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "outgroup_f3", "outgroup": outgroup,
                      "populations": ingroup})
        return StageResult(
            artifacts=[art],
            metrics={"outgroup": outgroup, "n_populations": len(ingroup)})
