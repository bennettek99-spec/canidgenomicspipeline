"""Pairwise F_ST between populations.

Uses Hudson's estimator (Bhatia et al. 2013), which is robust to unequal sample sizes: it
sums per-site numerators and denominators across all sites before taking the ratio. Output
is a symmetric population-by-population matrix.
"""

from __future__ import annotations

import allel
import numpy as np
import pandas as pd

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import (
    load_genotypes,
    load_sample_labels,
    population_indices,
)


class FstConfig(StageConfig):
    method: str = "hudson"
    min_samples_per_population: int = 1


@STAGES.register("fst")
class FstStage(Stage):
    name = "fst"
    config_model = FstConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "fst")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: FstConfig = self.config  # type: ignore[assignment]
        if cfg.method != "hudson":
            raise StageInputError(f"unsupported F_ST method: {cfg.method}")

        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, "genotypes").path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        groups = {
            pop: idx
            for pop, idx in population_indices(geno, labels).items()
            if len(idx) >= cfg.min_samples_per_population
        }
        if len(groups) < 2:
            raise StageInputError(
                f"F_ST needs >=2 populations with >={cfg.min_samples_per_population} "
                f"sample(s); got {list(groups)}"
            )

        pops = list(groups)
        allele_counts = {p: geno.calls.count_alleles(subpop=groups[p]) for p in pops}
        matrix = pd.DataFrame(0.0, index=pops, columns=pops)
        pair_values: dict[str, float] = {}
        for i, a in enumerate(pops):
            for b in pops[i + 1:]:
                fst = _hudson_fst(allele_counts[a], allele_counts[b])
                matrix.loc[a, b] = matrix.loc[b, a] = round(fst, 6)
                pair_values[f"{a}__{b}"] = round(fst, 6)

        out = ctx.datastore.path_for(self.name, "fst_matrix.csv")
        matrix.to_csv(out)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "fst", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "fst", "method": "hudson", "populations": pops,
                      "pairwise": pair_values},
        )
        finite = [v for v in pair_values.values() if np.isfinite(v)]
        return StageResult(
            artifacts=[art],
            metrics={"n_populations": len(pops),
                     "max_fst": max(finite) if finite else None,
                     "min_fst": min(finite) if finite else None},
        )


def _hudson_fst(ac1: allel.AlleleCountsArray, ac2: allel.AlleleCountsArray) -> float:
    num, den = allel.hudson_fst(ac1, ac2)
    denom = np.nansum(den)
    if denom <= 0:
        return float("nan")
    return float(np.nansum(num) / denom)
