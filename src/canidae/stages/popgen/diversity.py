"""Within-population genetic-diversity summaries.

For each population: nucleotide diversity (mean per-site pairwise difference, π), the number
of segregating sites, and mean observed heterozygosity across its samples. These are the
core comparative-diversity statistics; windowed/accessible-genome-scaled estimates and
between-population dxy are natural extensions.
"""

from __future__ import annotations

import allel
import numpy as np
import pandas as pd
from pydantic import Field

from canidae.core.errors import StageInputError
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import (
    load_genotypes,
    load_sample_labels,
    population_indices,
)


class DiversityConfig(StageConfig):
    min_samples_per_population: int = 2
    callable_sites: int = Field(
        default=0,
        ge=0,
        description="Verified callable bases; 0 keeps results explicitly panel-relative.",
    )


@STAGES.register("diversity")
class DiversityStage(Stage):
    name = "diversity"
    config_model = DiversityConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "diversity")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: DiversityConfig = self.config  # type: ignore[assignment]
        role = (
            "analysis_genotypes"
            if ctx.datastore.has(ArtifactKind.GENOTYPES, "analysis_genotypes")
            else "genotypes"
        )
        geno_art = ctx.datastore.get(ArtifactKind.GENOTYPES, role)
        geno = load_genotypes(geno_art.path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path
        )
        groups = {
            pop: idx
            for pop, idx in population_indices(geno, labels).items()
            if len(idx) >= cfg.min_samples_per_population
        }
        if not groups:
            raise StageInputError("no population met the minimum sample threshold")

        # Per-sample observed heterozygosity, computed once over all sites.
        is_het = geno.calls.is_het()
        is_called = geno.calls.is_called()
        sample_het = is_het.sum(axis=0) / np.maximum(is_called.sum(axis=0), 1)

        callable_sites = _callable_sites(cfg, geno_art.metadata)
        scope = "per_callable_site" if callable_sites else "panel_relative"
        rows = []
        for pop, idx in groups.items():
            ac = geno.calls.count_alleles(subpop=idx)
            mpd = allel.mean_pairwise_difference(ac, fill=0.0)
            panel_diversity = float(np.mean(mpd))
            diversity = float(mpd.sum() / callable_sites) if callable_sites else panel_diversity
            rows.append(
                {
                    "population": pop,
                    "n_samples": len(idx),
                    "n_segregating_sites": int(ac.is_segregating().sum()),
                    # ``pi`` remains for backwards-compatible table consumers; the explicit
                    # scope fields prevent a selected-SNP average being mistaken for genome-wide pi.
                    "pi": round(diversity, 8),
                    "panel_relative_diversity": round(panel_diversity, 8),
                    "diversity_scope": scope,
                    "callable_sites": callable_sites or None,
                    "mean_heterozygosity": round(float(np.mean(sample_het[idx])), 6),
                    "heterozygosity_scope": "called_panel_sites"
                    if not callable_sites
                    else "called_sites_with_callable_denominator",
                }
            )
        table = pd.DataFrame(rows).sort_values("population").reset_index(drop=True)

        out = ctx.datastore.path_for(self.name, "diversity.csv")
        table.to_csv(out, index=False)
        limitations = (
            "Panel-relative heterozygosity and diversity only; selected SNPs are not a "
            "callable-genome denominator."
            if not callable_sites
            else "Callable-site denominator supplied by configuration/source metadata."
        )
        small_groups = {
            pop: len(idx) for pop, idx in population_indices(geno, labels).items() if len(idx) < 2
        }
        if small_groups:
            limitations += (
                " Singleton population estimates are descriptive and do not represent "
                "population-level uncertainty."
            )
        art = Artifact(
            ArtifactKind.ANALYSIS_RESULT,
            "diversity",
            out,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "diversity",
                "populations": table["population"].tolist(),
                "scope": scope,
                "callable_sites": callable_sites or None,
                "limitations": limitations,
                "small_sample_populations": small_groups,
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={"n_populations": len(table), "mean_pi": round(float(table["pi"].mean()), 6)},
        )


def _callable_sites(cfg: DiversityConfig, metadata: dict) -> int:
    """Resolve an auditable callable-site denominator, never infer one from SNP spacing."""
    value = cfg.callable_sites or metadata.get("callable_sites") or 0
    try:
        return max(int(value), 0)
    except (TypeError, ValueError):
        return 0
