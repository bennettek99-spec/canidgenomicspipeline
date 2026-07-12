"""Demographic-summary stage: site-frequency spectrum and diversity-based estimators.

Per population, computes the folded SFS, Watterson's theta, nucleotide diversity (pi),
genome-wide Tajima's D, and — when a mutation rate is supplied — a diversity-based effective
population size (Ne = pi_per_site / 4*mu). These summaries are the input to, and a sanity
check on, full demographic inference (SMC++, momi2, Stairway Plot), which slot in as
backends. The relative pi/Ne ordering across populations recovers their relative sizes.
"""

from __future__ import annotations

import allel
import numpy as np
import pandas as pd
from pydantic import Field

from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import (
    load_genotypes,
    load_sample_labels,
    population_indices,
)


class DemographyConfig(StageConfig):
    mutation_rate: float = 0.0    # per-bp per-generation; 0 disables Ne
    sequence_length: int = 0      # legacy alias for a verified callable-site denominator
    callable_sites: int = Field(
        default=0,
        ge=0,
        description="Verified callable bases; required for whole-genome pi/theta/Tajima/Ne.",
    )


@STAGES.register("demography")
class DemographyStage(Stage):
    name = "demography"
    config_model = DemographyConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "demography")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: DemographyConfig = self.config  # type: ignore[assignment]
        geno_art = ctx.datastore.get(ArtifactKind.GENOTYPES, "genotypes")
        geno = load_genotypes(geno_art.path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        groups = population_indices(geno, labels)
        callable_sites = _callable_sites(cfg, geno_art.metadata)

        rows, sfs_rows = [], []
        for pop, idx in groups.items():
            ac = geno.calls.count_alleles(subpop=idx)
            n_chrom = 2 * len(idx)
            seg = int(ac.is_segregating().sum())
            a_n = float(np.sum(1.0 / np.arange(1, max(n_chrom, 2))))
            mpd = allel.mean_pairwise_difference(ac, fill=0.0)
            row = {
                "population": pop,
                "n_chromosomes": n_chrom,
                "n_segregating": seg,
                "panel_mean_pairwise_difference": round(float(np.mean(mpd)), 8),
                "statistical_scope": (
                    "per_callable_site" if callable_sites else "panel_relative_only"
                ),
                "callable_sites": callable_sites or None,
            }
            if callable_sites:
                theta_w = (seg / a_n) / callable_sites
                pi = float(mpd.sum()) / callable_sites
                tajd = allel.tajima_d(ac) if seg >= 3 else np.nan
                row.update({
                    "theta_w_per_site": round(theta_w, 8),
                    "pi_per_site": round(pi, 8),
                    "tajima_d": round(float(tajd), 5) if np.isfinite(tajd) else None,
                })
                if cfg.mutation_rate > 0:
                    row["Ne"] = int(pi / (4.0 * cfg.mutation_rate))
            else:
                # Deliberately omit values which would falsely turn a selected-SNP panel
                # into a per-base whole-genome estimate.
                row.update({"theta_w_per_site": None, "pi_per_site": None,
                            "tajima_d": None, "Ne": None})
            rows.append(row)

            folded = allel.sfs_folded(ac)
            for k, count in enumerate(folded):
                sfs_rows.append({"population": pop, "minor_allele_count": k,
                                 "n_sites": int(count)})

        summary = pd.DataFrame(rows)
        stage_dir = ctx.datastore.stage_dir(self.name)
        pd.DataFrame(sfs_rows).to_csv(stage_dir / "sfs.csv", index=False)
        out = stage_dir / "demography_summary.csv"
        summary.to_csv(out, index=False)

        limitations = (
            "Callable-site information unavailable: pi, theta, Tajima's D, and Ne are "
            "not reported as whole-genome quantities."
            if not callable_sites else "Verified callable-site denominator supplied."
        )
        art = Artifact(
            ArtifactKind.ANALYSIS_RESULT, "demography", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "demography", "callable_sites": callable_sites or None,
                      "mutation_rate": cfg.mutation_rate, "limitations": limitations,
                      "scope": "per_callable_site" if callable_sites else "panel_relative_only"})
        return StageResult(
            artifacts=[art],
            metrics={"n_populations": len(summary),
                     "mean_pi": _mean_or_none(summary["pi_per_site"])})


def _callable_sites(cfg: DemographyConfig, metadata: dict) -> int:
    """Use only an explicit callable denominator; SNP span is not one."""
    value = cfg.callable_sites or cfg.sequence_length or metadata.get("callable_sites") or 0
    try:
        return max(int(value), 0)
    except (TypeError, ValueError):
        return 0


def _mean_or_none(values: pd.Series) -> float | None:
    usable = pd.to_numeric(values, errors="coerce").dropna()
    return round(float(usable.mean()), 8) if not usable.empty else None
