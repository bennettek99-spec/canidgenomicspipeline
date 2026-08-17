"""Principal component analysis of the genotype matrix.

Uses scikit-allel's Patterson-scaled PCA (the standard for SNP data) over segregating
biallelic sites. Emits per-sample PC coordinates joined to taxon/population labels, plus the
explained-variance ratio (carried on the artifact metadata for the report).
"""

from __future__ import annotations

import allel
import numpy as np
import pandas as pd

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import align_labels, load_genotypes, load_sample_labels


class PCAConfig(StageConfig):
    n_components: int = 6
    scaler: str = "patterson"  # "patterson" | "standard" | None


@STAGES.register("pca")
class PCAStage(Stage):
    name = "pca"
    config_model = PCAConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "pca")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: PCAConfig = self.config  # type: ignore[assignment]
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        labels = align_labels(
            geno, load_sample_labels(
                ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        )

        # Restrict to segregating sites and retain ALT allele counts (n_variants, n_samples).
        ac = geno.calls.count_alleles()
        seg = ac.is_segregating()
        n_alt = geno.calls.to_n_alt()[seg]
        if n_alt.shape[0] < 2:
            raise StageInputError("PCA needs at least 2 segregating sites")

        max_components = int(min(cfg.n_components, geno.n_samples - 1, n_alt.shape[0] - 1))
        if max_components < 2:
            raise StageInputError(
                f"too few samples/sites for PCA (got n_components={max_components})"
            )

        scaler = None if cfg.scaler.lower() in {"none", ""} else cfg.scaler
        # SciPy is removing float16 SVD support; float32 is stable and still laptop-friendly.
        coords, model = allel.pca(
            np.asarray(n_alt, dtype=np.float32), n_components=max_components, scaler=scaler
        )

        pc_cols = [f"PC{i + 1}" for i in range(max_components)]
        table = pd.DataFrame(coords, columns=pc_cols)
        table.insert(0, "population", labels["population"].to_numpy())
        table.insert(0, "taxon", labels["taxon"].to_numpy())
        table.insert(0, "sample_id", geno.samples)

        out = ctx.datastore.path_for(self.name, "pca_coords.csv")
        table.to_csv(out, index=False)

        evr = [round(float(x), 5) for x in model.explained_variance_ratio_]
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "pca", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "pca",
                "n_components": max_components,
                "n_sites_used": int(n_alt.shape[0]),
                "explained_variance_ratio": evr,
                "pc_columns": pc_cols,
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={"n_components": max_components, "n_sites_used": int(n_alt.shape[0]),
                     "pc1_variance": evr[0] if evr else None},
        )
