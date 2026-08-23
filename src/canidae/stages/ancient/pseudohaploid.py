"""Ancient-DNA support: pseudohaploidization.

Low-coverage / ancient genomes are typically represented as *pseudohaploid* calls — one
randomly sampled read allele per site — to avoid false heterozygosity from damage and low
depth. This stage draws one allele per called genotype and emits a homozygous
(pseudohaploid) genotype set usable by the downstream analyses. Damage rescaling
(mapDamage) and genotype-likelihood workflows (ANGSD) are complementary backends for the
read-level steps; this operates on an existing genotype matrix.
"""

from __future__ import annotations

import allel
import numpy as np
import pandas as pd

from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import Genotypes, load_genotypes, save_genotypes


class AncientConfig(StageConfig):
    pass


@STAGES.register("pseudohaploid")
class PseudohaploidStage(Stage):
    name = "pseudohaploid"
    config_model = AncientConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "pseudohaploid"),
            ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "ancient"),
        ]

    def run(self, ctx: RunContext) -> StageResult:
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, "genotypes").path)
        gt = np.asarray(geno.calls)  # (n_sites, n_samples, 2)
        called = (gt[:, :, 0] >= 0) & (gt[:, :, 1] >= 0)

        rng = np.random.default_rng(ctx.config.seed)
        pick = rng.integers(0, 2, size=called.shape)
        chosen = np.where(pick == 0, gt[:, :, 0], gt[:, :, 1])
        pseudo = np.stack([chosen, chosen], axis=-1).astype(np.int8)
        pseudo[~called] = -1

        het_before = geno.calls.is_het().sum(axis=0)
        pseudo_geno = Genotypes(
            calls=allel.GenotypeArray(pseudo), pos=geno.pos, chrom=geno.chrom, samples=geno.samples
        )
        het_after = pseudo_geno.calls.is_het().sum(axis=0)

        stage_dir = ctx.datastore.stage_dir(self.name)
        npz = save_genotypes(stage_dir / "pseudohaploid.npz", pseudo_geno)
        summary = pd.DataFrame(
            {
                "sample_id": geno.samples,
                "n_called": called.sum(axis=0).astype(int),
                "het_before": het_before.astype(int),
                "het_after": het_after.astype(int),
            }
        )
        out = stage_dir / "pseudohaploid_summary.csv"
        summary.to_csv(out, index=False)

        geno_art = ctx.datastore.add(
            ArtifactKind.GENOTYPES,
            "pseudohaploid",
            npz,
            fmt=FileFormat.NPZ,
            produced_by=self.name,
            metadata={"n_samples": geno.n_samples, "n_variants": geno.n_variants},
        )
        summary_art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT,
            "ancient",
            out,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "pseudohaploid"},
        )
        return StageResult(
            artifacts=[geno_art, summary_art],
            metrics={"n_samples": geno.n_samples, "total_het_after": int(het_after.sum())},
        )
