"""Neighbor-joining tree stage with bootstrap support.

Builds an NJ tree of individuals from the genotype matrix (allele-difference distances) and
annotates internal branches with bootstrap support from site resampling. Emits a Newick
tree; the reporting stage draws it.
"""

from __future__ import annotations

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.phylogenetics.neighbor_joining import bootstrap_support, to_newick
from canidae.stages.popgen.store import allele_difference_matrix, load_genotypes


class NjTreeConfig(StageConfig):
    n_bootstrap: int = 100


@STAGES.register("nj_tree")
class NjTreeStage(Stage):
    name = "nj_tree"
    config_model = NjTreeConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.TREE, "nj")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: NjTreeConfig = self.config  # type: ignore[assignment]
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, "genotypes").path)
        if geno.n_samples < 3:
            raise StageInputError("NJ tree needs >= 3 individuals")

        ac = geno.calls.count_alleles()
        gn = geno.calls.to_n_alt()[ac.is_segregating()]  # (n_sites, n_samples)
        labels = [str(s) for s in geno.samples]

        tree, support = bootstrap_support(
            gn, labels, allele_difference_matrix,
            n_boot=cfg.n_bootstrap, seed=ctx.config.seed)
        newick = to_newick(tree, support)

        out = ctx.datastore.path_for(self.name, "nj_tree.nwk")
        out.write_text(newick + "\n", encoding="utf-8")
        art = ctx.datastore.add(
            ArtifactKind.TREE, "nj", out, fmt=FileFormat.NEWICK, produced_by=self.name,
            metadata={
                "method": "neighbor_joining",
                "n_bootstrap": cfg.n_bootstrap,
                "n_sites": int(gn.shape[0]),
                "n_tips": geno.n_samples,
                "max_support": round(max(support.values()), 3) if support else None,
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={"n_tips": geno.n_samples, "n_sites": int(gn.shape[0])},
        )
