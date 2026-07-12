"""TreeMix stage: population tree with migration edges (TreeMix backend).

TreeMix infers a maximum-likelihood population tree from allele frequencies and adds
migration edges to model gene flow. It is an external tool (bioconda, Linux); this stage
prepares its gzipped allele-count input (pure-Python, unit-tested) and runs it where the
``treemix`` binary is available. Not in the default example pipeline.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.runtime import ResourceSpec, ToolSpec
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import (
    Genotypes,
    load_genotypes,
    load_sample_labels,
    population_indices,
)

_TREEMIX = ToolSpec(name="treemix", version_args=("--version",))


class TreemixConfig(StageConfig):
    n_migrations: int = 2
    block_size: int = 500       # SNPs per block (-k), for the covariance jackknife
    root: str = ""              # optional outgroup population to root the tree


@STAGES.register("treemix")
class TreemixStage(Stage):
    name = "treemix"
    config_model = TreemixConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.TREE, "treemix")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: TreemixConfig = self.config  # type: ignore[assignment]
        ctx.runner.ensure(_TREEMIX)
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, "genotypes").path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        groups = population_indices(geno, labels)
        if len(groups) < 3:
            raise StageInputError("TreeMix needs >= 3 populations")

        stage_dir = ctx.datastore.stage_dir(self.name)
        infile = stage_dir / "treemix_input.gz"
        write_treemix_input(geno, groups, infile)

        args = ["-i", infile.name, "-m", str(cfg.n_migrations), "-k", str(cfg.block_size),
                "-o", "treemix_out"]
        if cfg.root:
            args += ["-root", cfg.root]
        ctx.runner.run(
            _TREEMIX, args, resources=ResourceSpec(cpus=ctx.config.resources.cpus),
            record=ctx.scratch.get("_record"), cwd=stage_dir,
            expect_outputs=[stage_dir / "treemix_out.treeout.gz"])

        with gzip.open(stage_dir / "treemix_out.treeout.gz", "rt") as fh:
            newick = fh.readline().strip()
        out = ctx.datastore.path_for(self.name, "treemix_tree.nwk")
        out.write_text(newick + "\n", encoding="utf-8")
        art = ctx.datastore.add(
            ArtifactKind.TREE, "treemix", out, fmt=FileFormat.NEWICK, produced_by=self.name,
            metadata={"method": "treemix", "n_migrations": cfg.n_migrations})
        return StageResult(artifacts=[art], metrics={"n_populations": len(groups)})


def write_treemix_input(geno: Genotypes, groups: dict[str, list[int]], path: Path) -> Path:
    """Write TreeMix's gzipped input: header of population names, then one line per SNP with
    ``ref,alt`` allele counts per population."""
    pops = list(groups)
    counts = {}
    for pop, idx in groups.items():
        ac = geno.calls.count_alleles(subpop=idx)
        counts[pop] = np.stack([ac[:, 0], ac[:, 1] if ac.shape[1] > 1
                                else np.zeros(ac.shape[0], int)], axis=1)
    n_sites = geno.n_variants
    with gzip.open(path, "wt") as fh:
        fh.write(" ".join(pops) + "\n")
        for i in range(n_sites):
            fh.write(" ".join(f"{counts[p][i, 0]},{counts[p][i, 1]}" for p in pops) + "\n")
    return path
