"""Maximum-likelihood tree stage (IQ-TREE 2 backend).

Converts the genotype matrix to an IUPAC nucleotide alignment and runs IQ-TREE 2 with model
selection and ultrafast bootstrap. IQ-TREE is an external tool (best via bioconda on Linux),
so this stage runs where the ``iqtree2`` binary is available; the alignment writer itself is
pure-Python and unit-tested. Not included in the default example pipeline.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.runtime import ResourceSpec, ToolSpec
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import Genotypes, load_genotypes

_IQTREE = ToolSpec(name="iqtree2", version_args=("--version",))
# ALT allele count -> IUPAC: 0=hom ref (A), 2=hom alt (T), 1=het (W = A/T), missing=N
_IUPAC = np.array(["A", "W", "T"], dtype="U1")


class MlTreeConfig(StageConfig):
    model: str = "MFP"          # ModelFinder Plus
    ufboot: int = 1000          # ultrafast bootstrap replicates
    segregating_only: bool = True


@STAGES.register("ml_tree")
class MlTreeStage(Stage):
    name = "ml_tree"
    config_model = MlTreeConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.TREE, "ml")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: MlTreeConfig = self.config  # type: ignore[assignment]
        ctx.runner.ensure(_IQTREE)
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        stage_dir = ctx.datastore.stage_dir(self.name)
        aln = stage_dir / "alignment.phy"
        write_phylip(geno, aln, segregating_only=cfg.segregating_only)

        ctx.runner.run(
            _IQTREE,
            ["-s", aln.name, "-m", cfg.model, "-B", str(cfg.ufboot), "-redo"],
            resources=ResourceSpec(cpus=ctx.config.resources.cpus),
            record=ctx.scratch.get("_record"), cwd=stage_dir,
            expect_outputs=[stage_dir / "alignment.phy.treefile"],
        )
        treefile = stage_dir / "alignment.phy.treefile"
        out = ctx.datastore.path_for(self.name, "ml_tree.nwk")
        out.write_text(treefile.read_text(encoding="utf-8"), encoding="utf-8")
        art = ctx.datastore.add(
            ArtifactKind.TREE, "ml", out, fmt=FileFormat.NEWICK, produced_by=self.name,
            metadata={"method": "iqtree2", "model": cfg.model, "ufboot": cfg.ufboot})
        return StageResult(artifacts=[art], metrics={"n_tips": geno.n_samples})


def write_phylip(geno: Genotypes, path: Path, *, segregating_only: bool = True) -> Path:
    """Write a relaxed-PHYLIP IUPAC alignment (rows = individuals, columns = SNP sites)."""
    n_alt = np.asarray(geno.calls.to_n_alt(fill=-1))  # (n_sites, n_samples)
    if segregating_only:
        seg = geno.calls.count_alleles().is_segregating()
        n_alt = n_alt[seg]
    n_sites, n_samples = n_alt.shape
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"{n_samples} {n_sites}\n")
        for s, name in enumerate(geno.samples):
            col = n_alt[:, s]
            chars = np.where(col < 0, "N", _IUPAC[np.clip(col, 0, 2)])
            fh.write(f"{name}  {''.join(chars.tolist())}\n")
    return path
