"""Structural-variation calling stage (Delly backend).

Wraps Delly to call deletions/duplications/inversions from the aligned CRAMs. Delly is an
external Linux/bioconda tool, so this stage runs where it is installed; the command builder
is pure-Python and unit-tested. A minimal placeholder for the SV module — Manta/GRIDSS and
genotyping/merging across samples (bcftools/SURVIVOR) are follow-ups behind the same stage.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pydantic import Field

from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.runtime import ResourceSpec, ToolSpec
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult

_BASH = ToolSpec(name="bash")


class SvConfig(StageConfig):
    reference: Path = Path("reference.fasta")
    sv_types: list[str] = Field(default_factory=lambda: ["DEL", "DUP", "INV"])
    threads: int = 4


@STAGES.register("sv")
class StructuralVariationStage(Stage):
    name = "sv"
    config_model = SvConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ALIGNMENT, "crams")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "sv")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: SvConfig = self.config  # type: ignore[assignment]
        ctx.runner.ensure(_BASH)
        ref = str(_resolve(cfg.reference, ctx.config.paths.root))
        crams = pd.read_csv(ctx.datastore.get(ArtifactKind.ALIGNMENT, "crams").path)
        stage_dir = ctx.datastore.stage_dir(self.name)
        out_vcf = str(stage_dir / "sv.vcf")
        cmd = sv_command(list(crams["cram"]), ref, out_vcf)
        ctx.runner.run(
            _BASH,
            ["-lc", cmd],
            resources=ResourceSpec(cpus=cfg.threads),
            record=ctx.scratch.get("_record"),
            cwd=stage_dir,
            expect_outputs=[Path(out_vcf)],
        )
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT,
            "sv",
            Path(out_vcf),
            fmt=FileFormat.VCF,
            produced_by=self.name,
            metadata={"caller": "delly", "sv_types": cfg.sv_types},
        )
        return StageResult(artifacts=[art], metrics={"n_samples": len(crams)})


def sv_command(crams: list[str], ref: str, out_vcf: str) -> str:
    inputs = " ".join(crams)
    return f"delly call -g {ref} -o sv.bcf {inputs} && bcftools view sv.bcf -Ov -o {out_vcf}"


def _resolve(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else (root / path)
