"""Read-processing stages: FASTQ -> CRAM -> callset.

Two calling routes share the ``align`` stage:
  * GATK gVCF workflow — ``call_variants`` (per sample) + ``joint_genotype`` (merge); scales
    to large, high-coverage cohorts.
  * ``bcftools_call`` — one multi-sample call over all CRAMs at once; much faster for small
    cohorts (the laptop/red-wolf default).

These wrap standard external tools (bwa-mem2, samtools, GATK/bcftools, GLnexus) that are
Linux/bioconda-native, so the stages run where those binaries are installed; the command
builders are pure-Python and unit-tested. Per-sample work is embarrassingly parallel — the
executor parallelizes across stages, and each tool is given the configured thread count
(tune to your core count; no GPU is used).

Artifacts flow as small manifests (CSV mapping sample -> file) so many-sample cohorts fit
the single-artifact-per-role model.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.runtime import ResourceSpec, ToolSpec
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult

_BASH = ToolSpec(name="bash")


# --------------------------------------------------------------------------------------
# Command builders (pure, testable)
# --------------------------------------------------------------------------------------


def align_command(
    sample: str, fq1: str, fq2: str, ref: str, cram: str, threads: int, *, layout: str = "paired"
) -> str:
    if layout == "paired":
        if not fq2:
            raise StageInputError("paired-end alignment requires fastq2")
        reads = f"{fq1} {fq2}"
    elif layout == "interleaved":
        reads = f"-p {fq1}"
    elif layout == "single":
        reads = fq1
    else:
        raise StageInputError(f"unknown FASTQ layout: {layout}")
    rg = rf"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA"
    return (
        f"bwa-mem2 mem -t {threads} -R '{rg}' {ref} {reads} "
        f"| samtools fixmate -m -u - - "
        f"| samtools sort -@ {threads} -u - "
        f"| samtools markdup -@ {threads} --reference {ref} -O CRAM - {cram} "
        f"&& samtools index {cram}"
    )


def call_command(sample: str, cram: str, ref: str, gvcf: str, backend: str, threads: int) -> str:
    if backend == "gatk":
        return (
            f"gatk HaplotypeCaller -R {ref} -I {cram} -O {gvcf} "
            f"-ERC GVCF --native-pair-hmm-threads {threads}"
        )
    if backend == "bcftools":
        return (
            f"bcftools mpileup -f {ref} --threads {threads} {cram} "
            f"| bcftools call -m -g 5 -Oz -o {gvcf} --threads {threads}"
        )
    if backend == "deepvariant":
        return (
            f"run_deepvariant --model_type=WGS --ref={ref} --reads={cram} "
            f"--output_gvcf={gvcf} --num_shards={threads}"
        )
    raise StageInputError(f"unknown calling backend: {backend}")


def joint_command(gvcfs: list[str], ref: str, out_vcf: str, backend: str, threads: int) -> str:
    if backend == "glnexus":
        joined = " ".join(gvcfs)
        return (
            f"glnexus_cli --config DeepVariantWGS --threads {threads} {joined} "
            f"| bcftools view -Ov -o {out_vcf}"
        )
    if backend == "gatk":
        variants = " ".join(f"-V {g}" for g in gvcfs)
        return (
            f"gatk CombineGVCFs -R {ref} {variants} -O combined.g.vcf.gz "
            f"&& gatk GenotypeGVCFs -R {ref} -V combined.g.vcf.gz -O {out_vcf}"
        )
    raise StageInputError(f"unknown joint-genotyping backend: {backend}")


def bcftools_joint_command(crams: list[str], ref: str, out_vcf: str, threads: int) -> str:
    """Single multi-sample call across all CRAMs at once (fast path for small cohorts)."""
    inputs = " ".join(crams)
    return (
        f"bcftools mpileup -f {ref} -a FORMAT/AD,FORMAT/DP --threads {threads} -Ou {inputs} "
        f"| bcftools call -m -v -Oz -o {out_vcf} --threads {threads} "
        f"&& bcftools index -t {out_vcf}"
    )


def read_fastq_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    required = {"sample_id", "fastq1"}
    if not required <= set(df.columns):
        raise StageInputError(f"fastq manifest needs columns {sorted(required)}")
    if "layout" not in df.columns:
        # Legacy manually-created manifests retain their paired-end interpretation.
        df["layout"] = "paired"
    if "fastq2" not in df.columns:
        df["fastq2"] = ""
    layouts = set(df["layout"].fillna("").astype(str))
    invalid = layouts - {"paired", "interleaved", "single"}
    if invalid:
        raise StageInputError(f"unknown FASTQ layout(s): {sorted(invalid)}")
    paired_missing = (df["layout"] == "paired") & df["fastq2"].fillna("").eq("")
    if paired_missing.any():
        raise StageInputError("paired FASTQ manifest rows need a non-empty fastq2")
    return df


# --------------------------------------------------------------------------------------
# Stages
# --------------------------------------------------------------------------------------


class AlignConfig(StageConfig):
    fastq_manifest: Path = Path("fastq_manifest.csv")
    reference: Path = Path("reference.fasta")
    threads: int = 4


@STAGES.register("align")
class AlignStage(Stage):
    name = "align"
    config_model = AlignConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.RAW_READS, "fastq_manifest", optional=True)]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ALIGNMENT, "crams")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: AlignConfig = self.config  # type: ignore[assignment]
        ctx.runner.ensure(_BASH)
        root = ctx.config.paths.root
        manifest_path = (
            ctx.datastore.get(ArtifactKind.RAW_READS, "fastq_manifest").path
            if ctx.datastore.has(ArtifactKind.RAW_READS, "fastq_manifest")
            else _resolve(cfg.fastq_manifest, root)
        )
        manifest = read_fastq_manifest(manifest_path)
        ref = str(_resolve(cfg.reference, root))
        stage_dir = ctx.datastore.stage_dir(self.name)
        rows = []
        for r in manifest.itertuples(index=False):
            cram = str(stage_dir / f"{r.sample_id}.cram")
            cmd = align_command(
                r.sample_id, r.fastq1, r.fastq2, ref, cram, cfg.threads, layout=r.layout
            )
            ctx.runner.run(
                _BASH,
                ["-lc", cmd],
                resources=ResourceSpec(cpus=cfg.threads),
                record=ctx.scratch.get("_record"),
                cwd=stage_dir,
                expect_outputs=[Path(cram)],
            )
            rows.append({"sample_id": r.sample_id, "cram": cram})
        return _emit_manifest(
            ctx, self.name, rows, "cram_manifest.csv", ArtifactKind.ALIGNMENT, "crams"
        )


class CallConfig(StageConfig):
    reference: Path = Path("reference.fasta")
    backend: str = "gatk"  # gatk | bcftools | deepvariant
    threads: int = 4


@STAGES.register("call_variants")
class CallVariantsStage(Stage):
    name = "call_variants"
    config_model = CallConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ALIGNMENT, "crams")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.GVCF, "gvcfs")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: CallConfig = self.config  # type: ignore[assignment]
        ctx.runner.ensure(_BASH)
        ref = str(_resolve(cfg.reference, ctx.config.paths.root))
        crams = pd.read_csv(ctx.datastore.get(ArtifactKind.ALIGNMENT, "crams").path)
        stage_dir = ctx.datastore.stage_dir(self.name)
        rows = []
        for r in crams.itertuples(index=False):
            gvcf = str(stage_dir / f"{r.sample_id}.g.vcf.gz")
            cmd = call_command(r.sample_id, r.cram, ref, gvcf, cfg.backend, cfg.threads)
            ctx.runner.run(
                _BASH,
                ["-lc", cmd],
                resources=ResourceSpec(cpus=cfg.threads),
                record=ctx.scratch.get("_record"),
                cwd=stage_dir,
                expect_outputs=[Path(gvcf)],
            )
            rows.append({"sample_id": r.sample_id, "gvcf": gvcf})
        return _emit_manifest(ctx, self.name, rows, "gvcf_manifest.csv", ArtifactKind.GVCF, "gvcfs")


class JointConfig(StageConfig):
    reference: Path = Path("reference.fasta")
    backend: str = "glnexus"  # glnexus | gatk
    threads: int = 4


@STAGES.register("joint_genotype")
class JointGenotypeStage(Stage):
    name = "joint_genotype"
    config_model = JointConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.GVCF, "gvcfs")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.CALLSET, "callset")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: JointConfig = self.config  # type: ignore[assignment]
        ctx.runner.ensure(_BASH)
        ref = str(_resolve(cfg.reference, ctx.config.paths.root))
        gvcfs = pd.read_csv(ctx.datastore.get(ArtifactKind.GVCF, "gvcfs").path)
        stage_dir = ctx.datastore.stage_dir(self.name)
        out_vcf = str(stage_dir / "joint.vcf")
        cmd = joint_command(list(gvcfs["gvcf"]), ref, out_vcf, cfg.backend, cfg.threads)
        ctx.runner.run(
            _BASH,
            ["-lc", cmd],
            resources=ResourceSpec(cpus=cfg.threads),
            record=ctx.scratch.get("_record"),
            cwd=stage_dir,
            expect_outputs=[Path(out_vcf)],
        )
        art = ctx.datastore.add(
            ArtifactKind.CALLSET,
            "callset",
            Path(out_vcf),
            fmt=FileFormat.VCF,
            produced_by=self.name,
            metadata={"backend": cfg.backend, "n_samples": len(gvcfs)},
        )
        return StageResult(artifacts=[art], metrics={"n_samples": len(gvcfs)})


class BcftoolsCallConfig(StageConfig):
    reference: Path = Path("reference.fasta")
    threads: int = 4


@STAGES.register("bcftools_call")
class BcftoolsCallStage(Stage):
    """Fast path: one multi-sample ``bcftools`` call over all CRAMs -> a joint callset.

    Replaces the per-sample ``call_variants`` + ``joint_genotype`` pair (GATK gVCF workflow)
    with a single command. Much faster for small cohorts; the GATK path scales better and
    gives higher-quality calls for large, high-coverage cohorts.
    """

    name = "bcftools_call"
    config_model = BcftoolsCallConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ALIGNMENT, "crams")]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.CALLSET, "callset")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: BcftoolsCallConfig = self.config  # type: ignore[assignment]
        ctx.runner.ensure(_BASH)
        ref = str(_resolve(cfg.reference, ctx.config.paths.root))
        crams = pd.read_csv(ctx.datastore.get(ArtifactKind.ALIGNMENT, "crams").path)
        stage_dir = ctx.datastore.stage_dir(self.name)
        out_vcf = str(stage_dir / "joint.vcf.gz")
        cmd = bcftools_joint_command(list(crams["cram"]), ref, out_vcf, cfg.threads)
        ctx.runner.run(
            _BASH,
            ["-lc", cmd],
            resources=ResourceSpec(cpus=cfg.threads),
            record=ctx.scratch.get("_record"),
            cwd=stage_dir,
            expect_outputs=[Path(out_vcf)],
        )
        art = ctx.datastore.add(
            ArtifactKind.CALLSET,
            "callset",
            Path(out_vcf),
            fmt=FileFormat.VCF,
            produced_by=self.name,
            metadata={"caller": "bcftools", "n_samples": len(crams)},
        )
        return StageResult(artifacts=[art], metrics={"n_samples": len(crams)})


def _resolve(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else (root / path)


def _emit_manifest(
    ctx: RunContext, stage: str, rows: list[dict], filename: str, kind: ArtifactKind, role: str
) -> StageResult:
    df = pd.DataFrame(rows)
    out = ctx.datastore.path_for(stage, filename)
    df.to_csv(out, index=False)
    art = ctx.datastore.add(
        kind, role, out, fmt=FileFormat.CSV, produced_by=stage, metadata={"n_samples": len(df)}
    )
    return StageResult(artifacts=[art], metrics={"n_samples": len(df)})
