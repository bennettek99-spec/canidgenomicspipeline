"""Ingest stage: bring a published callset + its sample sheet into the workspace.

This is the Phase-1 entry point that lets analyses start from an already-published VCF
(e.g. a Dog10K-derived callset) without running read alignment/calling first. It validates
the sample sheet, registers the (possibly large, referenced-in-place) VCF as a ``CALLSET``
artifact, and emits a normalized ``SAMPLE_SHEET`` artifact that downstream stages read for
per-sample metadata.

Full remote acquisition (SRA/ENA/Dog10K fetchers, integrity, versioning) plugs in here
behind the same stage as additional fetchers; this local-path ingest is the first one.
"""

from __future__ import annotations

from pathlib import Path

from pydantic import Field

from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.acquisition.sample_sheet import read_sample_sheet, to_cohort


class IngestConfig(StageConfig):
    """Configuration for :class:`IngestStage`."""

    sample_sheet: Path = Field(..., description="CSV of per-sample metadata.")
    callset: Path = Field(..., description="Published VCF/BCF to analyze.")
    dataset_id: str = "local"
    cohort_id: str = "cohort"
    reference_build: str = "unknown"
    panel_relative: bool = False
    panel_name: str = ""
    callable_sites: int = 0
    analysis_limitations: list[str] = Field(
        default_factory=list,
        description="Study-specific limitations displayed in the final report.",
    )


@STAGES.register("ingest")
class IngestStage(Stage):
    name = "ingest"
    config_model = IngestConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []  # entry point: reads from local paths named in config

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
            ArtifactSpec(ArtifactKind.CALLSET, "callset"),
        ]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: IngestConfig = self.config  # type: ignore[assignment]
        root = ctx.config.paths.root
        sheet_path = _resolve(cfg.sample_sheet, root)
        vcf_path = _resolve(cfg.callset, root)
        if not vcf_path.exists():
            from canidae.core.errors import StageInputError

            raise StageInputError(f"callset not found: {vcf_path}")

        df = read_sample_sheet(sheet_path)
        cohort = to_cohort(df, cohort_id=cfg.cohort_id, dataset_id=cfg.dataset_id)

        normalized = ctx.datastore.path_for(self.name, "sample_sheet.csv")
        df.to_csv(normalized, index=False)

        sheet_art = ctx.datastore.add(
            ArtifactKind.SAMPLE_SHEET, "sample_sheet", normalized,
            fmt=FileFormat.CSV, produced_by=self.name,
            metadata={"n_samples": len(df)},
        )
        callset_art = ctx.datastore.add(
            ArtifactKind.CALLSET, "callset", vcf_path,
            fmt=_guess_format(vcf_path), produced_by=self.name,
            metadata={
                "dataset_id": cfg.dataset_id,
                "source": str(vcf_path),
                "reference_build": cfg.reference_build,
                "panel_relative": cfg.panel_relative,
                "panel_name": cfg.panel_name or None,
                "callable_sites": cfg.callable_sites or None,
                "analysis_limitations": cfg.analysis_limitations,
            },
        )

        metrics = {
            "n_samples": cohort.size,
            "n_populations": len(cohort.populations),
            "taxa": {t.value: n for t, n in cohort.taxa.items()},
        }
        return StageResult(artifacts=[sheet_art, callset_art], metrics=metrics)


class SampleSheetConfig(StageConfig):
    """Configuration for :class:`SampleSheetStage`."""

    path: Path = Field(..., description="CSV of per-sample metadata.")


@STAGES.register("sample_sheet")
class SampleSheetStage(Stage):
    """Register a validated sample sheet without a callset.

    Used by the read-processing path, where the callset is produced by
    ``joint_genotype`` and only the per-sample metadata needs to be brought in for the
    downstream population-genomics stages.
    """

    name = "sample_sheet"
    config_model = SampleSheetConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: SampleSheetConfig = self.config  # type: ignore[assignment]
        df = read_sample_sheet(_resolve(cfg.path, ctx.config.paths.root))
        normalized = ctx.datastore.path_for(self.name, "sample_sheet.csv")
        df.to_csv(normalized, index=False)
        art = ctx.datastore.add(
            ArtifactKind.SAMPLE_SHEET, "sample_sheet", normalized, fmt=FileFormat.CSV,
            produced_by=self.name, metadata={"n_samples": len(df)})
        return StageResult(artifacts=[art], metrics={"n_samples": len(df)})


def _resolve(path: Path, root: Path) -> Path:
    return path if path.is_absolute() else (root / path)


def _guess_format(path: Path) -> FileFormat:
    name = path.name.lower()
    if name.endswith((".bcf",)):
        return FileFormat.BCF
    if name.endswith((".vcf", ".vcf.gz")):
        return FileFormat.VCF
    return FileFormat.OTHER
