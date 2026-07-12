"""Enforced, auditable callset quality control.

QC is deliberately an *analysis gate*, not a dashboard that downstream stages may ignore.
It evaluates sample and site call rates, removes failed records, and emits a compact filtered
VCF plus a filtered sample sheet.  The unfiltered source is never overwritten: the exclusion
table records every removal and its exact reason, while ``load_genotypes`` preferentially
consumes the ``qc_callset`` artifact.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import Field

from canidae.core.errors import StageInputError
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.acquisition.sample_sheet import read_sample_sheet
from canidae.stages.processing.vcf_io import read_biallelic_snps, write_minimal_vcf


class SampleQCConfig(StageConfig):
    """Thresholds for the mandatory sample- and site-level QC gate.

    ``min_call_rate`` remains as a backwards-compatible alias for the sample threshold.
    The default site threshold is intentionally stricter: a poor site affects every sample
    and should not survive simply because a few high-quality samples do.
    """

    min_call_rate: float = Field(default=0.5, ge=0.0, le=1.0)
    min_sample_call_rate: float | None = Field(default=None, ge=0.0, le=1.0)
    min_site_call_rate: float = Field(default=0.8, ge=0.0, le=1.0)
    min_site_maf: float = Field(default=0.0, ge=0.0, le=0.5)

    @property
    def sample_threshold(self) -> float:
        return (
            self.min_call_rate
            if self.min_sample_call_rate is None
            else self.min_sample_call_rate
        )


@STAGES.register("qc")
class SampleQCStage(Stage):
    """Filter failed samples/sites and retain an auditable exclusion report."""

    name = "qc"
    config_model = SampleQCConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.CALLSET, "callset"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.QC_TABLE, "sample_qc"),
            ArtifactSpec(ArtifactKind.QC_TABLE, "site_qc"),
            ArtifactSpec(ArtifactKind.QC_TABLE, "qc_exclusions"),
            ArtifactSpec(ArtifactKind.CALLSET, "qc_callset"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "qc_sample_sheet"),
        ]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: SampleQCConfig = self.config  # type: ignore[assignment]
        source = ctx.datastore.get(ArtifactKind.CALLSET, "callset")
        sheet = ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet")
        sites = read_biallelic_snps(source.path)
        labels = read_sample_sheet(sheet.path)
        if not len(sites.samples) or not len(sites.pos):
            raise StageInputError("QC needs at least one sample and one biallelic SNP")

        sample_table, sample_keep = _sample_metrics(
            sites.gt, sites.samples, set(labels["sample_id"].astype(str)), cfg.sample_threshold
        )
        if not sample_keep.any():
            raise StageInputError(
                "all samples failed QC; lower min_sample_call_rate or inspect qc_exclusions.csv"
            )

        retained_gt = sites.gt[:, sample_keep]
        site_table, site_keep = _site_metrics(
            retained_gt, sites.chrom, sites.pos, sites.ref, sites.alt, cfg
        )
        if not site_keep.any():
            raise StageInputError(
                "all sites failed QC; lower min_site_call_rate/min_site_maf or inspect "
                "qc_exclusions.csv"
            )

        retained_samples = sites.samples[sample_keep]
        retained_sheet = labels[labels["sample_id"].astype(str).isin(set(retained_samples))].copy()
        # Preserve VCF column order rather than the order inherited from an arbitrary CSV.
        retained_sheet["_qc_order"] = retained_sheet["sample_id"].map(
            {str(sample): i for i, sample in enumerate(retained_samples)}
        )
        retained_sheet = retained_sheet.sort_values("_qc_order").drop(columns="_qc_order")

        exclusions = pd.concat(
            [_sample_exclusions(sample_table), _site_exclusions(site_table)], ignore_index=True
        )
        stage_dir = ctx.datastore.stage_dir(self.name)
        callset_out = stage_dir / "qc_filtered.vcf"
        _write_filtered_vcf(
            callset_out,
            sites.chrom[site_keep], sites.pos[site_keep], sites.ref[site_keep],
            sites.alt[site_keep],
            retained_samples, retained_gt[site_keep],
        )
        sheet_out = _atomic_csv(retained_sheet, stage_dir / "qc_sample_sheet.csv")
        sample_out = _atomic_csv(sample_table, stage_dir / "sample_qc.csv")
        site_out = _atomic_csv(site_table, stage_dir / "site_qc.csv")
        exclusions_out = _atomic_csv(exclusions, stage_dir / "qc_exclusions.csv")

        source_metadata = dict(source.metadata)
        filtered_metadata = {
            **source_metadata,
            "qc_enforced": True,
            "n_samples_input": len(sites.samples),
            "n_samples_retained": int(sample_keep.sum()),
            "n_sites_input": len(sites.pos),
            "n_sites_retained": int(site_keep.sum()),
            "min_sample_call_rate": cfg.sample_threshold,
            "min_site_call_rate": cfg.min_site_call_rate,
            "min_site_maf": cfg.min_site_maf,
        }
        artifacts = [
            Artifact(ArtifactKind.QC_TABLE, "sample_qc", sample_out, FileFormat.CSV,
                     produced_by=self.name),
            Artifact(ArtifactKind.QC_TABLE, "site_qc", site_out, FileFormat.CSV,
                     produced_by=self.name),
            Artifact(ArtifactKind.QC_TABLE, "qc_exclusions", exclusions_out, FileFormat.CSV,
                     produced_by=self.name),
            Artifact(ArtifactKind.CALLSET, "qc_callset", callset_out, FileFormat.VCF,
                     produced_by=self.name, metadata=filtered_metadata),
            Artifact(ArtifactKind.SAMPLE_SHEET, "qc_sample_sheet", sheet_out, FileFormat.CSV,
                     produced_by=self.name,
                     metadata={"n_samples": len(retained_sheet), "qc_enforced": True}),
        ]
        return StageResult(
            artifacts=artifacts,
            metrics={
                "n_samples_input": len(sites.samples),
                "n_samples_retained": int(sample_keep.sum()),
                "n_samples_excluded": int((~sample_keep).sum()),
                "n_sites_input": len(sites.pos),
                "n_sites_retained": int(site_keep.sum()),
                "n_sites_excluded": int((~site_keep).sum()),
                "exclusion_report": str(exclusions_out),
            },
        )


def _sample_metrics(
    gt: np.ndarray, samples: np.ndarray, known_samples: set[str], threshold: float
) -> tuple[pd.DataFrame, np.ndarray]:
    called = np.all(gt >= 0, axis=2)
    n_variants = int(gt.shape[0])
    n_called = called.sum(axis=0)
    call_rate = n_called / max(n_variants, 1)
    is_het = (gt[:, :, 0] != gt[:, :, 1]) & called
    het = is_het.sum(axis=0) / np.maximum(n_called, 1)
    reasons: list[str] = []
    keep: list[bool] = []
    for sample, rate in zip(samples, call_rate, strict=True):
        reason_parts: list[str] = []
        if str(sample) not in known_samples:
            reason_parts.append("missing_sample_metadata")
        if rate < threshold:
            reason_parts.append("sample_call_rate_below_threshold")
        reasons.append(";".join(reason_parts))
        keep.append(not reason_parts)
    return pd.DataFrame({
        "sample_id": samples.astype(str),
        "n_variants": n_variants,
        "n_called": n_called.astype(int),
        "n_missing": (n_variants - n_called).astype(int),
        "call_rate": np.round(call_rate, 6),
        "heterozygosity": np.round(het, 6),
        "pass": keep,
        "exclusion_reason": reasons,
    }), np.asarray(keep, dtype=bool)


def _site_metrics(
    gt: np.ndarray,
    chrom: np.ndarray,
    pos: np.ndarray,
    ref: np.ndarray,
    alt: np.ndarray,
    cfg: SampleQCConfig,
) -> tuple[pd.DataFrame, np.ndarray]:
    called = np.all(gt >= 0, axis=2)
    n_samples = int(gt.shape[1])
    n_called = called.sum(axis=1)
    call_rate = n_called / max(n_samples, 1)
    allele_total = (called.sum(axis=1) * 2).astype(float)
    alt_count = np.where(called[:, :, None], np.maximum(gt, 0), 0).sum(axis=(1, 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        alt_freq = np.divide(alt_count, allele_total, out=np.zeros_like(alt_count, dtype=float),
                             where=allele_total > 0)
    maf = np.minimum(alt_freq, 1.0 - alt_freq)
    reasons: list[str] = []
    keep: list[bool] = []
    for rate, value in zip(call_rate, maf, strict=True):
        reason_parts: list[str] = []
        if rate < cfg.min_site_call_rate:
            reason_parts.append("site_call_rate_below_threshold")
        if value < cfg.min_site_maf:
            reason_parts.append("site_maf_below_threshold")
        reasons.append(";".join(reason_parts))
        keep.append(not reason_parts)
    return pd.DataFrame({
        "chrom": chrom.astype(str), "position": pos.astype(int), "ref": ref.astype(str),
        "alt": alt.astype(str), "n_samples": n_samples, "n_called": n_called.astype(int),
        "call_rate": np.round(call_rate, 6), "maf": np.round(maf, 6), "pass": keep,
        "exclusion_reason": reasons,
    }), np.asarray(keep, dtype=bool)


def _sample_exclusions(table: pd.DataFrame) -> pd.DataFrame:
    failed = table.loc[~table["pass"], ["sample_id", "exclusion_reason"]].copy()
    failed.insert(0, "entity_type", "sample")
    failed.rename(columns={"sample_id": "entity"}, inplace=True)
    return failed


def _site_exclusions(table: pd.DataFrame) -> pd.DataFrame:
    failed = table.loc[~table["pass"], ["chrom", "position", "ref", "alt", "exclusion_reason"]]
    entities = (
        failed["chrom"].astype(str) + ":" + failed["position"].astype(int).astype(str)
        + ":" + failed["ref"].astype(str) + ":" + failed["alt"].astype(str)
    )
    out = pd.DataFrame({
        "entity_type": ["site"] * len(failed),
        "entity": entities.to_numpy(),
        "exclusion_reason": failed["exclusion_reason"],
    })
    return out


def _atomic_csv(table: pd.DataFrame, destination: Path) -> Path:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    table.to_csv(temporary, index=False)
    temporary.replace(destination)
    return destination


def _write_filtered_vcf(
    destination: Path,
    chrom: np.ndarray,
    pos: np.ndarray,
    ref: np.ndarray,
    alt: np.ndarray,
    samples: np.ndarray,
    gt: np.ndarray,
) -> Path:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    write_minimal_vcf(temporary, chrom, pos, ref, alt, samples, gt)
    temporary.replace(destination)
    return destination
