"""Analysis-readiness gate for population-genomic interpretation.

The stage turns raw filtered genotypes into an explicit analysis panel: non-autosomal
contigs can be removed, tightly linked markers are pruned, population sizes and near-
duplicate samples are audited, and small-panel ascertainment is carried into provenance.
"""

from __future__ import annotations

import json

import allel
import numpy as np
from pydantic import Field

from canidae.core.errors import StageInputError
from canidae.core.model import Artifact, ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import (
    Genotypes,
    load_genotypes,
    load_sample_labels,
    population_indices,
    save_genotypes,
)


class AnalysisReadinessConfig(StageConfig):
    strict: bool = False
    autosomes_only: bool = True
    min_samples_per_population: int = Field(default=2, ge=1)
    ld_prune: bool = True
    ld_window_bp: int = Field(default=500_000, ge=1)
    ld_r2_threshold: float = Field(default=0.8, gt=0.0, le=1.0)
    duplicate_concordance: float = Field(default=0.98, gt=0.5, le=1.0)


@STAGES.register("analysis_readiness")
class AnalysisReadinessStage(Stage):
    name = "analysis_readiness"
    config_model = AnalysisReadinessConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes"),
            ArtifactSpec(ArtifactKind.QC_TABLE, "analysis_readiness"),
        ]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: AnalysisReadinessConfig = self.config  # type: ignore[assignment]
        source = ctx.datastore.get(ArtifactKind.GENOTYPES, "genotypes")
        geno = load_genotypes(source.path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path
        )
        groups = population_indices(geno, labels)
        issues: list[dict[str, str]] = []

        singleton = {
            pop: len(idx)
            for pop, idx in groups.items()
            if len(idx) < cfg.min_samples_per_population
        }
        if singleton:
            issues.append(
                {
                    "severity": "error" if cfg.strict else "warning",
                    "code": "small_population",
                    "detail": ", ".join(f"{pop}={n}" for pop, n in singleton.items()),
                }
            )

        keep: np.ndarray = np.ones(geno.n_variants, dtype=bool)
        non_autosomal = np.array([not _is_autosome(value) for value in geno.chrom])
        if cfg.autosomes_only and np.any(non_autosomal):
            keep &= ~non_autosomal
            issues.append(
                {
                    "severity": "info",
                    "code": "non_autosomal_removed",
                    "detail": str(int(non_autosomal.sum())),
                }
            )

        candidate = np.flatnonzero(keep)
        ld_removed = 0
        if cfg.ld_prune and candidate.size and geno.n_samples >= 4:
            retained = _ld_prune(
                geno,
                candidate,
                window_bp=cfg.ld_window_bp,
                r2_threshold=cfg.ld_r2_threshold,
            )
            keep[candidate] = False
            keep[retained] = True
            ld_removed = int(candidate.size - retained.size)
        elif cfg.ld_prune and geno.n_samples < 4:
            issues.append(
                {
                    "severity": "warning",
                    "code": "ld_pruning_underpowered",
                    "detail": f"n_samples={geno.n_samples}",
                }
            )

        duplicates = _near_duplicates(geno, cfg.duplicate_concordance)
        if duplicates:
            issues.append(
                {
                    "severity": "error" if cfg.strict else "warning",
                    "code": "near_duplicate_samples",
                    "detail": "; ".join(f"{a}/{b}={value:.4f}" for a, b, value in duplicates),
                }
            )
        if source.metadata.get("panel_relative"):
            issues.append(
                {
                    "severity": "warning",
                    "code": "ascertainment_panel_relative",
                    "detail": str(source.metadata.get("panel_name") or "selected SNP panel"),
                }
            )
        hard_call_filters = source.metadata.get("hard_call_filters") or {}
        if source.metadata.get("panel_relative") and not any(hard_call_filters.values()):
            issues.append(
                {
                    "severity": "warning",
                    "code": "hard_call_quality_unavailable",
                    "detail": (
                        "No GQ/DP thresholds were applied; interpret low-coverage calls cautiously."
                    ),
                }
            )
        if cfg.strict and any(issue["severity"] == "error" for issue in issues):
            detail = "; ".join(issue["detail"] for issue in issues if issue["severity"] == "error")
            raise StageInputError(f"analysis readiness failed in strict mode: {detail}")
        if not np.any(keep):
            raise StageInputError("analysis readiness removed every variant")

        ready = Genotypes(
            calls=allel.GenotypeArray(np.asarray(geno.calls)[keep]),
            pos=np.asarray(geno.pos)[keep],
            chrom=np.asarray(geno.chrom)[keep],
            samples=np.asarray(geno.samples),
        )
        stage_dir = ctx.datastore.stage_dir(self.name)
        backend = "npy_mmap" if source.path.is_dir() else "npz"
        output_name = (
            "analysis_genotypes.store" if backend == "npy_mmap" else "analysis_genotypes.npz"
        )
        output = save_genotypes(stage_dir / output_name, ready, backend=backend)
        audit_path = stage_dir / "analysis_readiness.json"
        status = "warning" if any(i["severity"] == "warning" for i in issues) else "ready"
        audit = {
            "schema_version": 1,
            "status": status,
            "n_samples": ready.n_samples,
            "n_variants_input": geno.n_variants,
            "n_variants_ready": ready.n_variants,
            "ld_removed": ld_removed,
            "population_sizes": {pop: len(idx) for pop, idx in groups.items()},
            "issues": issues,
        }
        audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
        metadata = {
            **source.metadata,
            "analysis_readiness": status,
            "ld_pruned": cfg.ld_prune and geno.n_samples >= 4,
            "autosomes_only": cfg.autosomes_only,
            "n_variants_before_readiness": geno.n_variants,
        }
        genotype_artifact = Artifact(
            ArtifactKind.GENOTYPES,
            "analysis_genotypes",
            output,
            FileFormat.OTHER if backend == "npy_mmap" else FileFormat.NPZ,
            produced_by=self.name,
            metadata=metadata,
        )
        audit_artifact = Artifact(
            ArtifactKind.QC_TABLE,
            "analysis_readiness",
            audit_path,
            FileFormat.JSON,
            produced_by=self.name,
            metadata={"status": status, "issues": issues},
        )
        return StageResult(
            [genotype_artifact, audit_artifact],
            metrics={
                "status": status,
                "n_variants_input": geno.n_variants,
                "n_variants_ready": ready.n_variants,
                "ld_removed": ld_removed,
                "near_duplicate_pairs": len(duplicates),
            },
        )


def _is_autosome(value: object) -> bool:
    contig = str(value).lower().removeprefix("chr")
    return contig.isdigit() and 1 <= int(contig) <= 38


def _ld_prune(
    geno: Genotypes,
    candidates: np.ndarray,
    *,
    window_bp: int,
    r2_threshold: float,
) -> np.ndarray:
    n_alt = np.asarray(geno.calls.to_n_alt(fill=-1), dtype=np.float32)
    retained: list[int] = []
    for index in candidates:
        chrom = str(geno.chrom[index])
        position = int(geno.pos[index])
        recent = [
            prior
            for prior in reversed(retained)
            if str(geno.chrom[prior]) == chrom and position - int(geno.pos[prior]) <= window_bp
        ]
        if any(_r2(n_alt[index], n_alt[prior]) >= r2_threshold for prior in recent):
            continue
        retained.append(int(index))
    return np.asarray(retained, dtype=int)


def _r2(left: np.ndarray, right: np.ndarray) -> float:
    called = (left >= 0) & (right >= 0)
    if called.sum() < 4 or np.std(left[called]) == 0 or np.std(right[called]) == 0:
        return 0.0
    value = float(np.corrcoef(left[called], right[called])[0, 1])
    return value * value if np.isfinite(value) else 0.0


def _near_duplicates(geno: Genotypes, threshold: float) -> list[tuple[str, str, float]]:
    n_alt = np.asarray(geno.calls.to_n_alt(fill=-1))
    values: list[tuple[str, str, float]] = []
    for left in range(geno.n_samples):
        for right in range(left + 1, geno.n_samples):
            called = (n_alt[:, left] >= 0) & (n_alt[:, right] >= 0)
            if called.sum() < 100:
                continue
            concordance = float(np.mean(n_alt[called, left] == n_alt[called, right]))
            if concordance >= threshold:
                values.append((str(geno.samples[left]), str(geno.samples[right]), concordance))
    return values


__all__ = ["AnalysisReadinessConfig", "AnalysisReadinessStage"]
