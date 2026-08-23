"""Two-source (coyote vs dog) mixture stage for hybrid-canid pedigree checks.

Consumes pre-called query allele-count tables and a small WGS reference genotype
JSON (as produced by the NYC coydog validation workflow). Emits a CSV of dog
fractions plus a JSON validation verdict. This is a bridge-locus diagnostic.
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

from pydantic import Field

from canidae.analysis.bridge_panel import bridge_sites
from canidae.analysis.reference_mixture import infer_dog_fraction
from canidae.core.errors import StageInputError
from canidae.core.logging import get_logger
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.hybrid.io_utils import (
    load_calls_csv,
    load_json,
    parse_locus_key,
    require_file,
    resolve_path,
    write_csv_dicts,
    write_json,
)

_log = get_logger("stage.reference_mixture")

_DEFAULT_LIMITATIONS = [
    "Cross-platform bridge loci rather than a jointly called whole genome.",
    "Dog ancestry is relative to the configured reference dogs, not a local dog panel.",
    "Use as a reproducible pedigree/phenotype check, not a legal, veterinary, "
    "or management diagnosis.",
]


class ReferenceMixtureConfig(StageConfig):
    """Configuration for :class:`ReferenceMixtureStage`."""

    bridge_vcf: Path = Field(..., description="Bridge VCF defining shared loci.")
    reference_genotypes: Path = Field(
        ..., description="JSON map of chrom:pos -> list of WGS genotype strings."
    )
    calls_dir: Path = Field(
        ..., description="Directory with <sample_id>_calls.csv allele-count tables."
    )
    query_samples: list[str] = Field(
        default_factory=list,
        description="Query sample IDs (must match calls CSV stems).",
    )
    wgs_samples: list[str] = Field(
        default_factory=list,
        description="Order of genotypes inside each reference_genotypes value list.",
    )
    coyote_samples: list[str] = Field(
        default_factory=list,
        description="Subset of wgs_samples treated as the coyote reference.",
    )
    min_depth: int = Field(default=8, ge=4)
    min_difference: float = Field(default=0.40, gt=0.0)
    min_called: int = Field(default=20, ge=1)
    # Optional pedigree-style pass criteria (NYC coydog defaults when left empty).
    expected_f1: str = ""
    expected_parent: str = ""
    expected_offspring: list[str] = Field(default_factory=list)
    f1_target: float = 0.46
    f1_tolerance: float = 0.20
    parent_max: float = 0.20
    offspring_min: float = 0.08
    offspring_max: float = 0.45


@STAGES.register("reference_mixture")
class ReferenceMixtureStage(Stage):
    name = "reference_mixture"
    config_model = ReferenceMixtureConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []  # entry-style: paths come from config (like ingest)

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "reference_mixture")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: ReferenceMixtureConfig = self.config  # type: ignore[assignment]
        root = ctx.config.paths.root
        bridge = require_file(resolve_path(cfg.bridge_vcf, root), "bridge VCF")
        ref_path = require_file(
            resolve_path(cfg.reference_genotypes, root), "reference genotypes JSON"
        )
        calls_dir = require_file(resolve_path(cfg.calls_dir, root), "calls directory")
        if not cfg.query_samples:
            raise StageInputError("stages.reference_mixture.query_samples is required")
        if not cfg.wgs_samples:
            raise StageInputError("stages.reference_mixture.wgs_samples is required")
        coyote = frozenset(cfg.coyote_samples) or frozenset(
            s for s in cfg.wgs_samples if "coyote" in s.lower()
        )
        if not coyote:
            raise StageInputError(
                "stages.reference_mixture.coyote_samples is empty and none could be inferred"
            )

        sites = bridge_sites(bridge)
        saved = load_json(ref_path)
        reference: dict[tuple[str, int], list[str]] = {
            parse_locus_key(key): [str(v) for v in values] for key, values in saved.items()
        }
        if len(reference) < 20:
            raise StageInputError(
                f"reference genotypes retain only {len(reference)} loci; need more shared sites"
            )

        rows: list[dict[str, object]] = []
        for sample_id in cfg.query_samples:
            call_path = calls_dir / f"{sample_id}_calls.csv"
            require_file(call_path, f"calls for {sample_id}")
            selected = {key: sites[key] for key in reference if key in sites}
            if len(selected) < 20:
                raise StageInputError(
                    f"{sample_id}: fewer than 20 bridge sites overlap the reference panel"
                )
            calls = load_calls_csv(call_path, selected, cfg.min_depth)
            dog_fraction, diagnostic_calls, rmse = infer_dog_fraction(
                calls,
                {k: reference[k] for k in selected},
                wgs_samples=cfg.wgs_samples,
                coyote_samples=coyote,
                min_difference=cfg.min_difference,
                min_called=cfg.min_called,
            )
            mean_depth = float(sum(depth for _, depth, _, _ in calls.values()) / max(len(calls), 1))
            rows.append(
                {
                    "sample_id": sample_id,
                    "dog_fraction": dog_fraction,
                    "diagnostic_sites_called": diagnostic_calls,
                    "mean_depth": round(mean_depth, 2),
                    "residual_rmse": rmse,
                }
            )
            _log.info(
                "%s dog_fraction=%s diagnostic_sites=%s",
                sample_id,
                dog_fraction,
                diagnostic_calls,
            )

        verdict = _pedigree_verdict(rows, cfg)
        out_csv = ctx.datastore.path_for(self.name, "reference_mixture.csv")
        write_csv_dicts(out_csv, rows)
        out_json = ctx.datastore.path_for(self.name, "reference_mixture.json")
        payload = {
            "results": rows,
            "validation": verdict,
            "limitations": list(_DEFAULT_LIMITATIONS),
            "method": {
                "model": "two-source least-squares dog fraction (coyote vs dog means)",
                "min_depth": cfg.min_depth,
                "min_difference": cfg.min_difference,
                "min_called": cfg.min_called,
                "n_reference_loci": len(reference),
                "wgs_samples": list(cfg.wgs_samples),
                "coyote_samples": sorted(coyote),
            },
        }
        write_json(out_json, payload)

        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT,
            "reference_mixture",
            out_csv,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "reference_mixture",
                "diagnostic": True,
                "validation_passed": verdict.get("passed"),
                "n_queries": len(rows),
                "n_reference_loci": len(reference),
                "limitations": list(_DEFAULT_LIMITATIONS),
                "manifest": str(out_json),
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={
                "n_queries": len(rows),
                "validation_passed": verdict.get("passed"),
            },
        )


def _pedigree_verdict(
    rows: list[dict[str, object]], cfg: ReferenceMixtureConfig
) -> dict[str, object]:
    """Optional pedigree-style pass/fail when F1/parent/offspring IDs are configured."""
    if not (cfg.expected_f1 and cfg.expected_parent and cfg.expected_offspring):
        return {
            "passed": None,
            "criterion": "no pedigree criteria configured",
            "interpretation": (
                "Dog fractions were estimated; configure expected_f1 / expected_parent / "
                "expected_offspring to enable an automatic pedigree-style check."
            ),
        }
    by_id = {str(row["sample_id"]): row for row in rows}
    missing = [
        name
        for name in [cfg.expected_f1, cfg.expected_parent, *cfg.expected_offspring]
        if name not in by_id
    ]
    if missing:
        return {
            "passed": False,
            "criterion": "configured pedigree sample IDs must be present in results",
            "interpretation": f"missing samples: {missing}",
        }
    f1 = by_id[cfg.expected_f1]["dog_fraction"]
    parent = by_id[cfg.expected_parent]["dog_fraction"]
    offspring = [by_id[name]["dog_fraction"] for name in cfg.expected_offspring]
    available = all(value is not None for value in [f1, parent, *offspring])
    passed = bool(
        available
        and abs(float(cast(float, f1)) - cfg.f1_target) <= cfg.f1_tolerance
        and float(cast(float, parent)) <= cfg.parent_max
        and all(cfg.offspring_min <= float(cast(float, v)) <= cfg.offspring_max for v in offspring)
        and all(float(cast(float, v)) < float(cast(float, f1)) for v in offspring)
    )
    return {
        "passed": passed,
        "criterion": (
            f"{cfg.expected_f1} near {cfg.f1_target}±{cfg.f1_tolerance}, "
            f"{cfg.expected_parent} <= {cfg.parent_max}, offspring in "
            f"[{cfg.offspring_min}, {cfg.offspring_max}] and below F1"
        ),
        "interpretation": (
            "Pass supports the expected hybrid pedigree pattern on this bridge panel. "
            "It is not a replacement for full local-ancestry analysis."
        ),
    }
