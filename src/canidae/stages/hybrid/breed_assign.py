"""Breed-panel scoring of dog ancestry on hybrid individuals.

Scores query dog components against breed allele-frequency panels built from a
full WGS genotype matrix at bridge loci, with leave-one-out calibration and a
pooled any-dog baseline for "no single breed" detection.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import cast

from pydantic import Field

from canidae.analysis.breed_panel import (
    INFERRED_VILLAGE_CODES,
    breed_of,
    classify_group,
)
from canidae.analysis.bridge_panel import bridge_sites
from canidae.analysis.reference_mixture import (
    SINGLE_BREED_GAP,
    allele_frequencies,
    calls_log_likelihood,
    leave_one_out,
)
from canidae.core.errors import StageInputError
from canidae.core.logging import get_logger
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.hybrid.io_utils import (
    load_calls_csv,
    load_json,
    load_wgs_panel_json,
    parse_locus_key,
    require_file,
    resolve_path,
    write_csv_dicts,
    write_json,
)

_log = get_logger("stage.breed_assign")

_DEFAULT_LIMITATIONS = [
    "Bridge-locus breed resolution is far below a genotype-array or WGS standard.",
    "The model scores which breed panel best explains dog-derived alleles; it cannot "
    "phase or directly observe the parent.",
    "A mixed-breed or unpanelled dog parent should appear as a small single-breed gap "
    "versus the pooled any-dog panel.",
]


class BreedAssignConfig(StageConfig):
    bridge_vcf: Path = Field(..., description="Bridge VCF defining shared loci.")
    wgs_genotypes: Path = Field(
        ..., description="JSON with samples + loci genotype lists (all-panel fetch)."
    )
    calls_dir: Path = Field(..., description="Directory with <sample>_calls.csv tables.")
    dog_fractions: Path = Field(
        ...,
        description=(
            "JSON or CSV providing per-query dog_fraction. JSON may be a "
            "reference_mixture manifest (results[]) or a plain sample->fraction map."
        ),
    )
    query_samples: list[str] = Field(default_factory=list)
    min_depth: int = Field(default=8, ge=4)
    min_group: int = Field(default=3, ge=2)
    top_k: int = Field(default=5, ge=1)
    single_breed_gap: float = Field(default=SINGLE_BREED_GAP, gt=0.0)
    # Optional locus restriction (e.g. previously retained NYC reference keys).
    locus_filter_json: Path | None = None


@STAGES.register("breed_assign")
class BreedAssignStage(Stage):
    name = "breed_assign"
    config_model = BreedAssignConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "reference_mixture", optional=True),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "breed_assign")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: BreedAssignConfig = self.config  # type: ignore[assignment]
        root = ctx.config.paths.root
        bridge = require_file(resolve_path(cfg.bridge_vcf, root), "bridge VCF")
        wgs_path = require_file(resolve_path(cfg.wgs_genotypes, root), "WGS genotypes JSON")
        calls_dir = require_file(resolve_path(cfg.calls_dir, root), "calls directory")
        frac_path = require_file(resolve_path(cfg.dog_fractions, root), "dog fractions")

        sites = bridge_sites(bridge)
        samples, records, keys = load_wgs_panel_json(wgs_path)
        if cfg.locus_filter_json is not None:
            filt = require_file(resolve_path(cfg.locus_filter_json, root), "locus filter")
            saved = load_json(filt)
            keys = sorted(
                (parse_locus_key(k) for k in saved),
                key=lambda item: (int(item[0].removeprefix("chr")), item[1]),
            )
            records = {key: records[key] for key in keys if key in records}
            keys = [k for k in keys if k in records]
        if len(keys) < 20:
            raise StageInputError(f"only {len(keys)} loci available for breed assignment")

        dog_fractions = _load_dog_fractions(frac_path, ctx)
        query_ids = cfg.query_samples or sorted(dog_fractions)
        if not query_ids:
            raise StageInputError("no query samples for breed assignment")

        groups: dict[str, list[int]] = defaultdict(list)
        for index, sample in enumerate(samples):
            groups[breed_of(sample)].append(index)
        categories = {group: classify_group(group) for group in groups}
        candidates = _candidate_groups(groups, categories, cfg.min_group)
        wild_groups = {
            group: indices for group, indices in groups.items() if categories[group] == "wild"
        }
        pooled_dogs = [
            index
            for group, indices in groups.items()
            if categories[group] != "wild"
            for index in indices
        ]
        coyote_indices = next(
            (indices for name, indices in wild_groups.items() if "coyote" in name.lower()),
            [],
        )
        if not coyote_indices:
            raise StageInputError("no coyote group available for the mixture model")

        coyote_panel = allele_frequencies(keys, records, coyote_indices)
        pooled_panel = allele_frequencies(keys, records, pooled_dogs)
        panels = {
            group: allele_frequencies(keys, records, indices)
            for group, indices in candidates.items()
        }

        loo_rows, loo_summary = leave_one_out(keys, records, samples, candidates, cfg.min_group)
        loo_path = ctx.datastore.path_for(self.name, "leave_one_out.csv")
        write_csv_dicts(loo_path, loo_rows)

        score_rows: list[dict[str, object]] = []
        per_sample: dict[str, dict[str, object]] = {}
        for sample_id in query_ids:
            dog_fraction = dog_fractions.get(sample_id)
            if dog_fraction is None:
                _log.warning("skipping %s: no dog_fraction", sample_id)
                continue
            call_path = calls_dir / f"{sample_id}_calls.csv"
            require_file(call_path, f"calls for {sample_id}")
            selected = {key: sites[key] for key in keys if key in sites}
            calls = load_calls_csv(call_path, selected, cfg.min_depth)
            entries: list[dict[str, object]] = []
            for breed, panel in panels.items():
                ll, used = calls_log_likelihood(
                    calls, keys, panel, coyote_panel, float(dog_fraction)
                )
                entries.append({"breed": breed, "log_likelihood": round(ll, 2), "loci_used": used})
            any_ll, any_used = calls_log_likelihood(
                calls, keys, pooled_panel, coyote_panel, float(dog_fraction)
            )
            entries.sort(key=lambda entry: cast(float, entry["log_likelihood"]), reverse=True)
            for rank, entry in enumerate(entries[: cfg.top_k], start=1):
                score_rows.append(
                    {
                        "sample_id": sample_id,
                        "rank": rank,
                        "breed": entry["breed"],
                        "log_likelihood": entry["log_likelihood"],
                        "gap_to_next": (
                            round(
                                cast(float, entries[rank - 1]["log_likelihood"])
                                - cast(float, entries[rank]["log_likelihood"]),
                                2,
                            )
                            if rank < len(entries)
                            else None
                        ),
                        "loci_used": entry["loci_used"],
                    }
                )
            best = entries[0]
            gap = cast(float, best["log_likelihood"]) - any_ll
            per_sample[sample_id] = {
                "dog_fraction": dog_fraction,
                "best_breed": best["breed"],
                "best_log_likelihood": best["log_likelihood"],
                "any_dog_log_likelihood": round(any_ll, 2),
                "single_breed_gap": round(gap, 2),
                "single_breed_supported": bool(gap > cfg.single_breed_gap),
                "loci_used": any_used,
            }
            _log.info("%s best=%s gap=%.2f", sample_id, best["breed"], gap)

        scores_path = ctx.datastore.path_for(self.name, "breed_scores.csv")
        write_csv_dicts(scores_path, score_rows)
        summary_path = ctx.datastore.path_for(self.name, "breed_assignment.json")
        write_json(
            summary_path,
            {
                "calibration": loo_summary,
                "results": per_sample,
                "limitations": list(_DEFAULT_LIMITATIONS),
                "sources": {
                    "n_wgs_samples": len(samples),
                    "n_candidate_groups": len(candidates),
                    "loci": len(keys),
                },
            },
        )

        # Primary tabular artifact for the report: one row per query (top breed).
        summary_rows = [
            {
                "sample_id": sample_id,
                "dog_fraction": info["dog_fraction"],
                "best_breed": info["best_breed"],
                "single_breed_gap": info["single_breed_gap"],
                "single_breed_supported": info["single_breed_supported"],
                "any_dog_log_likelihood": info["any_dog_log_likelihood"],
                "loci_used": info["loci_used"],
            }
            for sample_id, info in per_sample.items()
        ]
        out_csv = ctx.datastore.path_for(self.name, "breed_assign.csv")
        write_csv_dicts(out_csv, summary_rows)

        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT,
            "breed_assign",
            out_csv,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "breed_assign",
                "diagnostic": True,
                "calibration": loo_summary,
                "n_queries": len(summary_rows),
                "n_candidates": len(candidates),
                "limitations": list(_DEFAULT_LIMITATIONS),
                "scores_csv": str(scores_path),
                "manifest": str(summary_path),
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={
                "n_queries": len(summary_rows),
                "loo_top1_accuracy": loo_summary.get("top1_accuracy"),
            },
        )


def _candidate_groups(
    groups: dict[str, list[int]],
    categories: dict[str, str],
    min_group: int,
) -> dict[str, list[int]]:
    candidates: dict[str, list[int]] = {
        group: indices
        for group, indices in groups.items()
        if categories[group] == "breed" and len(indices) >= min_group
    }
    village_indices = [
        index
        for group, indices in groups.items()
        if categories[group] == "village"
        for index in indices
    ]
    if len(village_indices) >= min_group:
        candidates["VillageDog(all regions)"] = village_indices
    for code, region in INFERRED_VILLAGE_CODES.items():
        if len(groups.get(code, [])) >= min_group:
            candidates[f"VillageDog({region}, inferred)"] = groups[code]
    return candidates


def _load_dog_fractions(path: Path, ctx: RunContext) -> dict[str, float]:
    if path.suffix.lower() == ".csv":
        import pandas as pd

        df = pd.read_csv(path)
        if "sample_id" not in df.columns or "dog_fraction" not in df.columns:
            raise StageInputError("dog_fractions CSV must have sample_id and dog_fraction columns")
        return {
            str(row.sample_id): float(row.dog_fraction)
            for row in df.itertuples(index=False)
            if row.dog_fraction == row.dog_fraction  # not NaN
        }
    payload = load_json(path)
    if "results" in payload and isinstance(payload["results"], list):
        return {
            str(row["sample_id"]): float(row["dog_fraction"])
            for row in payload["results"]
            if row.get("dog_fraction") is not None
        }
    if all(isinstance(v, (int, float)) for v in payload.values()):
        return {str(k): float(v) for k, v in payload.items()}
    # Prefer upstream stage artifact when present.
    if ctx.datastore.has(ArtifactKind.ANALYSIS_RESULT, "reference_mixture"):
        art = ctx.datastore.get(ArtifactKind.ANALYSIS_RESULT, "reference_mixture")
        import pandas as pd

        df = pd.read_csv(art.path)
        return {
            str(row.sample_id): float(row.dog_fraction)
            for row in df.itertuples(index=False)
            if row.dog_fraction == row.dog_fraction
        }
    raise StageInputError(f"unrecognized dog_fractions format: {path}")
