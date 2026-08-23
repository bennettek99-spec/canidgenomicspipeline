"""Three-way coyote/wolf/dog ML admixture on bridge-panel queries.

Fits per-sample and cohort-level mixtures with locus bootstrap CIs. Designed for
eastern-coyote style diagnostics with western negative controls.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import cast

import numpy as np
from pydantic import Field

from canidae.analysis.admixture_ml import (
    BOOTSTRAP_N,
    GRID_STEP,
    MIN_DOG_FRACTION_FOR_BREED,
    MIN_PANEL_SEPARATION,
    SHRINKAGE_PSEUDOCOUNT,
    best_mixture,
    bootstrap_ci,
    calls_mixture_loglik,
    cohort_grid,
    panel_frequencies,
    two_way_fit,
    wgs_allele_count_matrix,
)
from canidae.analysis.breed_panel import (
    INFERRED_VILLAGE_CODES,
    breed_of,
    classify_group,
)
from canidae.analysis.bridge_panel import (
    BRIDGE_WGS_COLUMNS_BUG,
    BRIDGE_WGS_IDS,
    allele_count_vector,
    load_bridge_genotypes,
)
from canidae.core.errors import StageInputError
from canidae.core.logging import get_logger
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.hybrid.io_utils import (
    load_wgs_panel_json,
    require_file,
    resolve_path,
    write_csv_dicts,
    write_json,
)

_log = get_logger("stage.multiway_admixture")

_DEFAULT_LIMITATIONS = [
    "Bridge loci only; intervals are wide by WGS standards.",
    "Wolf and dog panels share derived alleles, so the wolf/dog split at this marker "
    "count is a coarse decomposition, not a definitive partition.",
    "Query and reference genotypes may come from different studies and platforms.",
]


class MultiwayAdmixtureConfig(StageConfig):
    bridge_vcf: Path = Field(..., description="Bridge VCF with query genotypes.")
    wgs_genotypes: Path = Field(
        ..., description="JSON with WGS panel samples + loci genotype lists."
    )
    western_ids: list[str] = Field(
        default_factory=list,
        description="Sample IDs treated as western/negative-control coyotes.",
    )
    western_id_prefix: str = ""
    western_id_range: list[int] = Field(
        default_factory=list,
        description="If length 2, expand prefix+zero-padded IDs in [start, end].",
    )
    seed: int = 20260816
    bootstrap_n: int = Field(default=BOOTSTRAP_N, ge=10)
    min_group: int = Field(default=3, ge=2)
    top_k: int = Field(default=5, ge=1)
    min_dog_fraction_for_breed: float = Field(default=MIN_DOG_FRACTION_FOR_BREED, ge=0.0)
    min_panel_separation: float = Field(default=MIN_PANEL_SEPARATION, gt=0.0)
    exclude_bridge_wgs_bug_ids: bool = True


@STAGES.register("multiway_admixture")
class MultiwayAdmixtureStage(Stage):
    name = "multiway_admixture"
    config_model = MultiwayAdmixtureConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return []

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "multiway_admixture")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: MultiwayAdmixtureConfig = self.config  # type: ignore[assignment]
        root = ctx.config.paths.root
        bridge = require_file(resolve_path(cfg.bridge_vcf, root), "bridge VCF")
        wgs_path = require_file(resolve_path(cfg.wgs_genotypes, root), "WGS genotypes JSON")
        rng = np.random.default_rng(cfg.seed)
        western_ids = _western_ids(cfg)

        wgs_samples, records, keys = load_wgs_panel_json(wgs_path)
        groups: dict[str, list[int]] = defaultdict(list)
        for index, sample in enumerate(wgs_samples):
            groups[breed_of(sample)].append(index)
        categories = {group: classify_group(group) for group in groups}

        coyote_indices = [
            index
            for group, indices in groups.items()
            if categories[group] == "wild" and "coyote" in group.lower()
            for index in indices
        ]
        gray_wolf_indices = [
            index for group in ("Wolf", "AlaskanWolf") for index in groups.get(group, [])
        ]
        eastern_wolf_indices = [
            index for group in ("AlgonquinWolf", "QuebecWolf") for index in groups.get(group, [])
        ]
        pooled_dogs = [
            index
            for group, indices in groups.items()
            if categories[group] != "wild"
            for index in indices
        ]
        if not (coyote_indices and gray_wolf_indices and pooled_dogs):
            raise StageInputError(
                f"reference panels incomplete: coyote={len(coyote_indices)}, "
                f"gray_wolf={len(gray_wolf_indices)}, dogs={len(pooled_dogs)}"
            )

        matrix = wgs_allele_count_matrix(keys, records)
        raw: dict[str, np.ndarray] = {}
        for name, indices in (
            ("coyote", coyote_indices),
            ("gray_wolf", gray_wolf_indices),
            ("eastern_wolf", eastern_wolf_indices or []),
            ("dog", pooled_dogs),
        ):
            panel_counts = matrix[indices]
            counts = np.nansum(panel_counts, axis=0)
            called = np.sum(~np.isnan(panel_counts), axis=0)
            freqs = (counts + 1.0) / (2.0 * called + 2.0)
            freqs[called == 0] = np.nan
            raw[name] = freqs
        pooled_prior = np.nanmean(
            np.vstack([raw["coyote"], raw["gray_wolf"], raw["eastern_wolf"], raw["dog"]]),
            axis=0,
        )
        p_coyote = panel_frequencies(matrix, coyote_indices, pooled_prior)
        p_gray = panel_frequencies(matrix, gray_wolf_indices, pooled_prior)
        p_dog = panel_frequencies(matrix, pooled_dogs, pooled_prior)
        informative = (
            (
                (np.abs(raw["gray_wolf"] - raw["coyote"]) >= cfg.min_panel_separation)
                | (np.abs(raw["dog"] - raw["coyote"]) >= cfg.min_panel_separation)
                | (np.abs(raw["eastern_wolf"] - raw["coyote"]) >= cfg.min_panel_separation)
            )
            & np.isfinite(p_coyote)
            & np.isfinite(p_gray)
            & np.isfinite(p_dog)
        )
        keep = np.where(informative)[0]
        if len(keep) < 10:
            raise StageInputError(f"only {len(keep)} informative loci after separation filter")
        kept_keys = [keys[i] for i in keep]
        p_coyote_k, p_gray_k, p_dog_k = p_coyote[keep], p_gray[keep], p_dog[keep]

        bridge_samples, bridge_calls = load_bridge_genotypes(bridge)
        radseq_ids = [
            s
            for s in bridge_samples
            if s not in wgs_samples
            and not (
                cfg.exclude_bridge_wgs_bug_ids and BRIDGE_WGS_COLUMNS_BUG and s in BRIDGE_WGS_IDS
            )
        ]
        if not radseq_ids:
            raise StageInputError("no query samples remain after excluding WGS/bridge IDs")

        candidates = {
            group: indices
            for group, indices in groups.items()
            if categories[group] == "breed" and len(indices) >= cfg.min_group
        }
        village_indices = [
            index
            for group, indices in groups.items()
            if categories[group] == "village"
            for index in indices
        ]
        if len(village_indices) >= cfg.min_group:
            candidates["VillageDog(all regions)"] = village_indices
        for code, region in INFERRED_VILLAGE_CODES.items():
            if len(groups.get(code, [])) >= cfg.min_group:
                candidates[f"VillageDog({region}, inferred)"] = groups[code]
        breed_panels = {
            breed: panel_frequencies(matrix, indices, pooled_prior)[keep]
            for breed, indices in candidates.items()
        }

        rows: list[dict[str, object]] = []
        breed_rows: list[dict[str, object]] = []
        for sample in radseq_ids:
            region = "western" if sample in western_ids else "eastern"
            counts = allele_count_vector(kept_keys, bridge_calls, sample)
            f_wolf, f_dog, ll = best_mixture(counts, p_coyote_k, p_gray_k, p_dog_k)
            ci = bootstrap_ci(
                counts,
                p_coyote_k,
                p_gray_k,
                p_dog_k,
                rng,
                bootstrap_n=cfg.bootstrap_n,
            )
            rows.append(
                {
                    "sample_id": sample,
                    "region": region,
                    "loci_called": int(np.sum(~np.isnan(counts))),
                    "f_wolf": round(f_wolf, 3),
                    "f_dog": round(f_dog, 3),
                    "f_wolf_ci_lo": ci["wolf"][0],
                    "f_wolf_ci_hi": ci["wolf"][1],
                    "f_dog_ci_lo": ci["dog"][0],
                    "f_dog_ci_hi": ci["dog"][1],
                    "loglik": round(ll, 1),
                }
            )
            if f_dog >= cfg.min_dog_fraction_for_breed and breed_panels:
                entries = [
                    (
                        calls_mixture_loglik(counts, p_coyote_k, p_gray_k, panel, f_wolf, f_dog),
                        breed,
                    )
                    for breed, panel in breed_panels.items()
                ]
                any_dog = calls_mixture_loglik(counts, p_coyote_k, p_gray_k, p_dog_k, f_wolf, f_dog)
                entries.sort(reverse=True)
                for rank, (score, breed) in enumerate(entries[: cfg.top_k], start=1):
                    breed_rows.append(
                        {
                            "sample_id": sample,
                            "rank": rank,
                            "breed": breed,
                            "log_likelihood": round(score, 2),
                            "gap_to_any_dog": round(score - any_dog, 2),
                        }
                    )

        out_csv = ctx.datastore.path_for(self.name, "multiway_admixture.csv")
        write_csv_dicts(out_csv, rows)
        if breed_rows:
            write_csv_dicts(ctx.datastore.path_for(self.name, "breed_scores.csv"), breed_rows)

        eastern = [r for r in rows if r["region"] == "eastern"]
        western = [r for r in rows if r["region"] == "western"]

        def _cohort(sample_ids: list[str]) -> dict[str, object]:
            if not sample_ids:
                return {"n_samples": 0}
            stacked = np.vstack(
                [allele_count_vector(kept_keys, bridge_calls, s) for s in sample_ids]
            )
            alt_counts = np.nansum(stacked, axis=0)
            chromosome_counts = 2 * np.sum(~np.isnan(stacked), axis=0)
            f_wolf, f_dog, ll = cohort_grid(
                alt_counts, chromosome_counts, p_coyote_k, p_gray_k, p_dog_k
            )
            draws = []
            use = chromosome_counts > 0
            n_use = int(use.sum())
            for _ in range(cfg.bootstrap_n):
                pick = rng.integers(0, n_use, n_use)
                fw_b, fd_b, _ = cohort_grid(
                    alt_counts[use][pick],
                    chromosome_counts[use][pick],
                    p_coyote_k[use][pick],
                    p_gray_k[use][pick],
                    p_dog_k[use][pick],
                )
                draws.append((fw_b, fd_b))
            draws_arr = np.array(draws)
            return {
                "n_samples": len(sample_ids),
                "f_wolf": round(f_wolf, 3),
                "f_dog": round(f_dog, 3),
                "f_wolf_ci": [round(v, 3) for v in np.percentile(draws_arr[:, 0], [2.5, 97.5])],
                "f_dog_ci": [round(v, 3) for v in np.percentile(draws_arr[:, 1], [2.5, 97.5])],
                "loglik": round(ll, 1),
            }

        cohort_eastern = _cohort([str(r["sample_id"]) for r in eastern])
        cohort_western = _cohort([str(r["sample_id"]) for r in western])

        discrimination: dict[str, object] = {}
        if eastern:
            stacked = np.vstack(
                [allele_count_vector(kept_keys, bridge_calls, str(r["sample_id"])) for r in eastern]
            )
            alt_e = np.nansum(stacked, axis=0)
            chrom_e = 2 * np.sum(~np.isnan(stacked), axis=0)
            f_dog_2w, ll_dog_2w = two_way_fit(alt_e, chrom_e, p_coyote_k, p_dog_k)
            f_wolf_2w, ll_wolf_2w = two_way_fit(alt_e, chrom_e, p_coyote_k, p_gray_k)
            discrimination = {
                "coyote_to_dog": {"f": round(f_dog_2w, 3), "loglik": round(ll_dog_2w, 1)},
                "coyote_to_wolf": {
                    "f": round(f_wolf_2w, 3),
                    "loglik": round(ll_wolf_2w, 1),
                },
                "loglik_difference": round(ll_dog_2w - ll_wolf_2w, 2),
            }

        def _stats(values: list[float]) -> dict[str, float]:
            if not values:
                return {"mean": float("nan"), "min": float("nan"), "max": float("nan")}
            return {
                "mean": round(float(np.mean(values)), 3),
                "min": round(float(np.min(values)), 3),
                "max": round(float(np.max(values)), 3),
            }

        group_summary = {
            "eastern_wolf_fraction": _stats([cast(float, r["f_wolf"]) for r in eastern]),
            "western_wolf_fraction": _stats([cast(float, r["f_wolf"]) for r in western]),
            "eastern_dog_fraction": _stats([cast(float, r["f_dog"]) for r in eastern]),
            "western_dog_fraction": _stats([cast(float, r["f_dog"]) for r in western]),
        }
        manifest = {
            "queries": {
                "n_queries": len(radseq_ids),
                "n_eastern": len(eastern),
                "n_western": len(western),
            },
            "references": {
                "coyote_n": len(coyote_indices),
                "gray_wolf_n": len(gray_wolf_indices),
                "eastern_wolf_n": len(eastern_wolf_indices),
                "dogs_n": len(pooled_dogs),
                "loci_total": len(keys),
                "loci_informative": len(keep),
            },
            "group_summary": group_summary,
            "cohort_estimates": {
                "eastern": cohort_eastern,
                "western_control": cohort_western,
                "wolf_vs_dog_discrimination": discrimination,
            },
            "method": {
                "model": "3-way maximum-likelihood admixture (coyote/gray wolf/dog)",
                "grid_step": GRID_STEP,
                "bootstrap_n": cfg.bootstrap_n,
                "shrinkage_pseudocount": SHRINKAGE_PSEUDOCOUNT,
                "min_panel_separation": cfg.min_panel_separation,
            },
            "limitations": list(_DEFAULT_LIMITATIONS),
        }
        manifest_path = ctx.datastore.path_for(self.name, "multiway_admixture.json")
        write_json(manifest_path, manifest)

        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT,
            "multiway_admixture",
            out_csv,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "multiway_admixture",
                "diagnostic": True,
                "n_queries": len(rows),
                "loci_informative": len(keep),
                "group_summary": group_summary,
                "cohort_estimates": manifest["cohort_estimates"],
                "limitations": list(_DEFAULT_LIMITATIONS),
                "manifest": str(manifest_path),
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={
                "n_queries": len(rows),
                "loci_informative": len(keep),
                "eastern_mean_f_dog": group_summary["eastern_dog_fraction"]["mean"],
                "western_mean_f_dog": group_summary["western_dog_fraction"]["mean"],
            },
        )


def _western_ids(cfg: MultiwayAdmixtureConfig) -> frozenset[str]:
    ids = set(cfg.western_ids)
    if cfg.western_id_prefix and len(cfg.western_id_range) == 2:
        start, end = cfg.western_id_range
        width = max(3, len(str(end)))
        for number in range(start, end + 1):
            ids.add(f"{cfg.western_id_prefix}{number:0{width}d}")
    return frozenset(ids)
