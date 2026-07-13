"""D-statistics (ABBA-BABA) stage with a windowed f_d scan.

For a designated outgroup, enumerates all ((P1,P2),P3),O trios among the remaining
populations (as Dsuite's Dtrios does) and computes Patterson's D with a block-jackknife
Z-score and p-value. A significant |D| for a trio whose assumed sisters (P1,P2) are the true
sisters indicates gene flow between one of them and P3; because all assignments are scored,
the single largest |Z| can instead reflect the tree topology, so interpret trios against the
known/estimated tree. For the max-|Z| trio it also writes a windowed f_d scan to localize
introgression along the genome.
"""

from __future__ import annotations

from itertools import combinations
from typing import Literal

import numpy as np
import pandas as pd
from pydantic import Field
from scipy.stats import norm

from canidae.core.errors import StageInputError
from canidae.core.logging import get_logger
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.introgression import fstats
from canidae.stages.popgen.store import (
    Genotypes,
    load_genotypes,
    load_sample_labels,
    population_indices,
)

_log = get_logger("stage.dstats")


class DStatsConfig(StageConfig):
    # Required: biological outgroups must never be inferred from PCA distance.
    outgroup: str = ""
    allow_auto_outgroup: bool = False
    # Physical blocks avoid treating adjacent sites from different chromosomes as a single
    # replicate. The site_count option retains the legacy equal-site-count jackknife.
    block_mode: Literal["fixed_mb", "chromosome", "site_count"] = "chromosome"
    block_size_mb: float = Field(default=5.0, gt=0)
    n_blocks: int = Field(default=20, ge=1)  # used only by legacy site_count mode
    min_effective_blocks: int = Field(default=2, ge=2)
    window_bp: int = 100_000
    min_window_sites: int = 10


@STAGES.register("dstats")
class DStatsStage(Stage):
    name = "dstats"
    config_model = DStatsConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "dstats")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: DStatsConfig = self.config  # type: ignore[assignment]
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        groups = population_indices(geno, labels)
        freqs = fstats.allele_frequencies(geno, groups)

        if not cfg.outgroup and not cfg.allow_auto_outgroup:
            raise StageInputError(
                "D-statistics require an explicit biological outgroup; set stages.dstats.outgroup"
            )
        outgroup = cfg.outgroup or _auto_outgroup(freqs)
        if outgroup not in groups:
            raise StageInputError(f"configured D-statistics outgroup is absent: {outgroup}")
        if len(groups) < 4:
            raise StageInputError(
                "D-statistics require at least four populations including outgroup"
            )
        ingroup = [p for p in groups if p != outgroup]
        fstats_block_mode = {
            "fixed_mb": "fixed_bp",
            "chromosome": "chromosome",
            "site_count": "site_count",
        }[cfg.block_mode]
        block_size_bp = max(1, round(cfg.block_size_mb * 1_000_000))

        rows = []
        for p3 in ingroup:
            for p1, p2 in combinations([p for p in ingroup if p != p3], 2):
                d = fstats.d_statistic(freqs[p1], freqs[p2], freqs[p3], freqs[outgroup],
                                       n_blocks=cfg.n_blocks, chrom=geno.chrom, pos=geno.pos,
                                       block_mode=fstats_block_mode,
                                       block_size_bp=block_size_bp)
                pval = float(2 * norm.sf(abs(d.z))) if np.isfinite(d.z) else float("nan")
                rows.append({"P1": p1, "P2": p2, "P3": p3, "O": outgroup,
                             "D": round(d.estimate, 5), "Z": round(d.z, 3),
                             "se": round(d.se, 5), "n_sites": d.n_sites,
                             "n_blocks": d.n_blocks, "block_mode": cfg.block_mode,
                             "block_size_bp": block_size_bp if cfg.block_mode == "fixed_mb"
                             else None,
                             "p_value": None if pval != pval else float(f"{pval:.3g}")})

        table = pd.DataFrame(rows)
        if table.empty:
            raise StageInputError("no valid D-statistic quartets were generated")
        if int(table["n_blocks"].min()) < cfg.min_effective_blocks:
            raise StageInputError(
                f"D-statistics produced fewer than {cfg.min_effective_blocks} effective blocks"
            )
        table["q_value_bh"] = _bh_qvalues(table["p_value"].to_numpy(dtype=float))
        out = ctx.datastore.path_for(self.name, "dstats.csv")
        table.to_csv(out, index=False)

        top = self._window_scan(ctx, geno, freqs, table, outgroup, cfg)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "dstats", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "dstats", "outgroup": outgroup,
                      "n_quartets": len(table), "top_quartet": top,
                      "block_mode": cfg.block_mode, "block_size_bp": block_size_bp,
                      "legacy_n_blocks": cfg.n_blocks,
                      "multiple_testing": "Benjamini-Hochberg",
                      "median_sites_per_block": float(
                          np.nanmedian(table["n_sites"] / table["n_blocks"])
                      )})
        return StageResult(
            artifacts=[art],
            metrics={"n_quartets": len(table), "outgroup": outgroup,
                     "block_mode": cfg.block_mode,
                     "n_blocks": int(table["n_blocks"].max()) if not table.empty else 0,
                     "max_abs_Z": None if table.empty else round(
                         float(table["Z"].abs().max()), 3)})

    def _window_scan(self, ctx, geno: Genotypes, freqs, table: pd.DataFrame,
                     outgroup: str, cfg: DStatsConfig):
        if table.empty:
            return None
        best = table.iloc[table["Z"].abs().argmax()]
        quartet = (best["P1"], best["P2"], best["P3"], outgroup)
        windows = fstats.f_d_windows(
            freqs, quartet, geno.chrom, geno.pos,
            window_bp=cfg.window_bp, min_sites=cfg.min_window_sites)
        if windows:
            pd.DataFrame(windows).to_csv(
                ctx.datastore.path_for(self.name, "fd_windows.csv"), index=False)
        return {"P1": best["P1"], "P2": best["P2"], "P3": best["P3"], "O": outgroup,
                "D": float(best["D"]), "Z": float(best["Z"])}


def _auto_outgroup(freqs: dict[str, np.ndarray]) -> str:
    """Pick the most genetically divergent population (largest mean allele-freq distance)."""
    pops = list(freqs)
    scores = {}
    for p in pops:
        diffs = [float(np.nanmean(np.abs(freqs[p] - freqs[q]))) for q in pops if q != p]
        scores[p] = float(np.mean(diffs)) if diffs else 0.0
    return max(scores, key=lambda p: scores[p])


def _bh_qvalues(values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, preserving NaNs and original order."""
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(values))
    if not finite.size:
        return result
    ordered = finite[np.argsort(values[finite])]
    adjusted = values[ordered] * finite.size / np.arange(1, finite.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result[ordered] = np.minimum(adjusted, 1.0)
    return result
