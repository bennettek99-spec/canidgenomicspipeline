"""Runs of homozygosity (ROH) per individual.

A window-based estimator: each contig is tiled into fixed-size windows; a window is
"homozygous" when its observed heterozygosity is at or below a threshold; consecutive
homozygous windows merge into ROH segments. Reports, per sample, the number of segments,
total ROH length, and F_ROH (ROH fraction of the covered genome) — a standard proxy for
inbreeding / small effective population size. A model-based HMM caller can replace this
behind the same stage later.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import Genotypes, load_genotypes


class RohConfig(StageConfig):
    window_bp: int = 100_000
    max_window_heterozygosity: float = 0.02
    min_segment_bp: int = 200_000


@STAGES.register("roh")
class RohStage(Stage):
    name = "roh"
    config_model = RohConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "roh")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: RohConfig = self.config  # type: ignore[assignment]
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        genome_bp = _covered_genome_bp(geno)
        is_het = geno.calls.is_het()  # (n_variants, n_samples)

        rows = []
        for s, sample in enumerate(geno.samples):
            total_bp, n_seg = _sample_roh(
                geno.chrom, geno.pos, is_het[:, s], cfg)
            rows.append({
                "sample_id": sample,
                "n_roh_segments": n_seg,
                "total_roh_bp": int(total_bp),
                "froh": round(total_bp / genome_bp, 6) if genome_bp else 0.0,
            })
        table = pd.DataFrame(rows)

        out = ctx.datastore.path_for(self.name, "roh.csv")
        table.to_csv(out, index=False)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "roh", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "roh", "window_bp": cfg.window_bp,
                      "genome_bp": int(genome_bp)},
        )
        return StageResult(
            artifacts=[art],
            metrics={"n_samples": geno.n_samples,
                     "mean_froh": round(float(table["froh"].mean()), 6)},
        )


def _covered_genome_bp(geno: Genotypes) -> int:
    total = 0
    for contig in np.unique(geno.chrom):
        p = geno.pos[geno.chrom == contig]
        total += int(p.max() - p.min()) if p.size else 0
    return total


def _sample_roh(chrom: np.ndarray, pos: np.ndarray, het: np.ndarray,
                cfg: RohConfig) -> tuple[int, int]:
    """Return (total_roh_bp, n_segments) for one sample across all contigs."""
    total_bp = 0
    n_segments = 0
    for contig in np.unique(chrom):
        m = chrom == contig
        cpos, chet = pos[m], het[m]
        if cpos.size == 0:
            continue
        start, end = int(cpos.min()), int(cpos.max())
        run_start = None
        for w0 in range(start, end + 1, cfg.window_bp):
            w1 = w0 + cfg.window_bp
            in_win = (cpos >= w0) & (cpos < w1)
            n = int(in_win.sum())
            het_rate = float(chet[in_win].mean()) if n else 0.0
            homozygous = n > 0 and het_rate <= cfg.max_window_heterozygosity
            if homozygous and run_start is None:
                run_start = w0
            elif not homozygous and run_start is not None:
                total_bp, n_segments = _close_run(run_start, w0, cfg, total_bp, n_segments)
                run_start = None
        if run_start is not None:
            total_bp, n_segments = _close_run(run_start, end, cfg, total_bp, n_segments)
    return total_bp, n_segments


def _close_run(start: int, stop: int, cfg: RohConfig, total_bp: int,
               n_segments: int) -> tuple[int, int]:
    length = stop - start
    if length >= cfg.min_segment_bp:
        return total_bp + length, n_segments + 1
    return total_bp, n_segments
