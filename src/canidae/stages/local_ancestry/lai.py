"""Local-ancestry inference stage.

For each individual in the target populations, paints the genome into windows of source
ancestry using the HMM in :mod:`hmm`, given reference source populations. Emits per-window
ancestry calls (for karyograms) and per-individual genome-wide ancestry fractions.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pydantic import Field

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.introgression import fstats
from canidae.stages.local_ancestry.hmm import infer_local_ancestry
from canidae.stages.popgen.store import load_genotypes, load_sample_labels, population_indices


class LocalAncestryConfig(StageConfig):
    sources: list[str] = Field(default_factory=list)   # reference source populations
    targets: list[str] = Field(default_factory=list)   # default: all non-source pops
    window_bp: int = 200_000
    switch_prob: float = 0.02
    min_sites: int = 5


@STAGES.register("local_ancestry")
class LocalAncestryStage(Stage):
    name = "local_ancestry"
    config_model = LocalAncestryConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "local_ancestry")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: LocalAncestryConfig = self.config  # type: ignore[assignment]
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, "genotypes").path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        groups = population_indices(geno, labels)

        sources = cfg.sources or _default_sources(groups)
        missing = [s for s in sources if s not in groups]
        if len(sources) < 2 or missing:
            raise StageInputError(
                f"local ancestry needs >=2 source populations present; missing {missing}")
        targets = cfg.targets or [p for p in groups if p not in sources]
        if not targets:
            raise StageInputError("no target populations to paint")

        source_freqs = {
            s: np.where(np.isfinite(f), f, 0.5)
            for s, f in fstats.allele_frequencies(geno, {s: groups[s] for s in sources}).items()
        }
        dosage = geno.calls.to_n_alt(fill=-1)  # (n_sites, n_samples)

        window_rows, summary_rows = [], []
        for pop in targets:
            for idx in groups[pop]:
                sample = str(geno.samples[idx])
                calls = infer_local_ancestry(
                    dosage[:, idx], geno.chrom, geno.pos, source_freqs,
                    window_bp=cfg.window_bp, switch_prob=cfg.switch_prob,
                    min_sites=cfg.min_sites)
                for c in calls:
                    window_rows.append({"sample_id": sample, "population": pop, **c})
                summary_rows.append(_fractions(sample, pop, calls, sources))

        summary = pd.DataFrame(summary_rows)
        stage_dir = ctx.datastore.stage_dir(self.name)
        pd.DataFrame(window_rows).to_csv(stage_dir / "local_ancestry_windows.csv",
                                         index=False)
        out = stage_dir / "local_ancestry_fractions.csv"
        summary.to_csv(out, index=False)

        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "local_ancestry", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "local_ancestry", "sources": sources,
                      "targets": targets, "window_bp": cfg.window_bp,
                      "chromosome_reset": True,
                      "processing_order": "sequential_chromosomes"})
        return StageResult(
            artifacts=[art],
            metrics={"n_targets": len(summary), "sources": sources})


def _default_sources(groups: dict[str, list[int]]) -> list[str]:
    # heuristic: the two largest populations act as reference panels
    return sorted(groups, key=lambda p: -len(groups[p]))[:2]


def _fractions(sample: str, pop: str, calls: list[dict], sources: list[str]) -> dict:
    total = len(calls)
    counts = {s: 0 for s in sources}
    for c in calls:
        counts[c["ancestry"]] = counts.get(c["ancestry"], 0) + 1
    row = {"sample_id": sample, "population": pop, "n_windows": total}
    for s in sources:
        row[f"frac_{s}"] = round(counts[s] / total, 4) if total else 0.0
    return row
