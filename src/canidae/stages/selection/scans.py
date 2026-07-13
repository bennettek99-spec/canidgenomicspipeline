"""Selection-scan stage: PBS and windowed Tajima's D.

Computes two complementary, dependency-light selection statistics in genomic windows:

* **PBS** (population branch statistic, Yi et al. 2010) — the branch length specific to a
  focal population from pairwise F_ST, highlighting lineage-specific allele-frequency change.
* **Tajima's D** per window — an allele-frequency-spectrum summary sensitive to sweeps
  (negative) and balancing selection / structure (positive).

Haplotype-based statistics (iHS/nSL/XP-EHH via scikit-allel or selscan) are natural
extensions on phased data. All CPU; no GPU used.
"""

from __future__ import annotations

import allel
import numpy as np
import pandas as pd
from pydantic import Field

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import load_genotypes, load_sample_labels, population_indices


class SelectionConfig(StageConfig):
    focal_population: str = ""       # PBS focal pop; auto (first) if empty
    reference_populations: list[str] = Field(default_factory=list)  # auto if empty
    window_bp: int = 100_000


@STAGES.register("selection")
class SelectionStage(Stage):
    name = "selection"
    config_model = SelectionConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.GENOTYPES, "genotypes"),
            ArtifactSpec(ArtifactKind.GENOTYPES, "analysis_genotypes", optional=True),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "selection")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: SelectionConfig = self.config  # type: ignore[assignment]
        role = "analysis_genotypes" if ctx.datastore.has(
            ArtifactKind.GENOTYPES, "analysis_genotypes"
        ) else "genotypes"
        geno = load_genotypes(ctx.datastore.get(ArtifactKind.GENOTYPES, role).path)
        labels = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
        groups = population_indices(geno, labels)

        focal, ref_b, ref_c = self._triad(cfg, groups)
        ac = {p: geno.calls.count_alleles(subpop=groups[p]) for p in (focal, ref_b, ref_c)}

        rows = []
        for contig in dict.fromkeys(geno.chrom):
            m = geno.chrom == contig
            pos = geno.pos[m]
            if pos.size < 2:
                continue
            rows.extend(self._contig_windows(str(contig), pos, {p: ac[p][m] for p in ac},
                                             focal, ref_b, ref_c, cfg.window_bp))
        table = pd.DataFrame(rows)
        if table.empty:
            raise StageInputError("no windows computed (too few sites)")

        out = ctx.datastore.path_for(self.name, "selection_windows.csv")
        table.to_csv(out, index=False)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT, "selection", out, fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={"analysis": "selection", "focal": focal,
                      "references": [ref_b, ref_c], "window_bp": cfg.window_bp,
                      "n_windows": len(table),
                      "max_pbs": round(float(np.nanmax(table["pbs"])), 4)})
        return StageResult(
            artifacts=[art],
            metrics={"n_windows": len(table),
                     "max_pbs": round(float(np.nanmax(table["pbs"])), 4)})

    def _triad(self, cfg: SelectionConfig, groups) -> tuple[str, str, str]:
        pops = list(groups)
        if len(pops) < 3:
            raise StageInputError("PBS needs >= 3 populations")
        focal = cfg.focal_population or pops[0]
        refs = cfg.reference_populations or [p for p in pops if p != focal][:2]
        if focal not in groups or any(r not in groups for r in refs) or len(refs) < 2:
            raise StageInputError("invalid focal/reference populations for PBS")
        return focal, refs[0], refs[1]

    def _contig_windows(self, contig, pos, ac, focal, ref_b, ref_c, window_bp):
        fst_ab = _windowed_fst(pos, ac[focal], ac[ref_b], window_bp)
        fst_ac = _windowed_fst(pos, ac[focal], ac[ref_c], window_bp)
        fst_bc = _windowed_fst(pos, ac[ref_b], ac[ref_c], window_bp)
        td, windows, _ = allel.windowed_tajima_d(pos, ac[focal], size=window_bp)
        rows = []
        for i, (start, stop) in enumerate(windows):
            pbs = _pbs(fst_ab[i], fst_ac[i], fst_bc[i])
            rows.append({"chrom": contig, "start": int(start), "end": int(stop),
                         "pbs": round(pbs, 5) if np.isfinite(pbs) else np.nan,
                         "tajima_d": round(float(td[i]), 5) if np.isfinite(td[i])
                         else np.nan})
        return rows


def _windowed_fst(pos, ac1, ac2, window_bp) -> np.ndarray:
    fst, _, _ = allel.windowed_hudson_fst(pos, ac1, ac2, size=window_bp)
    return np.asarray(fst, dtype=float)


def _pbs(fst_ab: float, fst_ac: float, fst_bc: float) -> float:
    def t(f: float) -> float:
        f = min(max(f, 0.0), 0.9999)
        return -np.log(1.0 - f)
    if not all(np.isfinite([fst_ab, fst_ac, fst_bc])):
        return float("nan")
    return (t(fst_ab) + t(fst_ac) - t(fst_bc)) / 2.0
