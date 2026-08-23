"""Geographic / spatial-genetic-structure stage.

Combines the genetic-distance matrix with sample localities to produce:

* an **isolation-by-distance** test (Mantel r + permutation p between geographic and
  genetic distance),
* an **isolation-by-distance** scatter table (pairwise geographic vs genetic distance),
* a **regional ancestry** summary (mean admixture proportions per region, when admixture
  results are available), and
* a **localities** table used to draw the sample map.

The map figure and IBD scatter are rendered by the reporting stage from these tables.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.geographic.distance import haversine_matrix, mantel_test


class GeographyConfig(StageConfig):
    permutations: int = 999


@STAGES.register("geography")
class GeographyStage(Stage):
    name = "geography"
    config_model = GeographyConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "distance"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
            ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "admixture", optional=True),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "geography")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: GeographyConfig = self.config  # type: ignore[assignment]
        ds = ctx.datastore
        stage_dir = ds.stage_dir(self.name)

        genetic = pd.read_csv(ds.get(ArtifactKind.ANALYSIS_RESULT, "distance").path, index_col=0)
        sheet = pd.read_csv(ds.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path, dtype=str)
        localities = _localities(sheet, list(genetic.index))
        localities.to_csv(stage_dir / "localities.csv", index=False)

        mantel_r, mantel_p, n_loc = self._isolation_by_distance(
            genetic, localities, stage_dir, cfg, seed=ctx.config.seed
        )
        regions = self._regional_ancestry(ds, localities, stage_dir)

        summary = pd.DataFrame(
            [
                {
                    "n_localities": n_loc,
                    "mantel_r": None if mantel_r != mantel_r else round(mantel_r, 4),
                    "mantel_p": None if mantel_p != mantel_p else round(mantel_p, 4),
                    "n_regions": 0 if regions is None else len(regions),
                }
            ]
        )
        out = stage_dir / "geography_summary.csv"
        summary.to_csv(out, index=False)

        art = ds.add(
            ArtifactKind.ANALYSIS_RESULT,
            "geography",
            out,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "geography",
                "mantel_r": None if mantel_r != mantel_r else round(mantel_r, 4),
                "mantel_p": None if mantel_p != mantel_p else round(mantel_p, 4),
                "n_localities": n_loc,
                "has_regional_ancestry": regions is not None,
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={
                "n_localities": n_loc,
                "mantel_r": None if mantel_r != mantel_r else round(mantel_r, 4),
            },
        )

    def _isolation_by_distance(
        self, genetic: pd.DataFrame, localities: pd.DataFrame, stage_dir, cfg, *, seed: int
    ) -> tuple[float, float, int]:
        loc = localities.dropna(subset=["latitude", "longitude"])
        ids = [s for s in genetic.index if s in set(loc["sample_id"])]
        if len(ids) < 3:
            return float("nan"), float("nan"), len(ids)

        loc = loc.set_index("sample_id").loc[ids]
        geo = haversine_matrix(loc["latitude"].to_numpy(float), loc["longitude"].to_numpy(float))
        gen = genetic.loc[ids, ids].to_numpy(float)
        r, p = mantel_test(gen, geo, permutations=cfg.permutations, seed=seed)

        iu = np.triu_indices(len(ids), k=1)
        pairs = pd.DataFrame({"geo_km": geo[iu], "genetic_distance": gen[iu]})
        pairs.to_csv(stage_dir / "ibd_pairs.csv", index=False)
        return r, p, len(ids)

    def _regional_ancestry(self, ds, localities: pd.DataFrame, stage_dir):
        if not ds.has(ArtifactKind.ANALYSIS_RESULT, "admixture"):
            return None
        q = pd.read_csv(ds.get(ArtifactKind.ANALYSIS_RESULT, "admixture").path)
        q_cols = [c for c in q.columns if c.startswith("Q")]
        merged = q.merge(localities[["sample_id", "region"]], on="sample_id", how="left")
        merged["region"] = merged["region"].fillna("unknown")
        regional = merged.groupby("region")[q_cols].mean().round(4).reset_index()
        regional.to_csv(stage_dir / "regional_ancestry.csv", index=False)
        return regional


def _localities(sheet: pd.DataFrame, order: list[str]) -> pd.DataFrame:
    cols = {c: c for c in ("sample_id", "population", "region", "latitude", "longitude")}
    for needed in cols:
        if needed not in sheet.columns:
            sheet[needed] = ""
    df = sheet[list(cols)].copy()
    df = df[df["sample_id"].isin(order)]
    for axis in ("latitude", "longitude"):
        df[axis] = pd.to_numeric(df[axis].replace("", np.nan), errors="coerce")
    # preserve the genetic-matrix order
    df["__order"] = df["sample_id"].map({s: i for i, s in enumerate(order)})
    return df.sort_values("__order").drop(columns="__order").reset_index(drop=True)
