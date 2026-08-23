"""Unsupervised clustering of individuals from the genetic-distance matrix.

Agglomerative (average-linkage) clustering over the precomputed distance matrix, sweeping
the number of clusters and choosing the value with the best silhouette. The resulting
clustering is scored against the sample sheet's population labels with the adjusted Rand
index (ARI) — a direct check that unsupervised structure matches known populations.
"""

from __future__ import annotations

import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import adjusted_rand_score, silhouette_score

from canidae.core.errors import StageInputError
from canidae.core.model import ArtifactKind, FileFormat
from canidae.core.registry import STAGES
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult
from canidae.stages.popgen.store import load_sample_labels


class ClusterConfig(StageConfig):
    k_min: int = 2
    k_max: int = 8
    linkage: str = "average"


@STAGES.register("cluster")
class ClusterStage(Stage):
    name = "cluster"
    config_model = ClusterConfig

    def required_inputs(self) -> list[ArtifactSpec]:
        return [
            ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "distance"),
            ArtifactSpec(ArtifactKind.SAMPLE_SHEET, "sample_sheet"),
        ]

    def produced_outputs(self) -> list[ArtifactSpec]:
        return [ArtifactSpec(ArtifactKind.ANALYSIS_RESULT, "cluster")]

    def run(self, ctx: RunContext) -> StageResult:
        cfg: ClusterConfig = self.config  # type: ignore[assignment]
        dist = pd.read_csv(
            ctx.datastore.get(ArtifactKind.ANALYSIS_RESULT, "distance").path, index_col=0
        )
        samples = list(dist.index)
        D = dist.to_numpy(dtype=float)
        n = len(samples)

        k_hi = min(cfg.k_max, n - 1)
        if k_hi < cfg.k_min:
            raise StageInputError(f"too few samples ({n}) to cluster in [{cfg.k_min}, {cfg.k_max}]")

        best = None
        for k in range(cfg.k_min, k_hi + 1):
            labels = AgglomerativeClustering(
                n_clusters=k, metric="precomputed", linkage=cfg.linkage
            ).fit_predict(D)
            score = float(silhouette_score(D, labels, metric="precomputed"))
            if best is None or score > best[1]:
                best = (k, score, labels)
        assert best is not None
        best_k, silhouette, labels = best

        sheet = load_sample_labels(
            ctx.datastore.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path
        )
        pops = [str(sheet.loc[s, "population"]) if s in sheet.index else "unknown" for s in samples]
        ari = float(adjusted_rand_score(pops, labels))

        table = pd.DataFrame({"sample_id": samples, "population": pops, "cluster": labels})
        out = ctx.datastore.path_for(self.name, "clusters.csv")
        table.to_csv(out, index=False)
        art = ctx.datastore.add(
            ArtifactKind.ANALYSIS_RESULT,
            "cluster",
            out,
            fmt=FileFormat.CSV,
            produced_by=self.name,
            metadata={
                "analysis": "cluster",
                "best_k": int(best_k),
                "silhouette": round(silhouette, 4),
                "adjusted_rand_index": round(ari, 4),
                "linkage": cfg.linkage,
            },
        )
        return StageResult(
            artifacts=[art],
            metrics={
                "best_k": int(best_k),
                "silhouette": round(silhouette, 4),
                "adjusted_rand_index": round(ari, 4),
            },
        )
