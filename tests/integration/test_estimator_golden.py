"""Golden regression for the core population-genomics estimators.

The existing integration tests assert biological ordering and detection. This
fixture also pins a compact set of numerical summaries so a refactor cannot
silently change several downstream estimators at once.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.model import ArtifactKind
from canidae.pipeline import run_pipeline

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_cohort import simulate_cohort

pytestmark = [pytest.mark.integration, pytest.mark.slow]

ROOT = Path(__file__).resolve().parents[2]
GOLDEN = ROOT / "tests" / "golden" / "population_genomics.json"
RESULT = ArtifactKind.ANALYSIS_RESULT
PIPELINE = [
    "ingest", "qc", "load_genotypes", "distance", "pca", "fst", "diversity",
    "admixture", "cluster", "roh", "geography",
]


def _round(value: object, places: int = 6) -> object:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (int, bool, str)):
        return value
    return round(float(value), places)


@pytest.fixture(scope="module")
def observed(tmp_path_factory) -> dict:
    workspace = tmp_path_factory.mktemp("population-golden")
    vcf, sheet = simulate_cohort(workspace / "input", seed=7)
    cfg = GlobalConfig.load(overrides={
        "project_name": "population-golden",
        "paths.root": str(workspace),
        "pipeline": PIPELINE,
        "executor.max_workers": 2,
        "logging.level": "WARNING",
        "stages.ingest.sample_sheet": str(sheet),
        "stages.ingest.callset": str(vcf),
        "stages.admixture.backend": "nmf",
    })
    report = run_pipeline(cfg)
    assert report.ok, report.failed
    store = DataStore(cfg.paths.data_root / "store")

    pca_art = store.get(RESULT, "pca")
    pca = pd.read_csv(pca_art.path)
    diversity = pd.read_csv(store.get(RESULT, "diversity").path)
    admixture_art = store.get(RESULT, "admixture")
    admixture = pd.read_csv(admixture_art.path)
    q_cols = [column for column in admixture.columns if column.startswith("Q")]
    roh = pd.read_csv(store.get(RESULT, "roh").path).sort_values("sample_id")
    distance = pd.read_csv(store.get(RESULT, "distance").path, index_col=0).to_numpy()
    upper = distance[np.triu_indices(distance.shape[0], k=1)]
    cluster = store.get(RESULT, "cluster")
    geography = store.get(RESULT, "geography")

    return {
        "distance": {
            "n_sites": store.get(RESULT, "distance").metadata["n_sites"],
            "mean_pairwise": _round(np.mean(upper)),
        },
        "fst": {
            pair: _round(value)
            for pair, value in sorted(store.get(RESULT, "fst").metadata["pairwise"].items())
        },
        "pca": {
            "explained_variance_ratio": pca_art.metadata["explained_variance_ratio"],
            "population_pc1": {
                population: _round(value)
                for population, value in pca.groupby("population")["PC1"].mean().items()
            },
        },
        "diversity": {
            row.population: {
                "pi": _round(row.pi, 8),
                "mean_heterozygosity": _round(row.mean_heterozygosity, 6),
            }
            for row in diversity.itertuples(index=False)
        },
        "admixture": {
            "best_k": admixture_art.metadata["best_k"],
            "population_q": {
                population: [_round(value) for value in values]
                for population, values in admixture.groupby("population")[q_cols].mean().iterrows()
            },
        },
        "cluster": {
            "best_k": cluster.metadata["best_k"],
            "silhouette": _round(cluster.metadata["silhouette"], 4),
            "adjusted_rand_index": _round(cluster.metadata["adjusted_rand_index"], 4),
        },
        "roh": {
            row.sample_id: {
                "n_roh_segments": int(row.n_roh_segments),
                "total_roh_bp": int(row.total_roh_bp),
                "froh": _round(row.froh, 6),
            }
            for row in roh.itertuples(index=False)
        },
        "geography": {
            "n_localities": geography.metadata["n_localities"],
            "mantel_r": _round(geography.metadata["mantel_r"], 4),
        },
    }


def test_population_estimators_match_golden(observed: dict) -> None:
    if os.environ.get("CANIDAE_UPDATE_GOLDEN"):
        GOLDEN.write_text(json.dumps(observed, indent=2) + "\n", encoding="utf-8")
        pytest.skip(f"golden snapshot rewritten: {GOLDEN}")
    assert GOLDEN.exists(), f"missing golden snapshot; regenerate {GOLDEN}"
    assert observed == json.loads(GOLDEN.read_text(encoding="utf-8"))
