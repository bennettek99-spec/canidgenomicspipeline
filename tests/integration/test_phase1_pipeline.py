"""End-to-end test: ingest -> qc -> load -> distance/pca/fst/diversity/admixture/cluster/
roh -> geography -> report.

Runs the full Phase-1 + Phase-2 population-genomics pipeline on an msprime-simulated cohort
whose structure is known, and asserts the pipeline recovers that ground truth. This is the
validation strategy the architecture calls for: simulate a known truth, confirm the
statistics reproduce it.
"""

from __future__ import annotations

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

PIPELINE = [
    "ingest", "qc", "load_genotypes", "distance", "pca", "fst", "diversity",
    "admixture", "cluster", "roh", "geography", "report",
]
R = ArtifactKind.ANALYSIS_RESULT


@pytest.fixture(scope="module")
def phase1_run(tmp_path_factory) -> tuple[GlobalConfig, object]:
    workspace = tmp_path_factory.mktemp("phase2")
    vcf, sheet = simulate_cohort(workspace / "input", seed=7)
    cfg = GlobalConfig.load(overrides={
        "project_name": "canis-test",
        "paths.root": str(workspace),
        "pipeline": PIPELINE,
        "executor.max_workers": 3,
        "logging.level": "WARNING",
        "stages.ingest.sample_sheet": str(sheet),
        "stages.ingest.callset": str(vcf),
        "stages.admixture.backend": "nmf",
    })
    report = run_pipeline(cfg)
    return cfg, report


def _store(cfg: GlobalConfig) -> DataStore:
    return DataStore(cfg.paths.data_root / "store")


def _csv(cfg, role) -> pd.DataFrame:
    return pd.read_csv(_store(cfg).get(R, role).path)


# -- Phase 1 (unchanged) ---------------------------------------------------------------


def test_pipeline_runs_all_stages(phase1_run) -> None:
    _cfg, report = phase1_run
    assert report.ok, report.failed
    assert set(report.executed) == set(PIPELINE)


def test_report_has_all_sections(phase1_run) -> None:
    cfg, _ = phase1_run
    html = _store(cfg).get(ArtifactKind.REPORT, "html").path.read_text(encoding="utf-8")
    for section in ("Principal component analysis", "Admixture", "Population clustering",
                    "Genetic diversity", "Runs of homozygosity", "Geography"):
        assert section in html, section


def test_fst_recovers_known_topology(phase1_run) -> None:
    cfg, _ = phase1_run
    fst = pd.read_csv(_store(cfg).get(R, "fst").path, index_col=0)
    assert fst.loc["coyote", "wolf"] > fst.loc["wolf", "dog"]
    assert fst.loc["coyote", "dog"] > fst.loc["wolf", "dog"]
    assert (fst.to_numpy() >= -0.02).all()


def test_pca_separates_coyote(phase1_run) -> None:
    cfg, _ = phase1_run
    pca_art = _store(cfg).get(R, "pca")
    means = pd.read_csv(pca_art.path).groupby("population")["PC1"].mean()
    assert abs(means["coyote"] - means["wolf"]) > abs(means["dog"] - means["wolf"])
    assert pca_art.metadata["explained_variance_ratio"][0] > 0.1


def test_diversity_positive(phase1_run) -> None:
    cfg, _ = phase1_run
    div = _csv(cfg, "diversity")
    assert set(div["population"]) == {"wolf", "dog", "coyote"}
    assert (div["pi"] > 0).all()


# -- Phase 2 ---------------------------------------------------------------------------


def test_admixture_separates_coyote(phase1_run) -> None:
    cfg, _ = phase1_run
    art = _store(cfg).get(R, "admixture")
    q = pd.read_csv(art.path)
    q_cols = [c for c in q.columns if c.startswith("Q")]
    assert np.allclose(q[q_cols].sum(axis=1), 1.0, atol=1e-6)   # proportions sum to 1
    assert 2 <= art.metadata["best_k"] <= 6
    q["dom"] = q[q_cols].to_numpy().argmax(axis=1)
    dom = q.groupby("population")["dom"].agg(lambda s: s.mode().iloc[0])
    assert dom["coyote"] != dom["wolf"]  # coyote's dominant ancestry differs from wolf's


def test_clustering_isolates_coyote(phase1_run) -> None:
    cfg, _ = phase1_run
    art = _store(cfg).get(R, "cluster")
    cl = pd.read_csv(art.path)
    coyote = set(cl.loc[cl["population"] == "coyote", "cluster"])
    others = set(cl.loc[cl["population"] != "coyote", "cluster"])
    assert len(coyote) == 1 and coyote.isdisjoint(others)
    assert art.metadata["adjusted_rand_index"] > 0.3


def test_roh_fraction_valid(phase1_run) -> None:
    cfg, _ = phase1_run
    roh = _csv(cfg, "roh")
    assert {"sample_id", "n_roh_segments", "total_roh_bp", "froh"} <= set(roh.columns)
    assert ((roh["froh"] >= 0) & (roh["froh"] <= 1)).all()
    assert len(roh) == 18


def test_geography_outputs(phase1_run) -> None:
    cfg, _ = phase1_run
    art = _store(cfg).get(R, "geography")
    assert art.metadata["n_localities"] == 18
    assert np.isfinite(art.metadata["mantel_r"])
    geo_dir = art.path.parent
    assert (geo_dir / "localities.csv").exists()
    assert (geo_dir / "ibd_pairs.csv").exists()
    assert (geo_dir / "regional_ancestry.csv").exists()


def test_resume_skips_second_run(phase1_run) -> None:
    cfg, _ = phase1_run
    second = run_pipeline(cfg)
    assert set(second.skipped) == set(PIPELINE)
    assert not second.executed
