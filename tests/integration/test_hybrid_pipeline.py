"""End-to-end hybrid-canid diagnostic recipes on a synthetic bridge panel.

Runs the two shipped recipes — NYC-style pedigree validation plus breed
assignment, and eastern-coyote three-way admixture — through the real executor,
with no network access and no study data. The assertions check that the stage
DAG wires up, that the estimators land on the fixture's known mixture
fractions, and that the report renders the hybrid sections together with their
limitations.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.model import ArtifactKind
from canidae.pipeline import run_pipeline

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_hybrid_panel import BUG_COLUMN, build_hybrid_panel

pytestmark = [pytest.mark.integration]
R = ArtifactKind.ANALYSIS_RESULT


def _store(cfg: GlobalConfig) -> DataStore:
    return DataStore(cfg.paths.data_root / "store")


def _manifest(store: DataStore, role: str) -> dict:
    art = store.get(R, role)
    return json.loads(Path(art.metadata["manifest"]).read_text(encoding="utf-8"))


# -- NYC-style recipe: reference_mixture -> breed_assign -> report ----------------------


@pytest.fixture(scope="module")
def nyc_run(tmp_path_factory):
    ws = tmp_path_factory.mktemp("hybrid_nyc")
    panel = build_hybrid_panel(ws / "panel")
    queries = [s for s in panel.two_source_query_ids() if s.startswith("NYFIXTURE")]
    cfg = GlobalConfig.load(overrides={
        "project_name": "hybrid-nyc", "paths.root": str(ws),
        "pipeline": ["reference_mixture", "breed_assign", "report"],
        "logging.level": "WARNING",
        "stages.reference_mixture.bridge_vcf": str(panel.bridge_vcf),
        "stages.reference_mixture.reference_genotypes": str(panel.reference_genotypes),
        "stages.reference_mixture.calls_dir": str(panel.calls_dir),
        "stages.reference_mixture.query_samples": queries,
        "stages.reference_mixture.wgs_samples": panel.reference_samples,
        "stages.reference_mixture.coyote_samples": panel.coyote_samples,
        "stages.reference_mixture.expected_f1": "NYFIXTURE_F1",
        "stages.reference_mixture.expected_parent": "NYFIXTURE_PARENT",
        "stages.reference_mixture.expected_offspring": ["NYFIXTURE_OFF1",
                                                        "NYFIXTURE_OFF2"],
        "stages.reference_mixture.f1_target": 0.50,
        "stages.breed_assign.bridge_vcf": str(panel.bridge_vcf),
        "stages.breed_assign.wgs_genotypes": str(panel.wgs_genotypes),
        "stages.breed_assign.calls_dir": str(panel.calls_dir),
        "stages.breed_assign.dog_fractions": str(panel.dog_fractions),
        "stages.breed_assign.locus_filter_json": str(panel.reference_genotypes),
        "stages.breed_assign.query_samples": queries,
        "stages.breed_assign.min_group": 3,
        "stages.report.title": "Hybrid fixture",
        "stages.report.exploratory": True,
        "stages.report.citations": ["nhgri_722g_wgs", "eastern_coyote_radseq"],
    })
    run_pipeline(cfg)
    return cfg, panel, queries


def test_nyc_recipe_produces_all_three_artifacts(nyc_run) -> None:
    cfg, _, _ = nyc_run
    store = _store(cfg)
    assert store.has(R, "reference_mixture")
    assert store.has(R, "breed_assign")
    assert store.has(ArtifactKind.REPORT, "html")


def test_nyc_recipe_recovers_the_known_pedigree(nyc_run) -> None:
    cfg, panel, queries = nyc_run
    table = pd.read_csv(_store(cfg).get(R, "reference_mixture").path)
    table = table.set_index("sample_id")
    for sample_id in queries:
        expected = panel.truth[sample_id].f_dog
        assert table.loc[sample_id, "dog_fraction"] == pytest.approx(expected, abs=0.06)
    # F1 above both offspring, offspring above the backcrossed parent.
    assert (
        table.loc["NYFIXTURE_F1", "dog_fraction"]
        > table.loc["NYFIXTURE_OFF1", "dog_fraction"]
        > table.loc["NYFIXTURE_PARENT", "dog_fraction"]
    )
    verdict = _manifest(_store(cfg), "reference_mixture")["validation"]
    assert verdict["passed"] is True


def test_breed_assign_finds_the_simulated_dog_source(nyc_run) -> None:
    cfg, _panel, _ = nyc_run
    store = _store(cfg)
    manifest = _manifest(store, "breed_assign")
    # Every query was drawn from a Beagle dog source.
    assert manifest["results"]["NYFIXTURE_F1"]["best_breed"] == "Beagle"
    # Only the sample with substantial dog ancestry clears the any-dog baseline.
    assert manifest["results"]["NYFIXTURE_F1"]["single_breed_supported"] is True
    assert manifest["results"]["NYFIXTURE_PARENT"]["single_breed_supported"] is False
    # Leave-one-out calibration is reported and strong on well-separated panels.
    assert manifest["calibration"]["top1_accuracy"] > 0.8
    assert manifest["limitations"]

    scores = pd.read_csv(store.get(R, "breed_assign").metadata["scores_csv"])
    assert set(scores["sample_id"]) <= set(manifest["results"])
    assert scores["rank"].min() == 1


def test_breed_assign_consumes_the_upstream_mixture_artifact(nyc_run) -> None:
    """breed_assign declares reference_mixture as an optional upstream input."""
    from canidae.stages.hybrid.breed_assign import BreedAssignStage

    specs = BreedAssignStage.required_inputs(BreedAssignStage.__new__(BreedAssignStage))
    assert any(s.role == "reference_mixture" and s.optional for s in specs)


def test_report_renders_the_hybrid_sections(nyc_run) -> None:
    cfg, _, _ = nyc_run
    html = _store(cfg).get(ArtifactKind.REPORT, "html").path.read_text(encoding="utf-8")
    assert "Hybrid-canid diagnostics" in html
    assert "Two-source mixture" in html
    assert "Breed assignment of dog component" in html
    assert "Pedigree-style validation: PASSED." in html
    # The honest-labelling requirement: never present this as WGS ancestry.
    assert "not whole-genome ancestry" in html.lower()
    # A hybrid-only recipe has no PCA/FST/diversity artifacts to plot.
    assert "PCA scatter" not in html


def test_report_cites_the_configured_data_sources(nyc_run) -> None:
    """A hybrid recipe produces no callset, so citations are its only sources block."""
    cfg, _, _ = nyc_run
    html = _store(cfg).get(ArtifactKind.REPORT, "html").path.read_text(encoding="utf-8")
    assert "Source accessions and citations" in html
    assert "Published data sources" in html
    assert "PRJNA448733" in html
    assert "PRJNA857904" in html
    assert "10.5061/dryad.7f9q2cd" in html
    # Unrecorded identifiers are labelled, never invented.
    assert "not recorded" in html


# -- Eastern-coyote recipe: multiway_admixture -> report --------------------------------


@pytest.fixture(scope="module")
def eastern_run(tmp_path_factory):
    ws = tmp_path_factory.mktemp("hybrid_eastern")
    panel = build_hybrid_panel(ws / "panel")
    cfg = GlobalConfig.load(overrides={
        "project_name": "hybrid-eastern", "paths.root": str(ws),
        "pipeline": ["multiway_admixture", "report"],
        "logging.level": "WARNING",
        "stages.multiway_admixture.bridge_vcf": str(panel.bridge_vcf),
        "stages.multiway_admixture.wgs_genotypes": str(panel.wgs_genotypes),
        "stages.multiway_admixture.western_ids": panel.query_ids(region="western"),
        "stages.multiway_admixture.bootstrap_n": 25,
        "stages.multiway_admixture.min_group": 3,
        "stages.report.title": "Eastern fixture",
        "stages.report.exploratory": True,
    })
    run_pipeline(cfg)
    return cfg, panel


def test_multiway_recovers_the_known_wolf_fraction(eastern_run) -> None:
    cfg, panel = eastern_run
    table = pd.read_csv(_store(cfg).get(R, "multiway_admixture").path).set_index("sample_id")
    eastern = [s for s in panel.query_ids(region="eastern") if s.startswith("EC")]
    observed = table.loc[eastern, "f_wolf"].mean()
    expected = panel.truth[eastern[0]].f_wolf
    assert observed == pytest.approx(expected, abs=0.08)
    # Pure-coyote western controls must not pick up spurious ancestry.
    western = panel.query_ids(region="western")
    assert table.loc[western, "f_wolf"].max() < 0.05
    assert table.loc[western, "f_dog"].max() < 0.05


def test_multiway_bootstrap_intervals_are_well_formed(eastern_run) -> None:
    cfg, _panel = eastern_run
    table = pd.read_csv(_store(cfg).get(R, "multiway_admixture").path).set_index("sample_id")
    for column in ("f_wolf", "f_dog"):
        low = table[f"{column}_ci_lo"]
        high = table[f"{column}_ci_hi"]
        assert (low <= high).all()
        assert (low >= 0.0).all() and (high <= 1.0).all()
    manifest = _manifest(_store(cfg), "multiway_admixture")
    for cohort in ("eastern", "western_control"):
        estimate = manifest["cohort_estimates"][cohort]
        assert len(estimate["f_wolf_ci"]) == 2
        assert estimate["f_wolf_ci"][0] <= estimate["f_wolf_ci"][1]
    assert manifest["method"]["bootstrap_n"] == 25


def test_multiway_excludes_the_mislabeled_bridge_columns(eastern_run) -> None:
    """The known-bad WGS columns in the bridge VCF must never enter the query set."""
    cfg, _ = eastern_run
    table = pd.read_csv(_store(cfg).get(R, "multiway_admixture").path)
    assert BUG_COLUMN not in set(table["sample_id"])


def test_multiway_keeps_them_when_the_guard_is_disabled(tmp_path) -> None:
    """Turning the guard off is the only way the mislabeled column appears."""
    from canidae.core.executor import build_context
    from canidae.core.provenance import ProvenanceWriter
    from canidae.core.runtime import LocalRunner
    from canidae.stages.hybrid.multiway_admixture import (
        MultiwayAdmixtureConfig,
        MultiwayAdmixtureStage,
    )

    panel = build_hybrid_panel(tmp_path / "panel")
    cfg = GlobalConfig.load(include_defaults=True, overrides={
        "project_name": "guard-off", "paths.root": str(tmp_path),
        "pipeline": ["multiway_admixture"], "logging.level": "WARNING",
    })
    store = DataStore(tmp_path / "store")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    ctx = build_context(
        cfg, store, LocalRunner(),
        ProvenanceWriter(run_dir, config_digest=cfg.digest(), seed=1), run_dir=run_dir,
    )
    result = MultiwayAdmixtureStage(MultiwayAdmixtureConfig(
        bridge_vcf=panel.bridge_vcf,
        wgs_genotypes=panel.wgs_genotypes,
        bootstrap_n=10,
        exclude_bridge_wgs_bug_ids=False,
    )).run(ctx)
    table = pd.read_csv(result.artifacts[0].path)
    assert BUG_COLUMN in set(table["sample_id"])


def test_eastern_report_renders_the_multiway_section(eastern_run) -> None:
    cfg, _ = eastern_run
    html = _store(cfg).get(ArtifactKind.REPORT, "html").path.read_text(encoding="utf-8")
    assert "Three-way coyote / wolf / dog admixture" in html
    assert "Hybrid-canid diagnostics" in html
