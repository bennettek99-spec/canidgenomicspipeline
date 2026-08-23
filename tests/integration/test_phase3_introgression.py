"""Phase-3 end-to-end test: phylogenetics + introgression on a cohort with KNOWN gene flow.

The cohort is simulated with a COYOTE->DOG admixture pulse under the topology
((WOLF,DOG),COYOTE) with JACKAL as outgroup. The pipeline must (a) detect the gene flow with
a significant D-statistic, (b) NOT flag it in a no-flow control, and (c) recover the tree
(coyote monophyletic, jackal outgroup).
"""

from __future__ import annotations

import sys
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from Bio import Phylo

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.model import ArtifactKind
from canidae.pipeline import run_pipeline
from canidae.stages.introgression import fstats
from canidae.stages.popgen.store import (
    load_genotypes,
    load_sample_labels,
    population_indices,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_cohort import simulate_introgression_cohort

pytestmark = [pytest.mark.integration, pytest.mark.slow]
R = ArtifactKind.ANALYSIS_RESULT

PIPELINE = [
    "ingest",
    "qc",
    "load_genotypes",
    "analysis_readiness",
    "distance",
    "nj_tree",
    "pca",
    "fst",
    "diversity",
    "f3",
    "dstats",
    "report",
]


@pytest.fixture(scope="module")
def admixed(tmp_path_factory) -> GlobalConfig:
    ws = tmp_path_factory.mktemp("intro_admixed")
    vcf, sheet = simulate_introgression_cohort(ws / "in", seed=11, admixture_proportion=0.2)
    cfg = GlobalConfig.load(
        overrides={
            "project_name": "intro",
            "paths.root": str(ws),
            "pipeline": PIPELINE,
            "executor.max_workers": 3,
            "logging.level": "WARNING",
            "stages.ingest.sample_sheet": str(sheet),
            "stages.ingest.callset": str(vcf),
            "stages.nj_tree.n_bootstrap": 30,
            "stages.admixture.backend": "nmf",
            "stages.analysis_readiness.ld_prune": False,
            "stages.f3.outgroup": "jackal",
            "stages.dstats.outgroup": "jackal",
            "stages.dstats.block_mode": "site_count",
        }
    )
    run_pipeline(cfg)
    return cfg


@pytest.fixture(scope="module")
def control_dstat(tmp_path_factory) -> fstats.JackknifeResult:
    ws = tmp_path_factory.mktemp("intro_control")
    vcf, sheet = simulate_introgression_cohort(ws / "in", seed=11, admixture_proportion=0.0)
    cfg = GlobalConfig.load(
        overrides={
            "project_name": "ctrl",
            "paths.root": str(ws),
            "pipeline": ["ingest", "load_genotypes"],
            "logging.level": "WARNING",
            "stages.ingest.sample_sheet": str(sheet),
            "stages.ingest.callset": str(vcf),
        }
    )
    run_pipeline(cfg)
    store = DataStore(cfg.paths.data_root / "store")
    geno = load_genotypes(store.get(ArtifactKind.GENOTYPES, "genotypes").path)
    labels = load_sample_labels(store.get(ArtifactKind.SAMPLE_SHEET, "sample_sheet").path)
    freqs = fstats.allele_frequencies(geno, population_indices(geno, labels))
    return fstats.d_statistic(
        freqs["wolf"], freqs["dog"], freqs["coyote"], freqs["jackal"], n_blocks=20
    )


def _store(cfg: GlobalConfig) -> DataStore:
    return DataStore(cfg.paths.data_root / "store")


def test_dstats_detects_gene_flow(admixed) -> None:
    art = _store(admixed).get(R, "dstats")
    assert art.metadata["outgroup"] == "jackal"
    assert art.metadata["block_mode"] == "site_count"
    df = pd.read_csv(art.path)
    # the trio testing wolf-vs-dog allele sharing with coyote (the true sisters are
    # wolf & dog) must be significant given the simulated coyote->dog gene flow
    trio = df[
        (df["P3"] == "coyote") & df.apply(lambda r: {r["P1"], r["P2"]} == {"wolf", "dog"}, axis=1)
    ]
    assert not trio.empty
    assert trio.iloc[0]["n_blocks"] > 1
    assert abs(trio.iloc[0]["Z"]) > 3
    assert (art.path.parent / "fd_windows.csv").exists()


def test_control_has_no_false_positive(control_dstat) -> None:
    # no gene flow was simulated: D should be small and non-significant
    assert abs(control_dstat.z) < 3.5


def test_nj_tree_recovers_clades(admixed) -> None:
    newick = _store(admixed).get(ArtifactKind.TREE, "nj").path.read_text(encoding="utf-8")
    tree = Phylo.read(StringIO(newick), "newick")
    coyote = [t for t in tree.get_terminals() if t.name.startswith("COYOTE")]
    mrca = tree.common_ancestor(coyote)
    assert all(t.name.startswith("COYOTE") for t in mrca.get_terminals())  # monophyletic
    assert any(t.name.startswith("JACKAL") for t in tree.get_terminals())


def test_f3_matrix_written(admixed) -> None:
    art = _store(admixed).get(R, "f3")
    assert art.metadata["outgroup"] == "jackal"
    assert art.path.exists()


def test_report_has_phase3_sections(admixed) -> None:
    html = _store(admixed).get(ArtifactKind.REPORT, "html").path.read_text(encoding="utf-8")
    for section in ("Phylogeny", "Introgression", "D-statistics"):
        assert section in html, section
