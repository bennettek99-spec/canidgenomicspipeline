"""Phase 4-6 end-to-end tests: harmonization, local ancestry, selection/demography/ancient.

All run on msprime-simulated cohorts and assert the pipeline recovers known structure. The
read-processing (align/call/joint) and SV stages are external-Linux-tool wrappers verified by
their command-builder unit tests, not here.
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
from canidae.stages.processing.vcf_io import read_biallelic_snps, write_minimal_vcf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_cohort import simulate_cohort, simulate_introgression_cohort

pytestmark = [pytest.mark.integration, pytest.mark.slow]
R = ArtifactKind.ANALYSIS_RESULT


def _store(cfg: GlobalConfig) -> DataStore:
    return DataStore(cfg.paths.data_root / "store")


# -- Phase 4: cross-study harmonization ------------------------------------------------


def _split_into_studies(vcf: Path, sheet: Path, out: Path):
    """Split a cohort into two 'studies' (by sample) sharing all sites."""
    sites = read_biallelic_snps(vcf)
    df = pd.read_csv(sheet)
    samples = list(sites.samples)
    a = [i for i, s in enumerate(samples) if s.startswith(("WOLF", "DOG"))]
    b = [i for i, s in enumerate(samples) if s.startswith("COYOTE")]
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, idx in (("A", a), ("B", b)):
        vpath = out / f"study{name}.vcf"
        write_minimal_vcf(vpath, sites.chrom, sites.pos, sites.ref, sites.alt,
                          sites.samples[idx], sites.gt[:, idx])
        keep = {samples[i] for i in idx}
        spath = out / f"study{name}.csv"
        df[df["sample_id"].isin(keep)].to_csv(spath, index=False)
        paths[name] = (vpath, spath, len(idx))
    return paths


@pytest.fixture(scope="module")
def harmonized(tmp_path_factory) -> GlobalConfig:
    ws = tmp_path_factory.mktemp("harmonize")
    vcf, sheet = simulate_cohort(ws / "raw", seed=7)
    parts = _split_into_studies(vcf, sheet, ws / "studies")
    cfg = GlobalConfig.load(overrides={
        "project_name": "harm", "paths.root": str(ws),
        "pipeline": ["harmonize", "load_genotypes", "pca"], "logging.level": "WARNING",
    })
    cfg = cfg.model_copy(update={"stages": {"harmonize": {"datasets": [
        {"vcf": str(parts["A"][0]), "sample_sheet": str(parts["A"][1]),
         "dataset_id": "studyA"},
        {"vcf": str(parts["B"][0]), "sample_sheet": str(parts["B"][1]),
         "dataset_id": "studyB"},
    ]}}})
    run_pipeline(cfg)
    return cfg


def test_harmonize_merges_datasets(harmonized) -> None:
    art = _store(harmonized).get(ArtifactKind.CALLSET, "callset")
    assert art.metadata["n_datasets"] == 2
    assert art.metadata["n_samples"] == 18            # all samples reunited
    assert art.metadata["n_shared_sites"] > 100       # all sites shared in this split
    # downstream PCA ran on the harmonized callset
    assert _store(harmonized).has(R, "pca")


# -- Phase 5: local ancestry -----------------------------------------------------------


@pytest.fixture(scope="module")
def lai(tmp_path_factory) -> GlobalConfig:
    ws = tmp_path_factory.mktemp("lai")
    vcf, sheet = simulate_introgression_cohort(ws / "in", seed=11, admixture_proportion=0.25)
    cfg = GlobalConfig.load(overrides={
        "project_name": "lai", "paths.root": str(ws),
        "pipeline": ["ingest", "load_genotypes", "local_ancestry"],
        "logging.level": "WARNING",
        "stages.ingest.sample_sheet": str(sheet), "stages.ingest.callset": str(vcf),
        "stages.local_ancestry.sources": ["wolf", "coyote"],
        "stages.local_ancestry.targets": ["wolf", "dog"],
        "stages.local_ancestry.window_bp": 200000,
    })
    run_pipeline(cfg)
    return cfg


def test_local_ancestry_detects_dog_coyote_segments(lai) -> None:
    frac = pd.read_csv(_store(lai).get(R, "local_ancestry").path)
    dog = frac[frac["population"] == "dog"]["frac_coyote"].mean()
    wolf = frac[frac["population"] == "wolf"]["frac_coyote"].mean()
    assert dog > 0.05          # coyote-ancestry segments detected in dogs
    assert dog > wolf          # ... and more than in (non-introgressed) wolves
    assert wolf < 0.15         # wolves are ~pure


# -- Phase 6: selection, demography, ancient -------------------------------------------


@pytest.fixture(scope="module")
def phase6(tmp_path_factory) -> GlobalConfig:
    ws = tmp_path_factory.mktemp("phase6")
    vcf, sheet = simulate_cohort(ws / "in", seed=7)
    cfg = GlobalConfig.load(overrides={
        "project_name": "p6", "paths.root": str(ws),
        "pipeline": ["ingest", "load_genotypes", "selection", "demography",
                     "pseudohaploid"],
        "logging.level": "WARNING",
        "stages.ingest.sample_sheet": str(sheet), "stages.ingest.callset": str(vcf),
        "stages.selection.window_bp": 100000,
        "stages.demography.mutation_rate": 1e-8,
        # The simulated 2 Mb chromosome has a known callable denominator; production
        # reduced panels intentionally leave this unset and suppress whole-genome estimates.
        "stages.demography.callable_sites": 2_000_000,
    })
    run_pipeline(cfg)
    return cfg


def test_selection_produces_windows(phase6) -> None:
    art = _store(phase6).get(R, "selection")
    assert art.metadata["n_windows"] > 0
    assert np.isfinite(art.metadata["max_pbs"])
    df = pd.read_csv(art.path)
    assert {"chrom", "start", "end", "pbs", "tajima_d"} <= set(df.columns)


def test_demography_recovers_ne_ordering(phase6) -> None:
    demo = pd.read_csv(_store(phase6).get(R, "demography").path).set_index("population")
    # simulated Ne: coyote (12k) > wolf (10k) > dog (8k) -> pi tracks it
    assert demo.loc["coyote", "pi_per_site"] > demo.loc["dog", "pi_per_site"]
    assert "Ne" in demo.columns and demo.loc["coyote", "Ne"] > demo.loc["dog", "Ne"]


def test_pseudohaploid_removes_heterozygosity(phase6) -> None:
    art = _store(phase6).get(R, "ancient")
    summary = pd.read_csv(art.path)
    assert (summary["het_after"] == 0).all()          # pseudohaploid => no hets
    assert (summary["het_before"] > 0).any()          # original had hets
    assert _store(phase6).has(ArtifactKind.GENOTYPES, "pseudohaploid")
