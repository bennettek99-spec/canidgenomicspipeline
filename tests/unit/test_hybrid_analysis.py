"""Unit tests for hybrid-canid analysis libraries and stages.

The numerical tests pin the mixture math against hand-computed values rather
than against whatever the code currently returns, so a scale or convention
change (the class of bug that silently doubled the multiway wolf fractions)
fails here instead of drifting into published results.
"""

from __future__ import annotations

import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

from canidae.analysis.admixture_ml import (
    best_mixture,
    bootstrap_ci,
    calls_mixture_loglik,
    cohort_grid,
    panel_frequencies,
    two_way_fit,
    wgs_allele_count_matrix,
)
from canidae.analysis.breed_panel import breed_of, classify_group
from canidae.analysis.bridge_panel import (
    BRIDGE_WGS_IDS,
    allele_count_vector,
    bridge_sites,
    load_bridge_genotypes,
)
from canidae.analysis.genotypes import allele_count, dosage, parse_gt_field
from canidae.analysis.reference_mixture import (
    allele_frequencies,
    binom2_logpmf,
    infer_dog_fraction,
    leave_one_out,
)
from canidae.core.config import GlobalConfig
from canidae.core.registry import STAGES
from canidae.pipeline import instantiate_stages
from canidae.stages import load_builtin_stages
from canidae.stages.hybrid.io_utils import load_calls_csv

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulated"))
from make_hybrid_panel import BUG_COLUMN, build_hybrid_panel

# -- genotype scale conventions --------------------------------------------------------


def test_allele_count_and_dosage_scales() -> None:
    """``allele_count`` is 0/1/2; ``dosage`` is the same value rescaled to [0, 1]."""
    assert [allele_count(gt) for gt in ("0/0", "0/1", "1|1")] == [0, 1, 2]
    assert [dosage(gt) for gt in ("0/0", "0/1", "1|1")] == [0.0, 0.5, 1.0]
    assert parse_gt_field("0/1:30", 0) == "0/1"
    assert parse_gt_field("./.:30", 0) is None
    assert parse_gt_field("0/1", 3) is None


def test_binom2_logpmf_matches_the_binomial_pmf() -> None:
    """Each of the three genotypes must hit its own branch, hom-ALT included."""
    p = 0.9
    assert binom2_logpmf(0, p) == pytest.approx(2 * math.log(1 - p))
    assert binom2_logpmf(1, p) == pytest.approx(math.log(2 * p * (1 - p)))
    assert binom2_logpmf(2, p) == pytest.approx(2 * math.log(p))
    # Hom-ALT and hom-REF must not collapse onto the het branch.
    assert binom2_logpmf(2, p) > binom2_logpmf(1, p) > binom2_logpmf(0, p)
    # Clamping keeps the degenerate panels finite.
    assert math.isfinite(binom2_logpmf(0, 0.0))
    assert math.isfinite(binom2_logpmf(2, 1.0))


def test_allele_frequencies_are_laplace_smoothed_on_the_count_scale() -> None:
    """A fixed-ALT group must approach 1.0, not the 0.5 a dosage-scale sum gives."""
    keys = [("chr1", 10)]
    all_alt = {("chr1", 10): ["1/1"] * 10}
    all_ref = {("chr1", 10): ["0/0"] * 10}
    all_het = {("chr1", 10): ["0/1"] * 10}
    assert allele_frequencies(keys, all_alt, list(range(10)))[0] == pytest.approx(21 / 22)
    assert allele_frequencies(keys, all_ref, list(range(10)))[0] == pytest.approx(1 / 22)
    assert allele_frequencies(keys, all_het, list(range(10)))[0] == pytest.approx(0.5)


def test_panel_frequencies_shrink_toward_the_pooled_prior() -> None:
    # Three samples, hom-REF at locus 0 and hom-ALT at locus 1.
    matrix = np.array([[0.0, 2.0], [0.0, 2.0], [0.0, 2.0]], dtype=float)
    prior = np.array([0.5, 0.5])
    freqs = panel_frequencies(matrix, [0, 1, 2], prior)
    assert freqs.shape == (2,)
    # (0 + 6*0.5) / (2*3 + 6) = 0.25 and (6 + 6*0.5) / 12 = 0.75.
    assert freqs[0] == pytest.approx(0.25)
    assert freqs[1] == pytest.approx(0.75)
    # A panel with no calls at a locus yields NaN rather than the prior.
    assert np.isnan(panel_frequencies(np.full((2, 1), np.nan), [0, 1], np.array([0.5]))[0])


def test_wgs_allele_count_matrix_shape_and_missingness() -> None:
    keys = [("chr1", 10), ("chr1", 20)]
    records = {("chr1", 10): ["0/0", "1/1"], ("chr1", 20): ["0/1", None]}
    matrix = wgs_allele_count_matrix(keys, records)
    assert matrix.shape == (2, 2)
    assert matrix[0].tolist() == [0.0, 1.0]
    assert matrix[1][0] == 2.0
    assert np.isnan(matrix[1][1])


# -- breed taxonomy --------------------------------------------------------------------


def test_classify_group_edge_cases() -> None:
    assert classify_group("Beagle") == "breed"
    # A domestic breed whose name contains a wild-canid substring.
    assert classify_group("IrishWolfhound") == "breed"
    assert classify_group("Coyote") == "wild"
    assert classify_group("AlaskanWolf") == "wild"
    assert classify_group("GoldenJackal") == "wild"
    assert classify_group("VillDog_Peru") == "village"
    assert classify_group("PER") == "village"  # inferred village code
    assert classify_group("CFA.") == "ambiguous"  # undocumented code
    assert classify_group("BC") == "ambiguous"  # too short to map safely
    assert classify_group("Helsinki_BC") == "ambiguous"
    assert classify_group("140447_S") == "ambiguous"  # contains digits
    assert classify_group("UnknownMix") == "mixed"
    assert breed_of("Beagle01") == "Beagle"
    assert breed_of("AlgonquinWolf13467") == "AlgonquinWolf"
    assert breed_of("VillDog_Peru02") == "VillDog_Peru"


# -- two-source mixture ----------------------------------------------------------------


def test_infer_dog_fraction_recovers_a_known_f1() -> None:
    wgs = ["Coy1", "Coy2", "Dog1", "Dog2"]
    coyote = frozenset({"Coy1", "Coy2"})
    # One diagnostic locus: coyotes ref, dogs alt; a het query is a clean F1.
    reference = {("chr1", 100): ["0/0", "0/0", "1/1", "1/1"]}
    calls = {("chr1", 100): ("0/1", 20, 10, 10)}
    est, called, rmse = infer_dog_fraction(
        calls, reference, wgs_samples=wgs, coyote_samples=coyote, min_called=1
    )
    assert called == 1
    assert est == pytest.approx(0.5)
    assert rmse == pytest.approx(0.0, abs=1e-9)


def test_infer_dog_fraction_returns_none_below_min_called() -> None:
    wgs = ["Coy1", "Dog1"]
    reference = {("chr1", 100): ["0/0", "1/1"]}
    calls = {("chr1", 100): ("0/1", 20, 10, 10)}
    est, called, rmse = infer_dog_fraction(
        calls, reference, wgs_samples=wgs, coyote_samples=frozenset({"Coy1"}),
        min_called=20,
    )
    assert est is None and rmse is None and called == 1


def test_infer_dog_fraction_skips_non_diagnostic_loci() -> None:
    """Loci whose coyote/dog means barely differ must not enter the fit."""
    wgs = ["Coy1", "Dog1"]
    reference = {
        ("chr1", 100): ["0/0", "1/1"],  # difference 1.0, diagnostic
        ("chr1", 200): ["0/1", "0/1"],  # difference 0.0, uninformative
    }
    calls = {
        ("chr1", 100): ("1/1", 20, 0, 20),
        ("chr1", 200): ("0/0", 20, 20, 0),
    }
    est, called, _ = infer_dog_fraction(
        calls, reference, wgs_samples=wgs, coyote_samples=frozenset({"Coy1"}),
        min_called=1,
    )
    assert called == 1
    assert est == pytest.approx(1.0)


def test_leave_one_out_separates_two_clean_panels() -> None:
    keys = [("chr1", 10), ("chr1", 20)]
    samples = ["A1", "A2", "A3", "B1", "B2", "B3"]
    records = {
        ("chr1", 10): ["0/0", "0/0", "0/0", "1/1", "1/1", "1/1"],
        ("chr1", 20): ["0/0", "0/0", "0/0", "1/1", "1/1", "1/1"],
    }
    groups = {"A": [0, 1, 2], "B": [3, 4, 5]}
    assert allele_frequencies(keys, records, groups["A"])[0] < 0.5
    assert allele_frequencies(keys, records, groups["B"])[0] > 0.5
    rows, summary = leave_one_out(keys, records, samples, groups, min_group=2)
    assert summary["n_tested"] == 6
    assert summary["top1_accuracy"] == 1.0
    assert all(row["correct"] for row in rows)
    assert {row["sample_id"] for row in rows} == set(samples)


def test_leave_one_out_skips_groups_below_min_group() -> None:
    keys = [("chr1", 10)]
    records = {("chr1", 10): ["0/0", "0/0", "1/1"]}
    rows, summary = leave_one_out(
        keys, records, ["A1", "A2", "B1"], {"A": [0, 1], "B": [2]}, min_group=2
    )
    assert {row["true_breed"] for row in rows} == {"A"}
    assert summary["n_tested"] == 2


# -- three-way mixture -----------------------------------------------------------------


def test_best_mixture_recovers_a_pure_wolf_sample() -> None:
    # Locus 0 separates wolf, locus 1 separates dog, locus 2 is uninformative.
    p_c = np.array([0.02, 0.02, 0.5])
    p_w = np.array([0.98, 0.02, 0.5])
    p_d = np.array([0.02, 0.98, 0.5])
    # Hom-ALT at the wolf locus, hom-REF at the dog locus.
    f_wolf, f_dog, ll = best_mixture(np.array([2.0, 0.0, 1.0]), p_c, p_w, p_d)
    assert f_wolf > 0.9
    assert f_dog < 0.1
    assert math.isfinite(ll)


def test_best_mixture_ignores_uncalled_loci() -> None:
    p_c, p_w, p_d = (np.array([0.02, 0.02]), np.array([0.98, 0.98]),
                     np.array([0.02, 0.02]))
    both = best_mixture(np.array([2.0, np.nan]), p_c, p_w, p_d)
    only = best_mixture(np.array([2.0]), p_c[:1], p_w[:1], p_d[:1])
    assert both[:2] == only[:2]
    assert both[2] == pytest.approx(only[2])


def test_calls_mixture_loglik_matches_the_hand_computed_value() -> None:
    p_c, p_w, p_d = np.array([0.1]), np.array([0.5]), np.array([0.9])
    f_wolf, f_dog = 0.2, 0.3
    p_eff = (1 - f_wolf - f_dog) * 0.1 + f_wolf * 0.5 + f_dog * 0.9
    for count, expected in (
        (0.0, 2 * math.log1p(-p_eff)),
        (1.0, math.log(2 * p_eff * (1 - p_eff))),
        (2.0, 2 * math.log(p_eff)),
    ):
        got = calls_mixture_loglik(np.array([count]), p_c, p_w, p_d, f_wolf, f_dog)
        assert got == pytest.approx(expected)


def test_cohort_grid_recovers_a_known_pooled_mixture() -> None:
    rng = np.random.default_rng(11)
    n_loci = 300
    p_c = rng.uniform(0.02, 0.2, n_loci)
    p_w = rng.uniform(0.8, 0.98, n_loci)
    p_d = rng.uniform(0.02, 0.2, n_loci)
    f_wolf_true = 0.30
    p_eff = (1 - f_wolf_true) * p_c + f_wolf_true * p_w
    chromosomes = np.full(n_loci, 40)
    alt = rng.binomial(chromosomes, p_eff)
    f_wolf, f_dog, ll = cohort_grid(alt, chromosomes, p_c, p_w, p_d)
    assert f_wolf == pytest.approx(f_wolf_true, abs=0.06)
    assert f_dog < 0.15
    assert math.isfinite(ll)


def test_two_way_fit_recovers_a_known_fraction() -> None:
    rng = np.random.default_rng(3)
    n_loci = 300
    p_from = rng.uniform(0.02, 0.2, n_loci)
    p_to = rng.uniform(0.8, 0.98, n_loci)
    f_true = 0.4
    chromosomes = np.full(n_loci, 40)
    alt = rng.binomial(chromosomes, (1 - f_true) * p_from + f_true * p_to)
    f, ll = two_way_fit(alt, chromosomes, p_from, p_to)
    assert f == pytest.approx(f_true, abs=0.05)
    assert math.isfinite(ll)


def test_bootstrap_ci_shape_and_bracketing() -> None:
    rng = np.random.default_rng(5)
    n_loci = 200
    p_c = np.full(n_loci, 0.05)
    p_w = np.full(n_loci, 0.95)
    p_d = np.full(n_loci, 0.05)
    counts = np.array([1.0] * n_loci)  # every locus het -> ~50% wolf
    ci = bootstrap_ci(counts, p_c, p_w, p_d, rng, bootstrap_n=40)
    assert set(ci) == {"wolf", "dog"}
    for bounds in ci.values():
        assert len(bounds) == 2
        low, high = bounds
        assert 0.0 <= low <= high <= 1.0
    point = best_mixture(counts, p_c, p_w, p_d)[0]
    assert ci["wolf"][0] <= point <= ci["wolf"][1]


def test_bootstrap_ci_is_deterministic_for_a_seeded_generator() -> None:
    p_c, p_w, p_d = (np.full(50, 0.05), np.full(50, 0.95), np.full(50, 0.05))
    counts = np.array([1.0] * 50)
    first = bootstrap_ci(counts, p_c, p_w, p_d, np.random.default_rng(7), bootstrap_n=20)
    second = bootstrap_ci(counts, p_c, p_w, p_d, np.random.default_rng(7), bootstrap_n=20)
    assert first == second


# -- bridge panel I/O ------------------------------------------------------------------


def test_bridge_sites_rejects_an_undersized_panel(tmp_path: Path) -> None:
    path = tmp_path / "tiny.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.write("chr1\t100\t.\tA\tG\t.\tPASS\t.\n")
    with pytest.raises(RuntimeError, match="unexpectedly small"):
        bridge_sites(path)


def test_bridge_panel_round_trip(tmp_path: Path) -> None:
    panel = build_hybrid_panel(tmp_path / "panel")
    sites = bridge_sites(panel.bridge_vcf)
    assert len(sites) >= 100
    assert sites[("chr1", 1000)] == ("A", "G")

    samples, calls = load_bridge_genotypes(panel.bridge_vcf)
    assert BUG_COLUMN in samples
    counts = allele_count_vector([("chr1", 1000)], calls, samples[0])
    assert counts.shape == (1,)
    assert counts[0] in (0.0, 1.0, 2.0)


def test_bridge_wgs_bug_ids_are_the_documented_mislabeled_columns() -> None:
    """The exclusion list must stay pinned; silently emptying it re-opens the bug."""
    assert frozenset(
        {"Coyote01", "Coyote02", "AlaskanWolf", "AlgonquinWolf13467",
         "AlgonquinWolf13470", "GoldenJackal01"}
    ) == BRIDGE_WGS_IDS


def test_load_calls_csv_masks_shallow_sites(tmp_path: Path) -> None:
    path = tmp_path / "Q1_calls.csv"
    path.write_text(
        "chrom,position,ref,alt,gt,depth,ref_count,alt_count\n"
        "chr1,100,A,G,0/1,20,10,10\n"
        "chr1,200,A,G,1/1,4,0,4\n"
        "chr1,300,A,G,1/1,20,1,19\n",
        encoding="utf-8",
    )
    sites = {("chr1", 100): ("A", "G"), ("chr1", 200): ("A", "G"),
             ("chr1", 300): ("A", "G")}
    calls = load_calls_csv(path, sites, min_depth=8)
    assert calls[("chr1", 100)][0] == "0/1"
    assert calls[("chr1", 200)][0] is None  # below min_depth
    assert calls[("chr1", 300)][0] == "1/1"  # 19/20 >= 0.85


def test_load_calls_csv_requires_full_site_coverage(tmp_path: Path) -> None:
    from canidae.core.errors import StageInputError

    path = tmp_path / "Q1_calls.csv"
    path.write_text(
        "chrom,position,ref,alt,gt,depth,ref_count,alt_count\n"
        "chr1,100,A,G,0/1,20,10,10\n",
        encoding="utf-8",
    )
    sites = {("chr1", 100): ("A", "G"), ("chr1", 999): ("A", "G")}
    with pytest.raises(StageInputError, match="do not cover"):
        load_calls_csv(path, sites, min_depth=8)


# -- registration and presets ----------------------------------------------------------


def test_hybrid_stages_register() -> None:
    load_builtin_stages()
    for name in ("reference_mixture", "breed_assign", "multiway_admixture"):
        assert name in STAGES


def test_nyc_preset_instantiates(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[2]
    cfg = GlobalConfig.load(
        root / "configs/examples/nyc_coydog_validation.yaml",
        overrides={"paths.root": str(tmp_path)},
    )
    stages = instantiate_stages(cfg)
    assert [s.name for s in stages] == ["reference_mixture", "breed_assign", "report"]


def test_eastern_preset_instantiates(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[2]
    cfg = GlobalConfig.load(
        root / "configs/examples/eastern_coyote_ancestry.yaml",
        overrides={"paths.root": str(tmp_path)},
    )
    stages = instantiate_stages(cfg)
    assert [s.name for s in stages] == ["multiway_admixture", "report"]


def test_reference_mixture_stage_recovers_fixture_truth(tmp_path: Path) -> None:
    """End-to-end on the synthetic panel: estimates must land on the known truth."""
    import pandas as pd

    from canidae.core.datastore import DataStore
    from canidae.core.executor import build_context
    from canidae.core.provenance import ProvenanceWriter
    from canidae.core.runtime import LocalRunner
    from canidae.stages.hybrid.reference_mixture import (
        ReferenceMixtureConfig,
        ReferenceMixtureStage,
    )

    panel = build_hybrid_panel(tmp_path / "panel")
    queries = panel.query_ids()
    cfg = GlobalConfig.load(
        include_defaults=True,
        overrides={
            "project_name": "hybrid-fixture",
            "paths.root": str(tmp_path),
            "pipeline": ["reference_mixture"],
        },
    )
    store = DataStore(tmp_path / "store")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    ctx = build_context(
        cfg, store, LocalRunner(),
        ProvenanceWriter(run_dir, config_digest=cfg.digest(), seed=1), run_dir=run_dir,
    )
    result = ReferenceMixtureStage(
        ReferenceMixtureConfig(
            bridge_vcf=panel.bridge_vcf,
            reference_genotypes=panel.reference_genotypes,
            calls_dir=panel.calls_dir,
            query_samples=queries,
            wgs_samples=panel.reference_samples,
            coyote_samples=panel.coyote_samples,
        )
    ).run(ctx)

    assert result.metrics["n_queries"] == len(queries)
    table = pd.read_csv(result.artifacts[0].path).set_index("sample_id")
    for sample_id in panel.two_source_query_ids():
        expected = panel.truth[sample_id].f_dog
        assert table.loc[sample_id, "dog_fraction"] == pytest.approx(expected, abs=0.06)

    # The model has no wolf axis, so a wolf component is absorbed into the dog
    # fraction. That is a documented limitation of the two-source recipe, and
    # the multiway stage is the one that separates the two.
    wolf_carrying = panel.wolf_carrying_query_ids()
    inflation = [
        table.loc[s, "dog_fraction"] - panel.truth[s].f_dog for s in wolf_carrying
    ]
    assert min(inflation) > 0.0

    manifest = json.loads(Path(result.artifacts[0].metadata["manifest"]).read_text())
    assert manifest["validation"]["passed"] is None  # no pedigree criteria configured
    assert manifest["limitations"]
