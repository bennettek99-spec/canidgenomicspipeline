"""f4 / D / f4-ratio estimators and the weighted block jackknife."""

from __future__ import annotations

import numpy as np
import pytest

from canidae.analysis.f_statistics import (
    JackknifeEstimate,
    d_statistic,
    f4,
    f4_ratio,
    population_frequencies,
    weighted_block_jackknife,
)


def test_population_frequencies_ignore_missing_calls() -> None:
    gt = np.array([[0, 2, -1], [1, -1, -1], [-1, -1, -1]], dtype=np.int8)
    freqs = population_frequencies(gt, ["s1", "s2", "s3"], {"p": ["s1", "s2"], "q": ["s3"]})
    np.testing.assert_allclose(freqs["p"][:2], [0.5, 0.5])
    assert np.isnan(freqs["p"][2])
    assert np.isnan(freqs["q"]).all()


def test_jackknife_constant_ratio_has_zero_se() -> None:
    num = np.arange(1.0, 101.0)
    est = weighted_block_jackknife(num, 2 * num, np.repeat(np.arange(10), 10))
    assert est.estimate == pytest.approx(0.5)
    assert est.se == pytest.approx(0.0, abs=1e-12)
    assert (est.n_sites, est.n_blocks) == (100, 10)


def test_jackknife_equal_blocks_matches_textbook_formula() -> None:
    rng = np.random.default_rng(1)
    num = rng.normal(size=200)
    den = np.ones(200)
    blocks: np.ndarray = np.repeat(np.arange(20), 10)
    est = weighted_block_jackknife(num, den, blocks)
    # With equal blocks and a mean, the jackknife SE is the SE of the block means.
    block_means = num.reshape(20, 10).mean(axis=1)
    assert est.estimate == pytest.approx(num.mean())
    assert est.se == pytest.approx(block_means.std(ddof=1) / np.sqrt(20))


def test_jackknife_degenerate_inputs() -> None:
    nan = weighted_block_jackknife(np.array([np.nan]), np.array([1.0]), np.array([0]))
    assert nan.n_sites == 0 and np.isnan(nan.estimate)
    one_block = weighted_block_jackknife(np.ones(5), np.ones(5), np.zeros(5))
    assert one_block.estimate == 1.0 and np.isnan(one_block.se)
    assert np.isnan(JackknifeEstimate(1.0, 0.0, 1, 1).z)


def test_f4_and_d_signs_follow_shared_drift() -> None:
    # P1 and P3 share the derived allele at every site; P2 and P4 do not.
    ones, zeros = np.ones(40), np.zeros(40)
    freqs = {"p1": ones, "p2": zeros, "p3": ones, "p4": zeros}
    blocks: np.ndarray = np.repeat(np.arange(4), 10)
    assert f4(freqs, "p1", "p2", "p3", "p4", blocks).estimate == pytest.approx(1.0)
    assert d_statistic(freqs, "p1", "p2", "p3", "p4", blocks).estimate == pytest.approx(1.0)
    assert d_statistic(freqs, "p2", "p1", "p3", "p4", blocks).estimate == pytest.approx(-1.0)


def _simulate_admixed_frequencies(
    alpha: float, seed: int
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    msprime = pytest.importorskip("msprime")
    demography = msprime.Demography()
    for name in ("O", "C", "A", "B", "X", "XB", "XC", "BC", "AB", "ROOT"):
        demography.add_population(name=name, initial_size=5_000)
    demography.add_admixture(
        time=50, derived="X", ancestral=["XB", "XC"], proportions=[alpha, 1 - alpha]
    )
    demography.add_population_split(time=1_000, derived=["XB"], ancestral="B")
    demography.add_population_split(time=1_000, derived=["XC"], ancestral="C")
    demography.add_population_split(time=3_000, derived=["A", "B"], ancestral="AB")
    demography.add_population_split(time=8_000, derived=["AB", "C"], ancestral="BC")
    demography.add_population_split(time=20_000, derived=["BC", "O"], ancestral="ROOT")
    demography.sort_events()
    samples = {name: 10 for name in ("O", "C", "A", "B", "X")}
    freqs: dict[str, list[np.ndarray]] = {name: [] for name in samples}
    blocks: list[np.ndarray] = []
    replicates = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=1_000_000,
        recombination_rate=1e-8,
        random_seed=seed,
        num_replicates=30,
    )
    for block, ts in enumerate(replicates):
        mts = msprime.sim_mutations(ts, rate=2e-8, random_seed=seed + block + 1)
        genotypes = mts.genotype_matrix()
        for pop_id, name in enumerate(samples):
            nodes = [n for n in mts.samples() if mts.node(n).population == pop_id]
            freqs[name].append(genotypes[:, nodes].mean(axis=1))
        blocks.append(np.full(genotypes.shape[0], block))
    return {k: np.concatenate(v) for k, v in freqs.items()}, np.concatenate(blocks)


@pytest.mark.parametrize("alpha", [0.25, 0.8])
def test_f4_ratio_recovers_simulated_admixture(alpha: float) -> None:
    freqs, blocks = _simulate_admixed_frequencies(alpha, seed=11)
    est = f4_ratio(freqs, a="A", o="O", x="X", b="B", c="C", blocks=blocks)
    assert est.n_blocks == 30
    assert abs(est.estimate - alpha) < max(4 * est.se, 0.06)
    # Unadmixed controls sit at the ends of the scale.
    assert f4_ratio(freqs, "A", "O", "C", "B", "C", blocks).estimate == pytest.approx(0.0)
    assert f4_ratio(freqs, "A", "O", "B", "B", "C", blocks).estimate == pytest.approx(1.0)
