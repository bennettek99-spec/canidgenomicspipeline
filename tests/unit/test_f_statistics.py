"""f4 / D / f4-ratio estimators and the weighted block jackknife."""

from __future__ import annotations

import numpy as np
import pytest

from canidae.analysis.f_statistics import (
    JackknifeEstimate,
    d_statistic,
    f4,
    f4_ratio,
    gl_population_frequencies,
    population_frequencies,
    qpadm,
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


def _phred(likelihoods: np.ndarray) -> np.ndarray:
    scaled = likelihoods / likelihoods.max(axis=-1, keepdims=True)
    return np.minimum(np.round(-10 * np.log10(np.maximum(scaled, 1e-30))), 255).astype(np.uint8)


def test_gl_frequencies_match_confident_calls_and_ignore_empty_pl() -> None:
    confident = {0: (0, 60, 255), 1: (60, 0, 60), 2: (255, 60, 0)}
    pl = np.array(
        [
            [confident[0], confident[1], confident[2], (0, 0, 0)],
            [(0, 0, 0)] * 4,
        ],
        dtype=np.uint8,
    )
    samples = ["a", "b", "c", "d"]
    freqs = gl_population_frequencies(pl, samples, {"p": samples, "empty": ["d"]})
    assert freqs["p"][0] == pytest.approx(0.5, abs=1e-4)
    assert np.isnan(freqs["p"][1]) and np.isnan(freqs["empty"]).all()


def test_gl_frequencies_are_unbiased_at_low_depth_where_hard_calls_are_not() -> None:
    rng = np.random.default_rng(3)
    n_sites, n_ind, depth, err, p_true = 3000, 20, 2, 0.01, 0.3
    genotypes = rng.binomial(2, p_true, size=(n_sites, n_ind))
    alt_reads = np.asarray(rng.binomial(depth, np.clip(genotypes / 2, err, 1 - err)))
    lik = np.stack(
        [
            (err**alt_reads) * ((1 - err) ** (depth - alt_reads)),
            np.full(alt_reads.shape, 0.5**depth),
            ((1 - err) ** alt_reads) * (err ** (depth - alt_reads)),
        ],
        axis=-1,
    )
    pl = _phred(lik)
    samples = [f"s{i}" for i in range(n_ind)]
    gl = gl_population_frequencies(pl, samples, {"p": samples})["p"]
    hard = np.argmin(pl, axis=-1).mean(axis=1) / 2
    assert gl.mean() == pytest.approx(p_true, abs=0.01)
    assert abs(hard.mean() - p_true) > 3 * abs(gl.mean() - p_true)


def _simulate_three_source(
    weights: tuple[float, float, float], seed: int
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """X = wW*W + wD*D + wC*C on the tree (O, ((C, C2), (A, (W, (D, D2)))))."""
    msprime = pytest.importorskip("msprime")
    demography = msprime.Demography()
    names = ("O", "C", "C2", "A", "W", "D", "D2", "X", "XW", "XD", "XC")
    for name in (*names, "DD", "WD", "AW", "CC", "ROOTI", "ROOT"):
        demography.add_population(name=name, initial_size=5_000)
    w_w, w_d, w_c = weights
    demography.add_admixture(
        time=40, derived="X", ancestral=["XW", "XD", "XC"], proportions=[w_w, w_d, w_c]
    )
    demography.add_population_split(time=500, derived=["XW"], ancestral="W")
    demography.add_population_split(time=500, derived=["XD"], ancestral="D")
    demography.add_population_split(time=500, derived=["XC"], ancestral="C")
    demography.add_population_split(time=1_500, derived=["D", "D2"], ancestral="DD")
    demography.add_population_split(time=1_500, derived=["C", "C2"], ancestral="CC")
    demography.add_population_split(time=3_000, derived=["W", "DD"], ancestral="WD")
    demography.add_population_split(time=5_000, derived=["A", "WD"], ancestral="AW")
    demography.add_population_split(time=10_000, derived=["AW", "CC"], ancestral="ROOTI")
    demography.add_population_split(time=25_000, derived=["ROOTI", "O"], ancestral="ROOT")
    demography.sort_events()
    sampled = ("O", "C", "C2", "A", "W", "D", "D2", "X")
    freqs: dict[str, list[np.ndarray]] = {name: [] for name in sampled}
    blocks: list[np.ndarray] = []
    replicates = msprime.sim_ancestry(
        samples={name: 10 for name in sampled},
        demography=demography,
        sequence_length=1_000_000,
        recombination_rate=1e-8,
        random_seed=seed,
        num_replicates=30,
    )
    pop_ids = {pop.name: pop.id for pop in demography.populations}
    for block, ts in enumerate(replicates):
        mts = msprime.sim_mutations(ts, rate=2e-8, random_seed=seed + block + 1)
        genotypes = mts.genotype_matrix()
        for name in sampled:
            nodes = [n for n in mts.samples() if mts.node(n).population == pop_ids[name]]
            freqs[name].append(genotypes[:, nodes].mean(axis=1))
        blocks.append(np.full(genotypes.shape[0], block))
    return {k: np.concatenate(v) for k, v in freqs.items()}, np.concatenate(blocks)


def test_qpadm_recovers_three_way_mixture_and_rejects_a_missing_source() -> None:
    truth = (0.5, 0.2, 0.3)
    freqs, blocks = _simulate_three_source(truth, seed=21)
    rights = ["O", "A", "D2", "C2"]
    fit = qpadm(freqs, "X", ["W", "D", "C"], rights, blocks)
    assert fit.weights.sum() == pytest.approx(1.0)
    for estimate, se, expected in zip(fit.weights, fit.se, truth, strict=True):
        assert abs(estimate - expected) < max(4 * se, 0.05)
    assert fit.dof == 1 and fit.p_value > 0.01
    # Dropping the dog source leaves allele sharing with D2 unexplained.
    wrong = qpadm(freqs, "X", ["W", "C"], rights, blocks)
    assert wrong.p_value < 0.01


def test_qpadm_rejects_unidentifiable_designs() -> None:
    ones = np.ones(20)
    freqs = {"x": ones, "s1": ones, "s2": ones, "r0": ones, "r1": ones}
    blocks: np.ndarray = np.repeat(np.arange(4), 5)
    with pytest.raises(ValueError, match="at least as many rights"):
        qpadm(freqs, "x", ["s1", "s2"], ["r0"], blocks)
    with pytest.raises(ValueError, match="rank-deficient"):
        qpadm(freqs, "x", ["s1", "s2"], ["r0", "r1"], blocks)
