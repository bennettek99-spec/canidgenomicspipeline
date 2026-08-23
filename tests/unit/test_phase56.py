from __future__ import annotations

import numpy as np

from canidae.stages.local_ancestry.hmm import genotype_loglik, infer_local_ancestry, viterbi
from canidae.stages.selection.scans import _pbs


def test_genotype_loglik_prefers_matching_frequency() -> None:
    n_alt = np.array([2])  # homozygous ALT
    high = genotype_loglik(n_alt, np.array([0.9]))[0]
    low = genotype_loglik(n_alt, np.array([0.1]))[0]
    assert high > low
    # missing contributes zero
    assert genotype_loglik(np.array([-1]), np.array([0.5]))[0] == 0.0


def test_viterbi_smooths_two_regimes() -> None:
    # first 5 windows favor state 0, last 5 favor state 1
    n = 10
    emissions = np.zeros((n, 2))
    emissions[:5, 0] = 1.0
    emissions[5:, 1] = 1.0
    path = viterbi(emissions, switch_prob=0.02)
    assert list(path) == [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]


def test_viterbi_single_state() -> None:
    assert list(viterbi(np.zeros((4, 1)), 0.02)) == [0, 0, 0, 0]


def test_local_ancestry_restarts_hmm_at_each_chromosome() -> None:
    # Each chromosome has only one weakly informative window. If an HMM transition were
    # allowed across the boundary, the high switch penalty would incorrectly retain "low"
    # on chromosome 2. Independent chromosome runs select the locally preferred source.
    calls = infer_local_ancestry(
        n_alt=np.array([0, 2]),
        chrom=np.array(["1", "2"]),
        pos=np.array([100, 100]),
        source_freqs={
            "low": np.array([0.4, 0.4]),
            "high": np.array([0.6, 0.6]),
        },
        window_bp=1_000,
        switch_prob=0.02,
        min_sites=1,
    )

    assert [(call["chrom"], call["ancestry"]) for call in calls] == [
        ("1", "low"),
        ("2", "high"),
    ]


def test_pbs_formula() -> None:
    # PBS = (T_ab + T_ac - T_bc)/2 with T = -log(1 - Fst)
    fst_ab, fst_ac, fst_bc = 0.3, 0.3, 0.1
    t = lambda f: -np.log(1 - f)  # noqa: E731
    expected = (t(fst_ab) + t(fst_ac) - t(fst_bc)) / 2
    assert abs(_pbs(fst_ab, fst_ac, fst_bc) - expected) < 1e-9
    # a focal population with more branch-specific drift has higher PBS
    assert _pbs(0.5, 0.5, 0.1) > _pbs(0.2, 0.2, 0.1)


def test_pbs_nan_on_missing() -> None:
    assert np.isnan(_pbs(float("nan"), 0.3, 0.1))
