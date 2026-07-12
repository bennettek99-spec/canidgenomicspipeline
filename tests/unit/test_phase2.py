from __future__ import annotations

from pathlib import Path

import allel
import numpy as np

from canidae.stages.geographic.distance import haversine_matrix, mantel_test
from canidae.stages.popgen.ancestry_nmf import (
    cross_validate_k,
    fit_admixture,
    select_k,
    weighted_nmf,
)
from canidae.stages.popgen.plink import decode_bed, write_plink_bed
from canidae.stages.popgen.store import Genotypes

# -- ancestry NMF ----------------------------------------------------------------------


def _two_group_matrix() -> np.ndarray:
    rng = np.random.default_rng(0)
    X = rng.random((10, 100)) * 0.2
    X[:5, :50] += 2.0     # group 1 carries alt alleles in the first half
    X[5:, 50:] += 2.0     # group 2 in the second half
    return X


def test_weighted_nmf_reconstructs() -> None:
    X = _two_group_matrix()
    mask = np.ones_like(X)
    W, H = weighted_nmf(X, mask, 2, seed=0)
    err = np.sqrt(np.mean((X - W @ H) ** 2))
    assert err < 0.3


def test_fit_admixture_separates_groups() -> None:
    fit = fit_admixture(_two_group_matrix(), 2, seed=0)
    assert np.allclose(fit.Q.sum(axis=1), 1.0, atol=1e-6)
    dom = fit.Q.argmax(axis=1)
    assert len(set(dom[:5])) == 1 and len(set(dom[5:])) == 1
    assert dom[0] != dom[5]


def test_cross_validate_and_select_k() -> None:
    errs = cross_validate_k(_two_group_matrix(), [1, 2, 3], seed=0)
    assert set(errs) == {1, 2, 3}
    assert select_k(errs) in {1, 2, 3}
    assert errs[2] < errs[1]  # 2 components beat 1 for two-group data


# -- PLINK bed round-trip --------------------------------------------------------------


def test_plink_bed_roundtrip(tmp_path: Path) -> None:
    gt = np.array(
        [
            [[0, 0], [0, 1], [1, 1], [-1, -1]],
            [[1, 1], [0, 0], [0, 1], [0, 0]],
        ],
        dtype="i1",
    )
    geno = Genotypes(
        calls=allel.GenotypeArray(gt),
        pos=np.array([100, 200], dtype=np.int64),
        chrom=np.array(["1", "1"]),
        samples=np.array(["A", "B", "C", "D"]),
    )
    fileset = write_plink_bed(geno, tmp_path / "cohort")
    assert fileset.bed.exists() and fileset.bim.exists() and fileset.fam.exists()
    decoded = decode_bed(fileset, n_samples=4)
    expected = np.asarray(geno.calls.to_n_alt(fill=-1))
    np.testing.assert_array_equal(decoded, expected)


# -- geographic distance + Mantel ------------------------------------------------------


def test_haversine_known_distance() -> None:
    d = haversine_matrix(np.array([0.0, 0.0]), np.array([0.0, 90.0]))
    assert d[0, 0] == 0.0
    assert abs(d[0, 1] - 10007.5) < 50.0  # quarter of Earth's circumference along equator
    assert np.allclose(d, d.T)


def test_mantel_identical_matrices() -> None:
    rng = np.random.default_rng(1)
    pts = rng.random((8, 2)) * 100
    from scipy.spatial.distance import pdist, squareform

    d = squareform(pdist(pts))
    r, p = mantel_test(d, d, permutations=199, seed=0)
    assert r > 0.99
    assert p <= 0.05
