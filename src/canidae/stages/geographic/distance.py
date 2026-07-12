"""Geographic distance and the Mantel test for isolation-by-distance."""

from __future__ import annotations

import numpy as np
from scipy.stats import pearsonr

_EARTH_RADIUS_KM = 6371.0088


def haversine_matrix(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Great-circle distance matrix (km) for arrays of latitudes/longitudes (degrees)."""
    lat_r = np.radians(np.asarray(lat, dtype=float))
    lon_r = np.radians(np.asarray(lon, dtype=float))
    dlat = lat_r[:, None] - lat_r[None, :]
    dlon = lon_r[:, None] - lon_r[None, :]
    a = np.sin(dlat / 2) ** 2 + (
        np.cos(lat_r[:, None]) * np.cos(lat_r[None, :]) * np.sin(dlon / 2) ** 2
    )
    return 2 * _EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def mantel_test(
    d1: np.ndarray, d2: np.ndarray, *, permutations: int = 999, seed: int = 0
) -> tuple[float, float]:
    """Mantel test: Pearson correlation of two distance matrices with a permutation p-value.

    Returns ``(r, p)``. ``p`` is the fraction of permutations whose \\|r\\| is at least the
    observed \\|r\\| (two-sided), with the +1 correction.
    """
    n = d1.shape[0]
    iu = np.triu_indices(n, k=1)
    x = d1[iu]
    y = d2[iu]
    if x.size < 2 or np.all(x == x[0]) or np.all(y == y[0]):
        return float("nan"), float("nan")
    r_obs = float(pearsonr(x, y)[0])

    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(permutations):
        perm = rng.permutation(n)
        y_perm = d2[np.ix_(perm, perm)][iu]
        if abs(float(pearsonr(x, y_perm)[0])) >= abs(r_obs):
            count += 1
    p = (count + 1) / (permutations + 1)
    return r_obs, p
