"""sNMF-style admixture estimation by weighted non-negative matrix factorization.

Model-based ancestry estimators (ADMIXTURE, sNMF) factor the ALT allele-count matrix
``X`` (samples x sites, entries in {0,1,2}) into ancestry proportions ``Q`` (samples x K)
and component ALT frequencies ``F`` (K x sites). This module provides a dependency-light,
fully in-Python estimator using weighted multiplicative-update NMF, plus a Wold-style
entrywise cross-validation to choose K. It is the default backend so admixture runs
anywhere; the ADMIXTURE binary is available as an alternative backend where installed.

Reference: Frichot et al. (2014), "Fast and efficient estimation of individual ancestry
coefficients" — sNMF is sparse NMF applied to exactly this factorization.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_EPS = 1e-9


@dataclass(slots=True)
class AdmixtureFit:
    """Result of an admixture factorization at a single K."""

    K: int
    Q: np.ndarray  # (n_samples, K), rows sum to 1
    F: np.ndarray  # (K, n_sites)
    reconstruction_error: float


def weighted_nmf(
    X: np.ndarray,
    mask: np.ndarray,
    K: int,
    *,
    n_iter: int = 250,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Weighted NMF via multiplicative updates: minimize ||mask * (X - W @ H)||_F.

    ``mask`` is 1 for observed entries, 0 for held-out/missing ones, so the same routine
    serves both fitting (all-ones mask) and cross-validation (held-out entries zeroed).
    """
    rng = np.random.default_rng(seed)
    n, m = X.shape
    W = rng.random((n, K)) + 1e-3
    H = rng.random((K, m)) + 1e-3
    MX = mask * X
    for _ in range(n_iter):
        WH = W @ H
        W *= (MX @ H.T) / ((mask * WH) @ H.T + _EPS)
        WH = W @ H
        H *= (W.T @ MX) / (W.T @ (mask * WH) + _EPS)
    return W, H


def _normalize_rows(W: np.ndarray) -> np.ndarray:
    totals = W.sum(axis=1, keepdims=True)
    return W / np.where(totals > 0, totals, 1.0)


def fit_admixture(X: np.ndarray, K: int, *, seed: int = 0, n_iter: int = 250) -> AdmixtureFit:
    """Fit ancestry proportions Q and component ALT frequencies F at a given K."""
    mask = np.ones_like(X, dtype=float)
    W, H = weighted_nmf(X, mask, K, n_iter=n_iter, seed=seed)
    err = float(np.sqrt(np.mean((X - W @ H) ** 2)))
    return AdmixtureFit(K=K, Q=_normalize_rows(W), F=H, reconstruction_error=err)


def cross_validate_k(
    X: np.ndarray,
    k_values: list[int],
    *,
    holdout: float = 0.1,
    n_iter: int = 200,
    seed: int = 0,
) -> dict[int, float]:
    """Entrywise (Wold) cross-validation error per K: hold out a random fraction of
    entries, fit on the rest, and measure reconstruction error on the held-out entries.

    Returns ``{K: mean_squared_holdout_error}``. The K with the lowest value is preferred.
    """
    rng = np.random.default_rng(seed)
    observed = rng.random(X.shape) >= holdout
    mask = observed.astype(float)
    test = (~observed).astype(float)
    n_test = max(float(test.sum()), 1.0)
    errors: dict[int, float] = {}
    for K in k_values:
        W, H = weighted_nmf(X, mask, K, n_iter=n_iter, seed=seed + K)
        pred = W @ H
        errors[K] = float((test * (X - pred) ** 2).sum() / n_test)
    return errors


def select_k(cv_errors: dict[int, float]) -> int:
    """Choose the K with the minimum cross-validation error."""
    return min(cv_errors, key=lambda k: cv_errors[k])
