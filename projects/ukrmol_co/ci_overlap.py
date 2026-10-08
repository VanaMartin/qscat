"""Frozen-core determinant-CI overlaps in a nonorthogonal orbital metric.

CI arrays use ascending integer-bitstring determinant order within each spin.
The overlap includes doubly occupied core orbitals and core/active cross terms.
"""

from __future__ import annotations

from itertools import combinations
from math import comb

import numpy as np


def _occupied(norb: int, electrons: int) -> np.ndarray:
    strings = sorted(combinations(range(norb), electrons), key=lambda s: sum(1 << i for i in s))
    return np.asarray(strings, dtype=int).reshape(comb(norb, electrons), electrons)


def _determinant_metric(metric: np.ndarray, electrons: int) -> np.ndarray:
    occupations = _occupied(len(metric), electrons)
    if not electrons:
        return np.ones((1, 1), dtype=metric.dtype)
    minors = metric[occupations[:, None, :, None], occupations[None, :, None, :]]
    return np.linalg.det(minors)


def frozen_core_overlaps(
    left: np.ndarray,
    right: np.ndarray,
    orbital_overlap: np.ndarray,
    ncore: int,
    active_electrons: tuple[int, int],
) -> np.ndarray:
    """Return the bra-left/ket-right root overlap matrix at fixed core occupation.

    ``left``/``right`` have shape (roots, alpha strings, beta strings).
    ``orbital_overlap`` orders core then active orbitals at both geometries.
    Each geometry's orbitals must already be orthonormal in its own AO metric.
    The core cross block must be nonsingular and well conditioned; no projection
    normalization is applied to lost overlap or to the finite root manifold.
    """
    left, right, metric = map(np.asarray, (left, right, orbital_overlap))
    if (
        metric.ndim != 2
        or metric.shape[0] != metric.shape[1]
        or not isinstance(ncore, int)
        or not 0 <= ncore < len(metric)
        or len(active_electrons) != 2
    ):
        raise ValueError("A square core-plus-active metric and valid occupations are required")
    ncas = len(metric) - ncore
    if any(not isinstance(n, int) or not 0 <= n <= ncas for n in active_electrons):
        raise ValueError("Active electron counts must fit the active orbital space")
    expected = tuple(comb(ncas, n) for n in active_electrons)
    if any(a.ndim != 3 or a.shape[1:] != expected or a.shape[0] < 1 for a in (left, right)):
        raise ValueError("CI arrays must contain roots in alpha/beta determinant order")
    if not all(np.all(np.isfinite(a)) for a in (left, right, metric)):
        raise ValueError("CI arrays and orbital overlaps must be finite")
    if ncore:
        core = metric[:ncore, :ncore]
        if np.linalg.cond(core) > 1e10:
            raise ValueError("Core cross block is singular or ill conditioned")
        effective = metric[ncore:, ncore:] - metric[ncore:, :ncore] @ np.linalg.solve(
            core, metric[:ncore, ncore:]
        )
        core_factor = np.linalg.det(core) ** 2
    else:
        effective = metric
        core_factor = 1.0
    alpha = _determinant_metric(effective, active_electrons[0])
    beta = _determinant_metric(effective, active_electrons[1])
    transformed = np.stack([alpha @ vector @ beta.T for vector in right])
    return core_factor * left.reshape(len(left), -1).conj() @ transformed.reshape(len(right), -1).T
