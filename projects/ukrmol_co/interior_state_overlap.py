"""Sparse physical determinant overlaps and conditional projection enclosures.

See docs/physics/co-interior-state-comparison.md. No native decoding, orbital
quadrature, state assignment or molecular qualification is implicit here.
"""

from __future__ import annotations

from itertools import combinations
from math import comb

import numpy as np


def _state(occupations, coefficients, orbitals):
    occupied, vector = np.asarray(occupations), np.asarray(coefficients)
    if (
        occupied.ndim != 2
        or not np.issubdtype(occupied.dtype, np.integer)
        or len(occupied) == 0
        or occupied.shape[1] > 32
        or np.any(occupied < 0)
        or np.any(occupied >= orbitals)
        or np.any(np.diff(occupied, axis=1) <= 0)
        or len(np.unique(occupied, axis=0)) != len(occupied)
    ):
        raise ValueError("Unique canonical integer occupations with at most32 electrons required")
    if vector.shape != (len(occupied),) or not np.isfinite(vector).all():
        raise ValueError("One finite coefficient per determinant required")
    return occupied, vector


def sparse_state_overlap(
    left_occupations,
    left_coefficients,
    right_occupations,
    right_coefficients,
    orbital_metric,
    *,
    max_pairs: int = 1_000_000,
    tile_size: int = 32,
) -> complex:
    """Evaluate bra-left/ket-right occupied minors, including explicit core orbitals.

    Occupations are arrays of sorted spin-orbital indices, not CSF addresses.
    The rectangular orbital cross metric uses the physical Hermitian product
    on a common declared domain. Neither coefficient array is renormalized.
    """
    if (
        type(max_pairs) is not int
        or not 1 <= max_pairs <= 4_000_000
        or type(tile_size) is not int
        or not 1 <= tile_size <= 64
    ):
        raise ValueError("Positive bounded pair cap and tile size required")
    metric = np.asarray(orbital_metric)
    if metric.ndim != 2 or min(metric.shape) < 1 or not np.isfinite(metric).all():
        raise ValueError("Finite rectangular orbital metric required")
    left, a = _state(left_occupations, left_coefficients, metric.shape[0])
    right, b = _state(right_occupations, right_coefficients, metric.shape[1])
    if left.shape[1] != right.shape[1]:
        raise ValueError("Both states must have the same electron count")
    if len(left) * len(right) > max_pairs:
        raise ValueError("Determinant pair cap exceeded before minor allocation")
    result = 0j
    for i in range(0, len(left), tile_size):
        rows = left[i : i + tile_size]
        for j in range(0, len(right), tile_size):
            columns = right[j : j + tile_size]
            minors = metric[rows[:, None, :, None], columns[None, :, None, :]]
            weights = np.linalg.det(minors)
            result += a[i : i + len(rows)].conj() @ weights @ b[j : j + len(columns)]
    if not np.isfinite(result):
        raise ValueError("Nonfinite determinant overlap")
    return complex(result)


def projected_overlap_interval(
    partial_overlap: complex,
    left_capture: float,
    right_capture: float,
    *,
    absolute_error: float = 0.0,
    left_complement_right_norm: float | None = None,
    right_complement_left_norm: float | None = None,
) -> tuple[float, float]:
    """Enclose squared overlap given justified orthogonal-projection captured norms.

    Captures are squared norms in self-orthonormal determinant frames of
    normalized physical states. The caller supplies any justified absolute
    numerical error; this function establishes only the residual algebra.
    """
    if (
        not np.isfinite(partial_overlap)
        or not np.isfinite([left_capture, right_capture, absolute_error]).all()
        or not 0 <= left_capture <= 1
        or not 0 <= right_capture <= 1
        or absolute_error < 0
    ):
        raise ValueError("Finite overlap, unit-interval captures and nonnegative error required")
    a, b = np.sqrt(left_capture), np.sqrt(right_capture)
    magnitude = abs(partial_overlap)
    if magnitude > a * b + absolute_error:
        raise ValueError("Partial overlap violates its physical Cauchy bound")
    ra, rb = np.sqrt(1 - left_capture), np.sqrt(1 - right_capture)
    radius = min(ra + a * rb, rb + b * ra, ra * b + a * rb + ra * rb) + absolute_error
    if (left_complement_right_norm is None) != (right_complement_left_norm is None):
        raise ValueError("Both complementary projected norms are required together")
    if left_complement_right_norm is not None:
        x, y = left_complement_right_norm, right_complement_left_norm
        if not np.isfinite([x, y]).all() or not 0 <= x <= b or not 0 <= y <= a:
            raise ValueError("Complementary projected norms must fit their captured norms")
        radius = min(radius, ra * x + rb * y + ra * rb + absolute_error)
    return float(max(0.0, magnitude - radius) ** 2), float(min(1.0, magnitude + radius) ** 2)


def bound_sector_projection(
    right: np.ndarray,
    orbital_metric: np.ndarray,
    ncore: int,
    active_electrons: tuple[int, int],
) -> np.ndarray:
    """Project into all left bound determinants using full occupied minors.

    Unlike the frozen-core Schur overlap routine, this includes core orbitals
    directly and permits a singular cross-core block. Self frames must still
    be orthonormal for the output norm to represent an orthogonal projection.
    """
    right, metric = map(np.asarray, (right, orbital_metric))
    if (
        metric.ndim != 2
        or metric.shape[0] != metric.shape[1]
        or type(ncore) is not int
        or not 0 <= ncore < len(metric) <= 32
        or len(active_electrons) != 2
        or any(type(n) is not int or not 0 <= n <= len(metric) - ncore for n in active_electrons)
        or not np.isfinite(metric).all()
    ):
        raise ValueError("Bound core/active orbital metric and occupations required")
    occupied = []
    for count in active_electrons:
        strings = comb(len(metric) - ncore, count)
        if strings**2 * max(1, (ncore + count) ** 2) > 20_000_000:
            raise ValueError("Full occupied-minor allocation cap exceeded")
        active = sorted(
            combinations(range(len(metric) - ncore), count), key=lambda c: sum(1 << i for i in c)
        )
        occupied.append(
            np.array([list(range(ncore)) + [ncore + i for i in c] for c in active], int)
        )
    expected = tuple(len(o) for o in occupied)
    if (
        right.ndim != 3
        or right.shape[0] < 1
        or right.shape[1:] != expected
        or not np.isfinite(right).all()
    ):
        raise ValueError("Finite bound-sector CI tensors in determinant order required")
    kernels = []
    for rows in occupied:
        minors = metric[rows[:, None, :, None], rows[None, :, None, :]]
        kernels.append(np.linalg.det(minors))
    return np.stack([kernels[0] @ vector @ kernels[1].T for vector in right])
