"""Small Python oracle for the one-continuum-electron overlap factorization.

No native conversion, orbital integration, normalization or state assignment.
See docs/physics/co-full-interior-metric.md; work caps deliberately exclude a
direct molecular FCI expansion.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise

import numpy as np

from projects.ukrmol_co.interior_state_overlap import sparse_state_overlap

DeterminantExpansion = dict[tuple[int, ...], complex]


@dataclass(frozen=True)
class InteriorExpansion:
    """Bound-only Q plus target/continuum coefficients A in a†(u) Phi order."""

    bound: DeterminantExpansion
    targets: tuple[DeterminantExpansion, ...]
    amplitudes: np.ndarray
    bound_orbitals: int
    target_electrons: int


def coefficient_removal(state: DeterminantExpansion, orbital: int) -> DeterminantExpansion:
    """Delete an occupied coordinate with its exterior-algebra permutation sign.

    This is canonical annihilation only in a self-orthonormal orbital frame.
    No complex conjugation or frozen-core exception belongs in this operation.
    """
    result = {}
    for occupied, coefficient in state.items():
        if orbital in occupied:
            index = occupied.index(orbital)
            reduced = occupied[:index] + occupied[index + 1 :]
            result[reduced] = result.get(reduced, 0j) + (-1) ** index * coefficient
    return result


def _validate_state(state, orbitals, electrons):
    if not isinstance(state, dict) or len(state) > 100:
        raise ValueError("A sparse determinant dictionary of at most100 terms is required")
    for occupied, coefficient in state.items():
        if (
            not isinstance(occupied, tuple)
            or len(occupied) != electrons
            or any(type(i) is not int or not 0 <= i < orbitals for i in occupied)
            or any(i >= j for i, j in pairwise(occupied))
            or not np.isfinite(coefficient)
        ):
            raise ValueError(
                "Finite canonical occupations with the declared electron count required"
            )


def _validate(expansion):
    if (
        type(expansion.bound_orbitals) is not int
        or not 1 <= expansion.bound_orbitals <= 10
        or type(expansion.target_electrons) is not int
        or not 0 <= expansion.target_electrons <= 5
        or not 1 <= len(expansion.targets) <= 6
    ):
        raise ValueError("Toy bound-orbital, electron or target cap exceeded")
    amplitudes = np.asarray(expansion.amplitudes)
    if (
        amplitudes.ndim != 2
        or amplitudes.shape[0] != len(expansion.targets)
        or not 1 <= amplitudes.shape[1] <= 6
        or not np.isfinite(amplitudes).all()
    ):
        raise ValueError("Finite target-by-continuum amplitudes within toy caps required")
    _validate_state(expansion.bound, expansion.bound_orbitals, expansion.target_electrons + 1)
    for target in expansion.targets:
        _validate_state(target, expansion.bound_orbitals, expansion.target_electrons)
    return amplitudes


def four_block_overlap(
    left: InteriorExpansion,
    right: InteriorExpansion,
    orbital_metric: np.ndarray,
    *,
    max_pairs: int = 1_000_000,
) -> dict[str, complex]:
    """Return coherent QQ, PQ, QP, PP and total amplitudes, without rescaling.

    The rectangular metric orders bound then continuum spin orbitals in each
    frame. Removal-state overlaps evaluate cofactors even for singular cross
    metrics. A shared pair budget bounds all determinant-minor contractions.
    """
    a, b = _validate(left), _validate(right)
    if left.target_electrons != right.target_electrons:
        raise ValueError("Both states must have the same electron count")
    if type(max_pairs) is not int or not 1 <= max_pairs <= 1_000_000:
        raise ValueError("A positive bounded determinant pair budget is required")
    na, nb = left.bound_orbitals, right.bound_orbitals
    metric = np.asarray(orbital_metric)
    if metric.shape != (na + a.shape[1], nb + b.shape[1]) or not np.isfinite(metric).all():
        raise ValueError("Finite bound-plus-continuum cross metric with matching shape required")
    bb, bc, cb, cc = metric[:na, :nb], metric[:na, nb:], metric[na:, :nb], metric[na:, nb:]
    remaining = max_pairs

    def inner(x, y):
        nonlocal remaining
        if not x or not y:
            return 0j
        remaining -= len(x) * len(y)
        if remaining < 0:
            raise ValueError("Total determinant pair budget exceeded")
        return sparse_state_overlap(
            np.asarray(list(x), int),
            np.asarray(list(x.values())),
            np.asarray(list(y), int),
            np.asarray(list(y.values())),
            bb,
        )

    remove_q_a = [coefficient_removal(left.bound, mu) for mu in range(na)]
    remove_q_b = [coefficient_removal(right.bound, nu) for nu in range(nb)]
    qq = inner(left.bound, right.bound)
    pq = qp = pp = 0j
    for i, phi_a in enumerate(left.targets):
        # a†(u) Phi on the bra: <Phi|D_nu Q> S_cb[u,nu].
        dyson = np.array([inner(phi_a, state) for state in remove_q_b])
        pq += a[i].conj() @ cb @ dyson
        removed_a = [coefficient_removal(phi_a, mu) for mu in range(na)]
        for j, phi_b in enumerate(right.targets):
            removed_b = [coefficient_removal(phi_b, nu) for nu in range(nb)]
            target = inner(phi_a, phi_b)
            removal = np.array([[inner(x, y) for y in removed_b] for x in removed_a])
            pp += a[i].conj() @ (cc * target - cb @ removal.T @ bc) @ b[j]
    for j, phi_b in enumerate(right.targets):
        # The reverse mixed block, expressed in this cross metric without a
        # real-coefficient or equal-geometry assumption.
        dyson = np.array([inner(state, phi_b) for state in remove_q_a])
        qp += dyson @ bc @ b[j]
    return {
        "qq": complex(qq),
        "pq": complex(pq),
        "qp": complex(qp),
        "pp": complex(pp),
        "total": complex(qq + pq + qp + pp),
    }
