"""Coherent spin-separated core-sector QQ/PQ/QP/PP amplitudes, without native IO."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from projects.ukrmol_co.core_sector_overlap import (
    CoreSectorCI,
    CoreSectorMetric,
    removal_overlaps,
    remove_spin_orbital,
)


@dataclass(frozen=True)
class TargetProducts:
    """One target spin-projection family in explicit a-dagger(continuum) Phi order."""

    targets: CoreSectorCI
    continuum_spin: int
    amplitudes: np.ndarray


@dataclass(frozen=True)
class CoreInteriorExpansion:
    """One bound-only Qroot plus up to two coherent continuum/target families."""

    bound: CoreSectorCI
    products: tuple[TargetProducts, ...]


def _validate(state, metric):
    if not isinstance(state, CoreInteriorExpansion) or not isinstance(state.products, tuple):
        raise ValueError("Explicit core-sector interior expansion required")
    q = metric._check(state.bound)
    if len(q) != 1 or not 1 <= len(state.products) <= 2:
        raise ValueError("One bound-only root and one or two target families required")
    q_counts = tuple(
        len(c) + n
        for c, n in zip(
            (state.bound.core_alpha, state.bound.core_beta),
            state.bound.active_electrons,
            strict=True,
        )
    )
    continuum = None
    spins = []
    for family in state.products:
        if not isinstance(family, TargetProducts):
            raise ValueError("Explicit target-product family required")
        c = metric._check(family.targets)
        spin = family.continuum_spin
        a = np.asarray(family.amplitudes)
        if (
            type(spin) is not int
            or spin not in (0, 1)
            or a.ndim != 2
            or a.shape[0] != len(c)
            or not 1 <= a.shape[1] <= 160
            or not np.isfinite(a).all()
        ):
            raise ValueError("Finite bounded target-by-continuum amplitudes and spin required")
        counts = [
            len(core) + n
            for core, n in zip(
                (family.targets.core_alpha, family.targets.core_beta),
                family.targets.active_electrons,
                strict=True,
            )
        ]
        counts[spin] += 1
        if tuple(counts) != q_counts:
            raise ValueError("Each target/continuum product must match Q spin-electron counts")
        if continuum is not None and continuum != a.shape[1]:
            raise ValueError("Each family must use the same spatial continuum inventory")
        continuum = a.shape[1]
        spins.append(spin)
    if len(spins) != len(set(spins)):
        raise ValueError("Combine targets with the same continuum spin into one family")
    return continuum


def core_sector_four_block_overlap(
    left: CoreInteriorExpansion,
    right: CoreInteriorExpansion,
    bound_metric: CoreSectorMetric,
    bound_continuum: np.ndarray,
    continuum_bound: np.ndarray,
    continuum_metric: np.ndarray,
) -> dict[str, complex]:
    """Return unnormalized coherent blocks in one declared spin-preserving domain.

    Spatial bb is owned by bound_metric; bc/cb/cc need not be zero, equal or
    Hermitian for different frames. Native phases belong in the unchanged input
    CI/amplitudes, not in a guessed Clebsch-Gordan conversion inside this kernel.
    """
    if not isinstance(bound_metric, CoreSectorMetric):
        raise ValueError("Explicit core-sector spatial bound cross metric required")
    na, nb = _validate(left, bound_metric), _validate(right, bound_metric)
    bc, cb, cc = map(np.asarray, (bound_continuum, continuum_bound, continuum_metric))
    bound = bound_metric.ncore + bound_metric.nactive
    if (
        bc.shape != (bound, nb)
        or cb.shape != (na, bound)
        or cc.shape != (na, nb)
        or not all(np.isfinite(v).all() for v in (bc, cb, cc))
    ):
        raise ValueError("Finite matching spatial bb/bc/cb/cc cross-metric blocks required")
    qq = bound_metric.overlaps(left.bound, right.bound)[0, 0]
    pq = qp = pp = 0j
    for family in left.products:
        a = np.asarray(family.amplitudes)
        left_projection = a.conj() @ cb
        for nu in range(bound):
            removed = remove_spin_orbital(right.bound, family.continuum_spin, nu)
            if removed is not None:
                dyson = bound_metric.overlaps(family.targets, removed)[:, 0]
                pq += left_projection[:, nu] @ dyson
        for other in right.products:
            b = np.asarray(other.amplitudes)
            if family.continuum_spin == other.continuum_spin:
                target = bound_metric.overlaps(family.targets, other.targets)
                pp += np.vdot(a, target @ b @ cc.T)
            removal = removal_overlaps(
                family.targets,
                other.targets,
                other.continuum_spin,
                bound_metric,
                right_spin=family.continuum_spin,
            )
            right_projection = b @ bc.T
            pp -= np.einsum(
                "in,minj,jm->", left_projection, removal, right_projection, optimize=True
            )
            del removal
    for family in right.products:
        b = np.asarray(family.amplitudes)
        right_projection = b @ bc.T
        for mu in range(bound):
            removed = remove_spin_orbital(left.bound, family.continuum_spin, mu)
            if removed is not None:
                dyson = bound_metric.overlaps(removed, family.targets)[0]
                qp += dyson @ right_projection[:, mu]
    return dict(
        qq=complex(qq),
        pq=complex(pq),
        qp=complex(qp),
        pp=complex(pp),
        total=complex(qq + pq + qp + pp),
    )
