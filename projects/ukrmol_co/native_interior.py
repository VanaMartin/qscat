"""Explicit native P/Q row conversion and retained-eigenbasis correspondence.

No Hamiltonian action, eigensolve, phase inference from energies alone or missing
coefficient recovery is performed. See docs/physics/co-native-interior-recovery.md.
"""

from __future__ import annotations

from itertools import pairwise
from math import comb

import numpy as np
from scipy.sparse import isspmatrix_csr

from projects.ukrmol_co.core_sector_full_overlap import CoreInteriorExpansion
from projects.ukrmol_co.core_sector_overlap import CoreSectorCI, _validate
from projects.ukrmol_co.native_target_products import native_target_products


def continuum_row_layout(
    orbital_irreps, flags, target_counts, continuum_irreps, continuum_counts
) -> tuple:
    """Explicit native MO flags -> per-sector zero-based continuum coordinates."""
    irreps, flags = np.asarray(orbital_irreps), np.asarray(flags)
    if (
        irreps.ndim != 1
        or not 1 <= len(irreps) <= 512
        or flags.shape != irreps.shape
        or not np.isin(irreps, range(1, 9)).all()
        or not np.isin(flags, [0, 1]).all()
        or not 1 <= np.count_nonzero(flags) <= 160
        or not isinstance(target_counts, tuple)
        or not 1 <= len(target_counts) <= 8
        or len(continuum_irreps) != len(target_counts)
        or len(continuum_counts) != len(target_counts)
        or any(type(n) is not int or not 1 <= n <= 64 for n in target_counts)
        or sum(target_counts) > 64
    ):
        raise ValueError("Bounded explicit MO/target/continuum inventories required")
    continuum = irreps[flags == 1]
    result = []
    for symmetry, count in zip(continuum_irreps, continuum_counts, strict=True):
        if (
            type(symmetry) is not int
            or not 0 <= symmetry < 8
            or type(count) is not int
            or not 1 <= count <= 160
        ):
            raise ValueError("Bounded zero-based native continuum symmetry/count required")
        indices = tuple(int(i) for i in np.flatnonzero(continuum == symmetry + 1))
        if len(indices) != count:
            raise ValueError("Executed continuum count disagrees with native MO flags")
        result.append(indices)
    return tuple(result)


def assemble_native_interior(
    sectors, continuum_indices, coefficients, bound_transform, ncontinuum
) -> CoreInteriorExpansion:
    """Map supplied complete P-then-Q coefficients, never substitute missing Q."""
    if not isinstance(sectors, tuple) or not sectors:
        raise ValueError("Explicit raw-target sectors required")
    first = sectors[0]
    _validate(first)
    total = sum(first.active_electrons)
    if total % 2:
        raise ValueError("Even-electron target required")
    electrons = (total // 2 + 1, total // 2)
    shape = tuple(comb(first.nactive, n) for n in electrons)
    if (
        not isspmatrix_csr(bound_transform)
        or bound_transform.shape[0] != np.prod(shape)
        or not 1 <= bound_transform.shape[1] <= 50_000
        or bound_transform.nnz > 2_000_000
        or not np.isfinite(bound_transform.data).all()
    ):
        raise ValueError("Qualified bounded sparse bound-CSF transform required")
    c = np.asarray(coefficients)
    if c.ndim != 1 or c.size > 100_000 or not np.isfinite(c).all():
        raise ValueError("Complete finite bounded native coefficient vector required")
    q_columns = bound_transform.shape[1]
    if len(c) <= q_columns:
        raise ValueError("Both supplied P and Q coefficient inventories required")
    products = native_target_products(sectors, continuum_indices, c[:-q_columns], ncontinuum)
    q = bound_transform @ c[-q_columns:]
    core = tuple(range(first.ncore))
    bound = CoreSectorCI(
        first.ncore, first.nactive, core, core, electrons, np.asarray(q).reshape(1, *shape)
    )
    _validate(bound)
    return CoreInteriorExpansion(bound, products)


def recover_retained_state(
    old_energies,
    old_continuum,
    matching,
    new_energies,
    new_vectors,
    *,
    energy_tolerance=1e-9,
    cluster_gap=1e-8,
    component_tolerance=1e-8,
    max_cluster=16,
) -> tuple[np.ndarray, dict]:
    """Transport accepted complex c through identifiable real spectral rotations.

    Roots are rows. For each energy cluster, old_P = R new_P and c_new = R.T c.
    Invisible/rank-deficient clusters reject; no coefficients are discarded.
    """
    old_e, old_p, c, new_e, new_u = map(
        np.asarray, (old_energies, old_continuum, matching, new_energies, new_vectors)
    )
    if (
        old_e.ndim != 1
        or not 1 <= len(old_e) <= 4096
        or new_e.shape != old_e.shape
        or old_p.ndim != 2
        or old_p.shape[0] != len(old_e)
        or not 1 <= old_p.shape[1] <= 4096
        or new_u.ndim != 2
        or new_u.shape[0] != len(old_e)
        or not old_p.shape[1] <= new_u.shape[1] <= 50_000
        or c.shape != old_e.shape
        or not all(np.isfinite(x).all() for x in (old_e, old_p, c, new_e, new_u))
        or any(not np.isrealobj(x) for x in (old_e, old_p, new_e, new_u))
        or np.any(np.diff(old_e) < 0)
        or np.any(np.diff(new_e) < 0)
        or type(max_cluster) is not int
        or not 1 <= max_cluster <= 64
        or any(
            not np.isfinite(t) or t <= 0
            for t in (energy_tolerance, cluster_gap, component_tolerance)
        )
    ):
        raise ValueError(
            "Bounded finite energy-ordered real spectral frames and complex matching c required"
        )
    energy_error = float(np.max(abs(new_e - old_e)))
    if energy_error > energy_tolerance:
        raise ValueError("Recovered energies disagree with retained native frame")
    boundaries = [0, *list(np.flatnonzero(np.diff(old_e) > cluster_gap) + 1), len(old_e)]
    if max(np.diff(boundaries)) > max_cluster:
        raise ValueError("Retained spectral cluster cap exceeded")
    result = np.zeros(new_u.shape[1], dtype=np.result_type(c, complex))
    maximum = 0.0
    clusters = []
    for start, stop in pairwise(boundaries):
        a = old_p[start:stop]
        b = new_u[start:stop, : old_p.shape[1]]
        visibility = np.linalg.svd(a, compute_uv=False)
        if visibility[-1] <= max(1e-12, visibility[0] * 1e-10):
            raise ValueError("Invisible or rank-deficient retained cluster cannot be aligned")
        left, _, right = np.linalg.svd(a @ b.T, full_matrices=False)
        rotation = left @ right
        error = float(np.max(abs(a - rotation @ b)))
        if error > component_tolerance:
            raise ValueError("Recovered continuum components fail spectral alignment")
        maximum = max(maximum, error)
        clusters.append(stop - start)
        result += new_u[start:stop].T @ (rotation.T @ c[start:stop])
    return result, dict(
        clusters=clusters,
        maximum_energy_error=energy_error,
        maximum_component_error=maximum,
        coefficients_discarded=0,
    )
