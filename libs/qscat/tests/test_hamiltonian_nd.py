"""Analytic validation of the N-dimensional Hamiltonian (V3) and the D=1
regression against the existing 1-D stack (V4).

The generality is EXERCISED, not asserted: every benchmark runs at D = 1, 2
and 3, with unequal per-axis extents and masses so a transposed-axis bug
cannot pass.
"""

from __future__ import annotations

import itertools
from collections.abc import Callable

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from qscat.dvr import (
    ElementSpec,
    FemDvrEcsGrid,
    GridSpec,
    TensorGrid,
    hamiltonian,
    hamiltonian_nd,
    kinetic,
    kinetic_nd,
    potential_nd,
)

Floats = npt.NDArray[np.float64]


def _box_grid(length: float, n_el: int, q: int, x_min: float = 0.0) -> FemDvrEcsGrid:
    """All-real grid tiling [x_min, x_min + length] with `n_el` equal elements."""
    return FemDvrEcsGrid(
        GridSpec(
            quadrature=q,
            elements=[ElementSpec(length / n_el) for _ in range(n_el)],
            x_min=x_min,
        )
    )


def _lowest_dense(H: sp.csr_matrix, k: int) -> Floats:
    """Lowest `k` eigenvalues by real part, via a dense eigensolve."""
    vals = np.linalg.eigvals(H.toarray())
    return np.asarray(np.sort(vals.real)[:k], dtype=np.float64)


def _lowest_shift_invert(H: sp.csr_matrix, k: int, sigma: float) -> Floats:
    """Lowest `k` eigenvalues near `sigma`, via sparse shift-invert."""
    vals = spla.eigs(H.tocsc(), k=k, sigma=sigma, return_eigenvectors=False)
    return np.asarray(np.sort(vals.real)[:k], dtype=np.float64)


def _analytic_box_levels(
    lengths: tuple[float, ...], masses: tuple[float, ...], k: int
) -> Floats:
    """Lowest `k` levels of sum_d n_d^2 pi^2 / (2 m_d L_d^2)."""
    per_axis = [
        [n**2 * np.pi**2 / (2.0 * m * L**2) for n in range(1, k + 2)]
        for L, m in zip(lengths, masses, strict=True)
    ]
    sums = [sum(c) for c in itertools.product(*per_axis)]
    return np.asarray(np.sort(np.array(sums))[:k], dtype=np.float64)


def _analytic_oscillator_levels(omegas: tuple[float, ...], k: int) -> Floats:
    """Lowest `k` levels of sum_d omega_d (n_d + 1/2)."""
    per_axis = [[w * (n + 0.5) for n in range(k + 1)] for w in omegas]
    sums = [sum(c) for c in itertools.product(*per_axis)]
    return np.asarray(np.sort(np.array(sums))[:k], dtype=np.float64)


# --------------------------------------------------------------------------
# V3a: particle in a D-dimensional box  (kinetic_nd alone -- no potential)
# --------------------------------------------------------------------------


def test_box_d1() -> None:
    L, m = (1.0,), (1.0,)
    tg = TensorGrid([_box_grid(L[0], 4, 8)])
    got = _lowest_dense(kinetic_nd(tg, m), 4)
    want = _analytic_box_levels(L, m, 4)
    # measured 3.4e-9 at this basis size (4 elements x 8-point quadrature)
    assert np.allclose(got, want, rtol=1e-8, atol=0)


def test_box_d2_unequal_extents_and_masses() -> None:
    L, m = (1.0, 1.3), (1.0, 2.0)
    tg = TensorGrid([_box_grid(L[0], 4, 8), _box_grid(L[1], 4, 8)])
    got = _lowest_dense(kinetic_nd(tg, m), 5)
    want = _analytic_box_levels(L, m, 5)
    # measured 3.8e-12 at this basis size
    assert np.allclose(got, want, rtol=1e-11, atol=0)


def test_box_d3_unequal_extents_and_masses() -> None:
    L, m = (1.0, 1.3, 0.9), (1.0, 2.0, 1.5)
    tg = TensorGrid([_box_grid(L[d], 2, 6) for d in range(3)])
    got = _lowest_dense(kinetic_nd(tg, m), 4)
    want = _analytic_box_levels(L, m, 4)
    # measured 5.4e-5 at this basis size (2 elements x 6-point quadrature per
    # axis is deliberately coarse to keep the D=3 dense eigensolve fast)
    assert np.allclose(got, want, rtol=1.5e-4, atol=0)


# --------------------------------------------------------------------------
# V3b: harmonic oscillator  (potential_nd + hamiltonian_nd)
# --------------------------------------------------------------------------


def _oscillator_V(
    masses: tuple[float, ...], omegas: tuple[float, ...]
) -> Callable[..., npt.ArrayLike]:
    """V = sum_d 0.5 m_d omega_d^2 x_d^2, for the broadcastable coords of points()."""

    def V(*coords: npt.NDArray[np.complex128]) -> npt.ArrayLike:
        total: npt.NDArray[np.complex128] = np.zeros((), dtype=np.complex128)
        for x, m, w in zip(coords, masses, omegas, strict=True):
            total = total + 0.5 * m * w**2 * x**2
        return total

    return V


def test_oscillator_d1() -> None:
    m, w = (1.0,), (1.0,)
    tg = TensorGrid([_box_grid(16.0, 8, 10, x_min=-8.0)])
    H = hamiltonian_nd(tg, m, _oscillator_V(m, w))
    got = _lowest_dense(H, 4)
    want = _analytic_oscillator_levels(w, 4)
    # measured 6.7e-9 at this basis size (8 elements x 10-point quadrature)
    assert np.allclose(got, want, rtol=2e-8, atol=0)


def test_oscillator_d2_unequal_frequencies() -> None:
    """Each axis gets a genuinely different grid (extent, element count AND
    quadrature order, so the two axes even have different point counts: 27 vs
    31) -- not just a different omega -- so a transposed-axis bug cannot pass
    by having both axes look identical.
    """
    m, w = (1.0, 1.0), (1.0, 1.7)
    tg = TensorGrid([_box_grid(14.0, 4, 8, x_min=-7.0), _box_grid(16.0, 4, 9, x_min=-8.0)])
    H = hamiltonian_nd(tg, m, _oscillator_V(m, w))
    got = _lowest_dense(H, 4)
    want = _analytic_oscillator_levels(w, 4)
    # measured 6.2e-4 at this basis size (finite box + finite basis)
    assert np.allclose(got, want, rtol=1.5e-3, atol=0)


def test_oscillator_d3_unequal_frequencies() -> None:
    """Each axis gets a distinct grid (extent, element count and quadrature
    order all differ -> distinct point counts 15, 17, 13), so a transposed-
    axis bug cannot pass by having any two axes look alike.
    """
    m, w = (1.0, 1.0, 1.0), (1.0, 1.3, 1.7)
    tg = TensorGrid(
        [
            _box_grid(12.0, 2, 9, x_min=-6.0),
            _box_grid(13.0, 2, 10, x_min=-6.5),
            _box_grid(11.0, 2, 8, x_min=-5.5),
        ]
    )
    H = hamiltonian_nd(tg, m, _oscillator_V(m, w))
    want = _analytic_oscillator_levels(w, 3)
    got = _lowest_shift_invert(H, 3, sigma=float(want[0]) - 0.3)
    # measured 8.3e-4 at this basis size (coarse boxes kept deliberately small
    # so the D=3 sparse shift-invert eigensolve stays fast)
    assert np.allclose(got, want, rtol=2e-3, atol=0)


# --------------------------------------------------------------------------
# V4: D = 1 reproduces the existing 1-D stack bit-for-bit
# --------------------------------------------------------------------------


def test_d1_reproduces_existing_1d_hamiltonian_with_ecs() -> None:
    """Makes every sub-project #1-#4 result a regression test on this code."""
    els = [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0), ElementSpec(2.0, 35.0)]
    g = FemDvrEcsGrid(GridSpec(quadrature=9, elements=els, x_min=0.0))
    mass = 12766.36

    def V(z: npt.NDArray[np.complex128]) -> npt.ArrayLike:
        return 1.0 / (1.0 + z**2)

    H_dense = hamiltonian(g, V, mass)
    H_nd = hamiltonian_nd(TensorGrid([g]), [mass], V).toarray()

    scale = np.abs(H_dense).max()
    assert np.abs(H_nd - H_dense).max() <= 1e-13 * scale


def test_d1_kinetic_nd_matches_dense_kinetic() -> None:
    g = FemDvrEcsGrid(
        GridSpec(quadrature=7, elements=[ElementSpec(1.0) for _ in range(3)])
    )
    dense = kinetic(g, 2.5)
    got = kinetic_nd(TensorGrid([g]), [2.5]).toarray()
    assert np.abs(got - dense).max() <= 1e-13 * np.abs(dense).max()


# --------------------------------------------------------------------------
# Structural: ECS makes H complex symmetric, never Hermitian
# --------------------------------------------------------------------------


def test_hamiltonian_nd_is_complex_symmetric_under_ecs() -> None:
    """Two DIFFERENT ECS grids (different real-region extent -- 2.0 vs 4.5 --
    and a different number of complex tail elements -- 1 vs 2, hence
    different sizes: 20 vs 34 basis functions) and an asymmetric potential
    (different coefficients on each coordinate, different mass per axis), so
    a transposed-axis bug could not masquerade as a still-symmetric matrix.
    """
    els_a = [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0)]
    g_a = FemDvrEcsGrid(GridSpec(quadrature=8, elements=els_a))
    els_b = [
        ElementSpec(1.5),
        ElementSpec(1.5),
        ElementSpec(1.5),
        ElementSpec(1.0, 35.0),
        ElementSpec(1.0, 35.0),
    ]
    g_b = FemDvrEcsGrid(GridSpec(quadrature=8, elements=els_b))
    tg = TensorGrid([g_a, g_b])
    H = hamiltonian_nd(tg, [1.0, 2.0], lambda a, b: 1.0 / (1.0 + 2.0 * a**2 + 3.0 * b**2))
    assert abs(H - H.T).max() < 1e-10
    assert abs(H - H.conj().T).max() > 1e-6


def test_potential_nd_preserves_complex_points_on_the_ecs_tail() -> None:
    """A potential that coerced to float would silently kill the ECS tail."""
    els = [ElementSpec(1.0), ElementSpec(1.0), ElementSpec(2.0, 35.0)]
    g = FemDvrEcsGrid(GridSpec(quadrature=8, elements=els))
    tg = TensorGrid([g])
    vals = potential_nd(tg, lambda z: z**2)
    assert np.abs(vals.imag).max() > 1e-3
