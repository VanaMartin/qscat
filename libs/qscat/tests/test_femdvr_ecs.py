"""Analytic benchmark tests for the promoted FEM-DVR-ECS grid (`qscat.dvr`).

Originally developed as a standalone toy model, then ported here once the
modules were promoted into `qscat.dvr`. Four benchmarks:

  B1 (particle-in-box, theta=0) -- exact oracle E_n = n^2 pi^2 / (2 m L^2).
     Extremely sensitive to any assembly error (bridge-weight normalization,
     Dirichlet trim, scatter bookkeeping), so it is the primary arbiter of
     correctness for `kinetic.kinetic`. Includes a spectral-convergence check.
  B2 (harmonic oscillator, theta=0) -- E_n = omega*(n + 1/2). Additionally
     exercises the diagonal-potential DVR approximation.
  B3 (ECS continuum rotation) -- free-particle mid-spectrum eigenvalues
     cluster near arg(E) ~ -2*theta (R0 << Lt asymptote of the exact
     Z_eff = R0 + Lt*e^{i theta} box-quantization formula; `docs/physics/
     femdvr-ecs.md` derives it and explains why this grid's R0/Lt ratio is
     deliberately lopsided).
  B4 (bound-state theta-independence) -- a square-well bound state's energy
     must not depend on the ECS rotation angle theta.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec, eigen, hamiltonian


def _box(L: float = 1.0, nel: int = 4, nq: int = 10) -> FemDvrEcsGrid:
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=[ElementSpec(L / nel)] * nel))


def test_B1_particle_in_box() -> None:
    L, m = 1.0, 1.0
    g = _box(L=L, nel=4, nq=12)
    H = hamiltonian(g, lambda z: 0.0 * z, mass=m)  # V = 0
    E, _ = eigen(H)
    exact = np.array([n**2 * np.pi**2 / (2 * m * L**2) for n in range(1, 6)])
    assert np.allclose(E[:5].real, exact, rtol=1e-6), (E[:5].real, exact)
    assert np.allclose(E[:5].imag, 0.0, atol=1e-9)


def test_B1_spectral_convergence() -> None:
    # nq=(6, 9, 12) as originally specified saturates the double-precision
    # noise floor: measured errors are 1.5e-10 (nq=6), 5.3e-14 (nq=9),
    # 1.2e-13 (nq=12) -- nq=9 and nq=12 are both already at ~1e-13/1e-14,
    # indistinguishable from roundoff, so the monotonicity assertion becomes
    # a coin flip on numerical noise rather than a check of the assembly.
    # Sampling (4, 5, 6) instead stays inside the genuine pre-saturation
    # exponential-convergence regime (errors ~2e-5, ~7e-8, ~1.5e-10).
    L, m = 1.0, 1.0
    err = []
    for nq in (4, 5, 6):
        g = _box(L=L, nel=3, nq=nq)
        E, _ = eigen(hamiltonian(g, lambda z: 0.0 * z, mass=m))
        err.append(abs(E[0].real - np.pi**2 / (2 * m * L**2)))
    assert err[0] > err[1] > err[2]  # error falls as order rises


def test_B2_harmonic_oscillator() -> None:
    m, omega, L = 1.0, 1.0, 20.0
    xc = L / 2
    g = FemDvrEcsGrid(GridSpec(quadrature=10, elements=[ElementSpec(L / 10)] * 10))
    H = hamiltonian(g, lambda z: 0.5 * m * omega**2 * (z - xc) ** 2, mass=m)
    E, _ = eigen(H)
    exact = np.array([omega * (n + 0.5) for n in range(5)])
    assert np.allclose(E[:5].real, exact, rtol=1e-6), (E[:5].real, exact)


def _b3_grid(
    theta_deg: float,
    real_len: float = 1.0,
    nreal: int = 2,
    tail_len: float = 20.0,
    ncomplex: int = 10,
    nq: int = 8,
) -> FemDvrEcsGrid:
    """Short real region (R0=real_len) + long complex tail (Lt=tail_len),
    R0/Lt = 0.05, so the mid-spectrum eigenvalues sit close to the
    R0->0 asymptote arg(E) = -2*theta (see module docstring)."""
    els = [ElementSpec(real_len / nreal)] * nreal + [
        ElementSpec(tail_len / ncomplex, theta_deg)
    ] * ncomplex
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=els))


def _b4_grid(
    theta_deg: float,
    real_len: float = 12.0,
    nreal: int = 4,
    tail_len: float = 12.0,
    ncomplex: int = 6,
    nq: int = 10,
) -> FemDvrEcsGrid:
    """Well edge a=3.0 lands on an element boundary: real_len/nreal = 3.0
    (boundaries at 0, 3, 6, 9, 12); see module docstring."""
    els = [ElementSpec(real_len / nreal)] * nreal + [
        ElementSpec(tail_len / ncomplex, theta_deg)
    ] * ncomplex
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=els))


def test_B3_continuum_rotation() -> None:
    # Free particle on an ECS grid: mid-spectrum eigenvalues cluster near
    # arg(E) ~ -2*theta (see module docstring for the exact Z_eff derivation
    # and why R0 << Lt is required to approach the -2*theta asymptote).
    theta = 30.0
    g = _b3_grid(theta)
    E, _ = eigen(hamiltonian(g, lambda z: 0.0 * z, mass=1.0))
    # pick mid-spectrum eigenvalues with sizeable |E| (avoid ~0 and the
    # top-of-grid numerical-junk states that always appear at the edge of
    # a finite DVR/FEM basis)
    mag = np.abs(E)
    sel = E[(mag > 0.2) & (mag < 5.0)]
    ang = np.degrees(np.angle(sel))
    near = np.abs(ang - (-2 * theta)) < 5.0
    assert sel.size >= 5, "too few mid-spectrum eigenvalues selected"
    assert near.mean() > 0.5, (np.median(ang), -2 * theta)


def test_B4_bound_state_theta_independence() -> None:
    # Square well V=-V0 on [0,a], deep enough for a bound state; energy
    # invariant under theta. a=3.0 lands on an element boundary (see module
    # docstring) so the diagonal-potential DVR represents the discontinuity
    # cleanly.
    a, V0 = 3.0, 5.0

    def Vwell(z: npt.NDArray[np.complex128]) -> npt.NDArray[np.complex128]:
        return np.where(np.real(z) <= a, -V0, 0.0).astype(np.complex128)

    Eb = []
    for theta in (20.0, 35.0):
        g = _b4_grid(theta)
        E, _ = eigen(hamiltonian(g, Vwell, mass=1.0))
        bound = E[E.real < 0].real
        assert bound.size >= 1, "expected a bound state"
        Eb.append(bound.min())
    # Design spec requires rtol <= 1e-6; measured agreement is ~2.8e-14
    # (machine precision), so 1e-8 is a real regression guard with ample
    # headroom rather than a rubber-stamp tolerance.
    assert abs(Eb[0] - Eb[1]) < 1e-8, Eb


def test_spec_rejects_bent_tail_with_multiple_distinct_angles() -> None:
    # more than one distinct nonzero angle_deg among the tail elements is a
    # bent/graded ECS contour, rejected as unvalidated
    elements = [ElementSpec(1.0), ElementSpec(1.0, 30.0), ElementSpec(1.0, 45.0)]
    with pytest.raises(ValueError):
        GridSpec(quadrature=6, elements=elements)
