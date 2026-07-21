"""Analytic benchmark tests for the FEM-DVR kinetic operator + Hamiltonian/eigen helpers.

B1 (particle-in-box) is an exact oracle: E_n = n^2 pi^2 / (2 m L^2). It is extremely
sensitive to any assembly error (bridge-weight normalization, Dirichlet trim, scatter
bookkeeping), so it is the primary arbiter of correctness for `kinetic.kinetic`.

B2 (harmonic oscillator) additionally exercises the diagonal-potential DVR approximation.

See .superpowers/sdd/femdvr-ecs-extraction.md and task-2-brief.md.
"""

import numpy as np
from spec import ElementSpec, GridSpec
from grid import FemDvrEcsGrid
from operators import hamiltonian, eigen


def _box(L=1.0, nel=4, nq=10):
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=[ElementSpec(L / nel)] * nel))


def test_B1_particle_in_box():
    L, m = 1.0, 1.0
    g = _box(L=L, nel=4, nq=12)
    H = hamiltonian(g, lambda z: 0.0 * z, mass=m)   # V = 0
    E, _ = eigen(H)
    exact = np.array([n**2 * np.pi**2 / (2 * m * L**2) for n in range(1, 6)])
    assert np.allclose(E[:5].real, exact, rtol=1e-6), (E[:5].real, exact)
    assert np.allclose(E[:5].imag, 0.0, atol=1e-9)


def test_B1_spectral_convergence():
    # nq=(6, 9, 12) as originally specified saturates the double-precision
    # noise floor: measured errors are 1.5e-10 (nq=6), 5.3e-14 (nq=9),
    # 1.2e-13 (nq=12) -- nq=9 and nq=12 are both already at ~1e-13/1e-14,
    # indistinguishable from roundoff, so the monotonicity assertion becomes
    # a coin flip on numerical noise rather than a check of the assembly.
    # Sampling (4, 5, 6) instead stays inside the genuine pre-saturation
    # exponential-convergence regime (errors ~2e-5, ~7e-8, ~1.5e-10), which
    # is what this test is meant to exercise. This is the only test-spec
    # deviation from the brief; the exact-oracle B1 test and its rtol=1e-6
    # are untouched.
    L, m = 1.0, 1.0
    err = []
    for nq in (4, 5, 6):
        g = _box(L=L, nel=3, nq=nq)
        E, _ = eigen(hamiltonian(g, lambda z: 0.0 * z, mass=m))
        err.append(abs(E[0].real - np.pi**2 / (2 * m * L**2)))
    assert err[0] > err[1] > err[2]                 # error falls as order rises


def test_B2_harmonic_oscillator():
    m, omega, L = 1.0, 1.0, 20.0
    xc = L / 2
    g = FemDvrEcsGrid(GridSpec(quadrature=10, elements=[ElementSpec(L / 10)] * 10))
    H = hamiltonian(g, lambda z: 0.5 * m * omega**2 * (z - xc) ** 2, mass=m)
    E, _ = eigen(H)
    exact = np.array([omega * (n + 0.5) for n in range(5)])
    assert np.allclose(E[:5].real, exact, rtol=1e-6), (E[:5].real, exact)
