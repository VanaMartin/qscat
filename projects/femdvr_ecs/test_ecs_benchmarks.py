"""ECS validation benchmarks for the FEM-DVR-ECS grid.

B3 (continuum rotation) and B4 (bound-state theta-independence) validate the
exterior-complex-scaling (ECS) machinery in `grid.py`/`kinetic.py` -- the
`e^{i*theta}` Jacobian in the complex half-lengths `hz` and its propagation
through the kinetic-energy assembly. B1/B2 (see test_kinetic_benchmarks.py)
already trust the theta=0 machinery to high precision; this file exercises
theta != 0.

See .superpowers/sdd/task-3-brief.md and femdvr-ecs-extraction.md.

--- B3 derivation (why the grid parameters below, not the brief's literal
    equal-split example) ---

For a free particle (V=0) on an ECS grid with a real region of length R0
followed by a single complex tail of length Lt at angle theta, and Dirichlet
(psi=0) at both the origin and the far tail edge, the exact eigenvalue
problem is analytically solvable (match psi and psi' at z=R0, where the
tail's local coordinate is z = R0 + rho*e^{i*theta}):

    psi_real(z) = sin(k z),                       z in [0, R0]
    psi_tail(rho) = cos(k R0) sin(k e^{i theta} rho) + sin(k R0) cos(k e^{i theta} rho)
                  = sin(k R0 + k e^{i theta} rho)

Dirichlet at rho = Lt gives the quantization condition

    k (R0 + Lt e^{i theta}) = n*pi   =>   E_n = n^2 pi^2 / (2 m Z_eff^2),
    Z_eff := R0 + Lt * e^{i theta}

This is EXACT (not asymptotic) for every low-lying box mode: arg(E_n) =
-2*arg(Z_eff) for ALL n, independent of n, because Z_eff does not depend on n.

The textbook "continuum rotates to arg(E) = -2*theta" picture is the R0=0 (or
R0 << Lt) limit of this formula, where arg(Z_eff) -> theta. With an EQUAL
real/complex split (R0 == Lt, e.g. the brief's literal nreal=ncomplex=8
suggestion), arg(Z_eff) = theta/2 exactly (since 1 + e^{i theta} =
2*cos(theta/2)*e^{i theta/2} independent of the common length), giving
arg(E) = -theta, NOT -2*theta -- verified numerically: that grid clusters
mid-spectrum eigenvalues tightly at exactly -theta, and the test as given in
the brief would fail its own `near.mean() > 0.5` assertion (fraction was
0.0, all mass sitting at -theta instead of -2*theta).

This is a test-setup/physics artifact, not an ECS Jacobian bug: the "pure
rotation" sanity check (all elements at the same theta, no real region --
see the derivation check performed during development) reproduces the exact
box-with-scaled-length analytic spectrum E_n = n^2 pi^2 e^{-2i theta} /
(2 m L^2) to machine precision, confirming `kinetic.py`'s `e^{i theta}`
handling is correct.

So B3 below uses a SHORT real region and a LONG complex tail (R0/Lt = 1/20 =
0.05), for which arg(Z_eff) = 28.63 deg (vs the theta=30 deg limit), i.e.
arg(E) = -57.25 deg (vs the -60 deg asymptote) -- well within the brief's
+/-5 deg window, and confirmed numerically to place ~100% of the
mid-spectrum window's eigenvalues within that window (see development notes
in task-3-report.md).

--- B4 note ---

The square well's discontinuity at r=a must land on an element boundary for
the diagonal-potential DVR to represent it cleanly (the note in the brief
and femdvr-ecs-extraction.md): a=3.0 is placed at a real-element edge by
choosing real_len=12.0, nreal=4 (element length 12/4=3.0, so boundaries are
at 0, 3, 6, 9, 12 -- a=3.0 sits exactly on the first internal boundary).
"""

import numpy as np
from spec import ElementSpec, GridSpec
from grid import FemDvrEcsGrid
from operators import hamiltonian, eigen


def _b3_grid(theta_deg, real_len=1.0, nreal=2, tail_len=20.0, ncomplex=10, nq=8):
    """Short real region (R0=real_len) + long complex tail (Lt=tail_len),
    R0/Lt = 0.05, so the mid-spectrum eigenvalues sit close to the
    R0->0 asymptote arg(E) = -2*theta (see module docstring)."""
    els = (
        [ElementSpec(real_len / nreal)] * nreal
        + [ElementSpec(tail_len / ncomplex, theta_deg)] * ncomplex
    )
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=els))


def _b4_grid(theta_deg, real_len=12.0, nreal=4, tail_len=12.0, ncomplex=6, nq=10):
    """Well edge a=3.0 lands on an element boundary: real_len/nreal = 3.0
    (boundaries at 0, 3, 6, 9, 12); see module docstring."""
    els = (
        [ElementSpec(real_len / nreal)] * nreal
        + [ElementSpec(tail_len / ncomplex, theta_deg)] * ncomplex
    )
    return FemDvrEcsGrid(GridSpec(quadrature=nq, elements=els))


def test_B3_continuum_rotation():
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


def test_B4_bound_state_theta_independence():
    # Square well V=-V0 on [0,a], deep enough for a bound state; energy
    # invariant under theta. a=3.0 lands on an element boundary (see module
    # docstring) so the diagonal-potential DVR represents the discontinuity
    # cleanly.
    a, V0 = 3.0, 5.0

    def Vwell(z):
        return np.where(np.real(z) <= a, -V0, 0.0).astype(complex)

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
    assert abs(Eb[0] - Eb[1]) < 1e-8, Eb        # theta-independent
