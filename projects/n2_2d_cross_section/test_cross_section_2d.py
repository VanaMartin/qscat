"""The exact 2-D driven-equation VE cross section (sub-project #6, crux).

Validated WITHOUT reference data: the free-particle limit and the first Born
limit together pin the normalization, the ECS masking, the DVR coefficient
convention and the T-matrix; S-matrix reciprocity/unitarity (this file's
`test_s_matrix_is_reciprocal_and_conserves_flux`) is the test that genuinely
exercises `H_2D`, the LU-based propagator and the ECS boundary condition
together -- see that test's docstring for why the two Born-limit tests
cannot do this on their own (`T(lam) = lam*T1 + lam^2*T2` is exactly
quadratic in `lam_scale` regardless of whether `H_2D` is right or wrong).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import pytest
import scipy.sparse as sp
from qscat.dvr import TensorGrid
from qscat.linalg import SparseLU, c_product
from scipy.special import spherical_jn

from projects.n2_2d_cross_section.cross_section_2d import (
    channel_vector,
    ve_cross_section_2d,
)
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import ELL, MU, build_h2d, interaction_diag
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

# Small but physically sane: the interaction lives at r < ~3 bohr, so a modest
# box still supports a meaningful (if unconverged) T-matrix. Task 4 converges it.
TG = TensorGrid(
    [
        n2_electronic_grid(r_max=16.0, order=7, n_complex=5),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)


def _first_principles_channel_vector(
    tgrid: TensorGrid, k: float, chi_v: npt.NDArray[np.complex128], l: int
) -> npt.NDArray[np.complex128]:
    """`F_{E,l}(r) chi_v(R)`, built from scratch WITHOUT calling `channel_vector`
    (the solver's own constructor) or any other solver-side helper.

    This exists only to give `test_weak_coupling_matches_first_born` a
    reference that is actually independent of the solver. An earlier version
    of that test built its "independent" Born reference with `channel_vector`
    itself, which meant any error in `sqrt(w)`, the `sqrt(2k/pi)` energy
    normalization, the mask, the c-product convention, or the `4 pi^3/k^2`
    prefactor would appear identically in both `sigma` and `sigma_born` and
    cancel -- confirmed by monkeypatching `channel_vector` to scale its
    return by 1.5x and observing the old test's deviation was bit-identical
    with and without the patch (see review report). Every factor is written
    out explicitly here instead: the energy-normalized regular free radial
    solution `sqrt(2k/pi) r j_l(kr)` (`scipy.special.spherical_jn`, on REAL
    points only -- this function is never evaluated on the complex ECS tail),
    the grid's `sqrt(w)` DVR-coefficient factor, the c-product vibrational
    normalization, and the real-region mask.
    """
    g_r = tgrid.grids[0]
    r = g_r.real_points
    f_vals = np.sqrt(2.0 * k / np.pi) * r * spherical_jn(l, k * r)
    sqrt_w_r = tgrid.sqrt_weights()[0].ravel()
    f_coeff = f_vals * sqrt_w_r

    chi = np.asarray(chi_v, dtype=np.complex128)
    chi = chi / np.sqrt(c_product(chi, chi))

    psi = tgrid.outer([f_coeff, chi])
    psi[~tgrid.real_mask()] = 0.0
    return psi


def test_zero_driving_gives_exactly_zero() -> None:
    """`lam_scale=0` removes the driving/vertex term, so there is nothing to
    scatter off -- even though `H_2D` itself always carries the FULL
    interaction (`ve_cross_section_2d`'s docstring: `lam_scale` scales only
    `V_int`, never the propagator). Renamed from
    `test_free_particle_limit_gives_exactly_zero`, whose old name/docstring
    wrongly implied this exercised a free-particle Hamiltonian; it only ever
    tested zero driving.
    """
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [0, 1, 2], 0.2, lam_scale=0.0)
    assert np.all(sigma == 0.0)


def test_weak_coupling_matches_first_born() -> None:
    """As lam_scale -> 0, T -> <Phi_f|V_int|Psi_i>: the first Born amplitude.

    Two fixes relative to an earlier version of this test:

    1. The reference amplitude is built by `_first_principles_channel_vector`
       (see its docstring), never `channel_vector` -- so a bug in
       `channel_vector` is no longer guaranteed to cancel between `sigma`
       and `sigma_born`.
    2. The assertion compares the *ratio* `sigma[0] / sigma_born` against
       1.0 (an O(1) quantity) with `abs=0`, rather than
       `sigma[0] == pytest.approx(sigma_born, rel=1e-3)` directly. At this
       magnitude (`sigma`, `sigma_born` ~ 1e-17) the direct form is
       dominated by `pytest.approx`'s default `abs=1e-12` floor, not
       `rel=1e-3` -- true even for `sigma[0] == 0.0` (confirmed: with the
       old scale=1e-4, `0.0 == pytest.approx(7.716e-13, rel=1e-3)` is
       `True`). The ratio form has no such floor.

    `lam_scale=1e-6` (not the earlier 1e-4) is chosen for headroom: the
    residual second-Born contamination in the exact quadratic
    `T(lam) = lam*T1 + lam^2*T2` is O(lam_scale) and measured directly here
    at 4.95e-3 @ 1e-4, 4.93e-4 @ 1e-5, 4.93e-5 @ 1e-6. At 1e-4 the ratio
    fails the `rel=1e-3` gate outright (4.95e-3 > 1e-3) -- genuine
    second-Born contamination, not a tolerance artifact -- and even 1e-5
    leaves only ~2x headroom; 1e-6 gives ~20x headroom, matching
    `test_sigma_scales_as_lambda_squared_in_the_born_regime` below.
    """
    scale = 1e-6
    E, vp = 0.2, 1
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [vp], E, lam_scale=scale)

    # First Born, computed independently of ANY solver-side helper.
    k = np.sqrt(2.0 * E)
    e_tot = E + EPS[0]
    kp = np.sqrt(2.0 * (e_tot - EPS[vp]))
    psi_i = _first_principles_channel_vector(TG, k, CHI[0], ELL)
    phi_f = _first_principles_channel_vector(TG, kp, CHI[vp], ELL)
    v_int = scale * interaction_diag(TG)
    t_born = c_product(phi_f, v_int * psi_i)
    sigma_born = 4.0 * np.pi**3 * abs(t_born) ** 2 / (2.0 * E)

    assert sigma[0] / sigma_born == pytest.approx(1.0, rel=1e-3, abs=0.0)


def test_sigma_scales_as_lambda_squared_in_the_born_regime() -> None:
    """The ratio test compares an O(1) quantity, so pytest.approx's default
    abs=1e-12 tolerance does not mask residual second-Born contamination the
    way it could in a direct comparison (there, |T| is itself ~1e-13-1e-17,
    so an absolute floor can dominate regardless of lam_scale -- see
    `test_weak_coupling_matches_first_born`). Here the O(lam_scale) cross
    term in the exact quadratic `T(lam) = lam*T1 + lam^2*T2` (`H_2D` always
    carries the FULL, un-scaled interaction; only the driving/vertex
    `V_int` is lam_scale-scaled) is resolved directly.

    `lam_scale=1e-6` (moved from an earlier 1e-5, which had the ratio at
    4.0198 down to 4.0002 as lam_scale shrank -- see history in the review
    report -- leaving only ~2x headroom under `rel=1e-3`). Measured here:
    ratio = 4.000197 at lam_scale=1e-6 (rel_dev from 4.0 = 4.93e-5, ~20x
    headroom under `rel=1e-3`), matching the per-decade linear shrink of the
    second-Born contamination documented in
    `test_weak_coupling_matches_first_born`.
    """
    a = ve_cross_section_2d(TG, EPS, CHI, 0, [1], 0.2, lam_scale=1e-6)[0]
    b = ve_cross_section_2d(TG, EPS, CHI, 0, [1], 0.2, lam_scale=2e-6)[0]
    assert b / a == pytest.approx(4.0, rel=1e-3)


def test_s_matrix_is_reciprocal_and_conserves_flux() -> None:
    """`S = 1 - 2 pi i T` over the OPEN channels must be reciprocal
    (`S = S^T`, an exact algebraic consequence of `H_2D = H_2D^T` under ECS)
    and conserve flux (`S^dagger S = I`, exact for a closed multichannel
    system with no loss to channels outside the tracked set).

    Unlike the two Born-limit tests above, this genuinely exercises `H_2D`,
    the LU-based propagator, and the ECS boundary condition together: with
    `lam_scale` fixed at 1 (the real interaction), `T` here is the FULL
    non-perturbative T-matrix, not a series truncation, so a wrong `H_2D`
    (wrong mass on an axis, wrong sign of the kinetic term, a broken
    potential split), a wrong propagator (a transposed or mis-scaled LU
    solve), or a wrong boundary condition (wrong ECS angle/pivot, an
    unmasked channel function) generically breaks reciprocity and/or flux
    conservation -- whereas both Born tests above pass unconditionally
    regardless of whether `H_2D` is right, because `lam_scale` only scales
    the driving/vertex term, never the propagator (Critical-3 in the review
    that added this test).

    Energy chosen so that EXACTLY two vibrational channels are open (v=0,
    elastic, and v=1) and v=2 is the first closed channel: `e_tot` a hair
    below `EPS[2]`. This is deliberate, not incidental -- the E=0.2 collision
    energy used by the other tests in this file corresponds to
    `e_tot = E + EPS[0] ~= -0.545`, which (measured separately, using a
    15-state `vibrational_states` call) leaves at least 15 vibrational
    channels open, not the 4 tracked by this module's `EPS`/`CHI`. Building
    a "full open-channel" S-matrix out of only 4 of >=15 truly open channels
    measures real physical leakage into the untracked channels, not a
    solver defect: `|S^dagger S - I|` there is ~0.24, and -- tellingly --
    IDENTICAL to 3+ significant figures across several independent grid
    refinements (bigger electronic r_max, more DVR order, more ECS tail
    elements, finer nuclear quadrature), which rules out "just needs a
    bigger grid" and confirms it is a converged, genuine physical effect of
    channel truncation, not numerical error. Restricting to an energy where
    the tracked 4-state `EPS` truly is the complete open-channel set removes
    that confound.

    Measured on this module's (deliberately small, unconverged -- see the
    module-level `TG` comment) grid, at `e_tot=-0.727`:
    `max|S - S^T| = 3.25e-19` (reciprocity -- an exact algebraic identity of
    the complex-symmetric LU solve, at machine precision) and
    `max|S^dagger S - I| = 1.02e-6` (flux conservation, limited by the
    finite box / ECS-tail masking of the channel functions -- shrinks with
    grid refinement, e.g. 2.4e-7 measured with a larger electronic r_max).
    Tolerances below sit roughly an order of magnitude above each measured
    value.
    """
    e_tot = -0.727
    open_ = [v for v in range(len(EPS)) if e_tot - EPS[v] > 0.0]
    assert open_ == [0, 1]  # pins the intended open-channel count at this energy

    H = build_h2d(TG)
    v_diag = interaction_diag(TG)
    ident = sp.identity(TG.size, format="csc", dtype=np.complex128)
    lu = SparseLU((e_tot * ident - H).tocsc())

    phi = {}
    psi_plus = {}
    for v in open_:
        k = float(np.sqrt(2.0 * (e_tot - EPS[v])))
        pv = channel_vector(TG, k, CHI[v])
        phi[v] = pv
        psi_plus[v] = pv + lu.solve(v_diag * pv)

    n_open = len(open_)
    T = np.zeros((n_open, n_open), dtype=np.complex128)
    for a, vf in enumerate(open_):
        for b, vi in enumerate(open_):
            T[a, b] = c_product(phi[vf], v_diag * psi_plus[vi])

    S = np.eye(n_open, dtype=np.complex128) - 2j * np.pi * T
    reciprocity = np.max(np.abs(S - S.T))
    flux = np.max(np.abs(S.conj().T @ S - np.eye(n_open)))

    assert reciprocity < 1e-14  # measured 3.25e-19
    assert flux < 1e-5  # measured 1.02e-6


def test_sigma_is_real_and_non_negative() -> None:
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [0, 1, 2, 3], 0.2)
    assert sigma.dtype == np.float64
    assert np.all(sigma >= 0.0)


def test_closed_channels_are_zero() -> None:
    """At E below a channel's threshold that channel cannot be populated."""
    e_small = 0.005
    sigma = ve_cross_section_2d(TG, EPS, CHI, 0, [0, 1, 2, 3], e_small)
    open_ = (e_small + EPS[0] - EPS) > 0.0
    assert np.all(sigma[~open_[:4]] == 0.0)


def test_channel_vector_is_masked_to_the_unscaled_region() -> None:
    """A channel projection on the complex-scaled tail is meaningless."""
    psi = channel_vector(TG, 0.6, CHI[0])
    assert np.all(psi[~TG.real_mask()] == 0.0)
    assert np.abs(psi[TG.real_mask()]).max() > 0.0


@pytest.mark.slow
def test_array_of_energies_matches_scalar_calls() -> None:
    energies = [0.1, 0.2]
    both = ve_cross_section_2d(TG, EPS, CHI, 0, [1], energies)
    assert both.shape == (2, 1)
    for i, e in enumerate(energies):
        assert both[i, 0] == pytest.approx(
            ve_cross_section_2d(TG, EPS, CHI, 0, [1], e)[0], rel=1e-12
        )


@pytest.mark.slow
def test_reuse_swept_equals_per_energy_calls() -> None:
    """V2 gate: analyze-once/refactor-per-energy reuse must not change the
    physics. An array `E` builds ONE `SparseLU` and `refactor`s it per open
    energy; a scalar `E` builds a fresh `SparseLU`. The swept result must equal
    the stack of the individual scalar calls to round-off.

    The leading `E=0.0` (below threshold) exercises the lazy-init path: no
    factorization happens there, and the solver is built at the first `E > 0`
    (0.1), then refactored at 0.15 and 0.2. On the scipy path both routes run
    `splu` (bit-identical); on the MUMPS path the array route reuses the
    symbolic analysis while the scalar route re-analyzes each time, so they
    agree only to ~1e-9 -- hence the modest `rel=1e-9` tolerance rather than
    exact equality.
    """
    energies = [0.0, 0.1, 0.15, 0.2]
    vprimes = [0, 1, 2]
    swept = ve_cross_section_2d(TG, EPS, CHI, 0, vprimes, energies)
    assert swept.shape == (len(energies), len(vprimes))

    per_energy = np.stack([ve_cross_section_2d(TG, EPS, CHI, 0, vprimes, e) for e in energies])
    assert swept == pytest.approx(per_energy, rel=1e-9)
