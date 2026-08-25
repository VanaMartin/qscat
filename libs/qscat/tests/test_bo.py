"""Tests for `qscat.core.bo` -- the Born-Oppenheimer reference states.

The oracle is again the SEPARABLE LIMIT, and it is unusually sharp here. With
`V(r, R) = v_el(r) + v_nuc(R)` the frozen-nucleus electronic problem at `R` is

    -1/2 d^2/dr^2 + v_el(r) + v_nuc(R),

whose eigenvalues are `eps_el_j + v_nuc(R)` EXACTLY and whose eigenvectors are
`phi_j(r)`, independent of `R`. So every quantity this module produces has a
closed form on the same grid:

- `electronic_curves(...).energies[j]` is `eps_el_j + v_nuc(R)`, to solver
  precision rather than discretization accuracy;
- every column of `states[j]` is the same `phi_j` (which also makes phase
  alignment testable: without it the columns differ by random phases);
- `bo_basis`'s levels are `eps_el_j + (levels of v_nuc)`;
- its product states are `phi_j (x) chi_v` exactly.

A coarse grid does not weaken any of this, because the oracle is evaluated on
the identical discretization.
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest
from qscat.core import anion_electronic_states, vibrational_states
from qscat.core.bo import (
    BoBasis,
    admissible_levels,
    basis_covers,
    bo_basis,
    electronic_curves,
    n_eff,
    resonance_curve,
)
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import eigen, kinetic
from qscat.exceptions import GridError
from qscat.linalg import c_product
from qscat.model import N2

_R_FIXED = 2.02
_MU = 12528.0
_D_MORSE, _A_MORSE, _RE_MORSE = 0.30, 1.0, 2.02


# A BINDING electronic well. It has to bind: `electronic_curves` selects the
# lowest-Re(E) eigenvalues, and on an ECS grid a potential with no bound state
# offers only rotated continuum there -- `bo_basis` then (correctly) refuses to
# build vibrational levels in a complex curve. N2's own surface at fixed R is
# such a potential: its resonance is a SHAPE resonance with nothing bound below
# it, which is why it serves `resonance_curve` below and not this fixture.
# Depth and width are chosen to bind five states, so `n_curves > 1` is testable.
def _v_el(r):
    x = np.asarray(r, dtype=np.complex128)
    return np.asarray(-5.0 * np.exp(-0.05 * x * x), dtype=np.complex128)


def _v_nuc(R):
    x = 1.0 - np.exp(-_A_MORSE * (np.asarray(R) - _RE_MORSE))
    return _D_MORSE * x * x - _D_MORSE


class _SeparableModel:
    """`V(r, R) = v_el(r) + v_nuc(R)` -- the analytic oracle."""

    charge = 0
    mu = _MU
    ell = 0

    def surface(self, r, R):
        return _v_el(r) + _v_nuc(R)

    def v0(self, R):
        return np.asarray(_v_nuc(R), dtype=np.complex128)


# Small enough to stay in the fast suite; the oracle is exact on this grid.
def _elec(theta: float = 35.0):
    return electronic_grid(r_max=14.0, order=6, n_complex=5, angle_deg=theta)


def _nuc(theta: float = 25.0):
    return nuclear_grid(r_max=16.0, quadrature=6, n_complex=3, angle_deg=theta)


@pytest.fixture(scope="module")
def grids():
    return _elec(), _nuc()


@pytest.fixture(scope="module")
def electronic_reference(grids):
    """`(eps_el, phi)` of the frozen electronic problem, computed densely."""
    g_r, _ = grids
    w, v = eigen(kinetic(g_r, 1.0) + np.diag(_v_el(g_r.points)))
    return w, v


# --- electronic_curves -------------------------------------------------------


def test_curves_are_the_analytic_sum(grids, electronic_reference):
    """`energies[j, k] == eps_el_j + v_nuc(R_k)` -- exact, not approximate."""
    g_r, g_R = grids
    eps_el, _ = electronic_reference
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=4)
    expected = eps_el[:4, None] + _v_nuc(g_R.points)[None, :]
    assert np.allclose(cur.energies, expected, rtol=0, atol=1e-10)


def test_curves_without_states_carry_no_states(grids):
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=2)
    assert not cur.has_states
    assert cur.states.size == 0
    assert cur.n_curves == 2


def test_phase_alignment_makes_the_columns_identical(grids):
    """Every column of a curve is the SAME phi_j here, so alignment shows up.

    The eigenvector is `R`-independent in the separable limit, but its PHASE is
    not: `eigen` fixes it arbitrarily at each `R`. Alignment by continuity is
    what makes the columns agree, and without it this assertion fails on sign
    flips alone -- which is exactly the failure mode that silently zeroes an
    overlap against a smooth partner.
    """
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=3, with_states=True)
    assert cur.has_states
    for j in range(3):
        first = cur.states[j, :, 0]
        for k in range(1, g_R.n):
            assert np.allclose(cur.states[j, :, k], first, rtol=0, atol=1e-10)


def test_curves_reject_impossible_counts(grids):
    g_r, g_R = grids
    with pytest.raises(GridError):
        electronic_curves(_SeparableModel(), g_r, g_R, n_curves=0)
    with pytest.raises(GridError):
        electronic_curves(_SeparableModel(), g_r, g_R, n_curves=g_r.n + 1)


# --- bo_basis ----------------------------------------------------------------


def test_levels_are_the_electronic_energy_plus_the_vibrational_ladder(grids, electronic_reference):
    """`energies[j, v] == eps_el_j + eps_vib_v`, the separable-limit oracle."""
    g_r, g_R = grids
    eps_el, _ = electronic_reference
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=2)
    basis = bo_basis(cur, g_R, _MU, n_vib=3)

    vib = vibrational_states(g_R, _MU, 3, lambda R: np.asarray(_v_nuc(R), dtype=np.complex128))
    for j in range(2):
        expected = float(eps_el[j].real) + np.asarray(vib.eps, dtype=np.float64)
        assert np.allclose(basis.energies[j], expected, rtol=0, atol=1e-9)


def test_product_states_are_the_kronecker_product(grids):
    """`psi_(j,v) == phi_j (x) chi_v`, up to the c-product normalization."""
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=2, with_states=True)
    basis = bo_basis(cur, g_R, _MU, n_vib=2)

    vib = vibrational_states(g_R, _MU, 2, lambda R: np.asarray(_v_nuc(R), dtype=np.complex128))
    for (j, v), state in basis.items():
        raw = np.outer(cur.states[j, :, 0], vib.chi[v]).ravel()
        raw = raw / np.sqrt(complex(c_product(raw, raw)))
        # Sign is not fixed by the c-product normalization; compare magnitudes
        # of the overlap, which is the quantity anything downstream uses.
        ov = abs(complex(c_product(raw, state.psi)))
        assert ov == pytest.approx(1.0, abs=1e-8)
        assert state.curve == j and state.vib == v


def test_product_states_are_c_normalized(grids):
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=2, with_states=True)
    for _, state in bo_basis(cur, g_R, _MU, n_vib=2).items():
        assert complex(c_product(state.psi, state.psi)) == pytest.approx(1.0, abs=1e-10)


def test_basis_without_states_still_gives_the_level_table(grids):
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=2)
    basis = bo_basis(cur, g_R, _MU, n_vib=2)
    assert not basis.has_states
    assert np.isfinite(basis.energies).all()
    assert basis.levels() == [(0, 0), (0, 1), (1, 0), (1, 1)]


def test_flat_orders_levels_by_energy(grids):
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=2)
    e, keys = bo_basis(cur, g_R, _MU, n_vib=3).flat()
    assert len(keys) == e.size
    assert np.all(np.diff(e) >= 0)


def test_bo_basis_rejects_a_grid_mismatch(grids):
    """Curves tabulated on one grid may not be laid onto another."""
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=1)
    other = nuclear_grid(r_max=16.0, quadrature=6, n_complex=5, angle_deg=25.0)
    assert other.n != g_R.n
    with pytest.raises(GridError, match="same grid"):
        bo_basis(cur, other, _MU, n_vib=1)
    with pytest.raises(GridError):
        bo_basis(cur, g_R, _MU, n_vib=0)


def test_allow_partial_pads_a_curve_that_runs_out(grids):
    """A request past a curve's capacity raises strictly and pads permissively.

    Curves do not share one vibrational capacity, so `n_vib` is a request. The
    strict default is what tells a caller its rectangular assumption is wrong.
    """
    g_r, g_R = grids
    cur = electronic_curves(_SeparableModel(), g_r, g_R, n_curves=1)
    with pytest.raises(GridError):
        bo_basis(cur, g_R, _MU, n_vib=400)
    padded = bo_basis(cur, g_R, _MU, n_vib=400, allow_partial=True)
    finite = np.isfinite(padded.energies[0])
    assert finite.any(), "allow_partial should still find the levels that exist"
    assert not finite.all(), "this curve cannot support 400 clean levels"
    # Padding is a suffix: the levels found are the LOWEST ones, contiguous.
    assert np.all(finite[: int(finite.sum())])


# --- the Rydberg / closed-channel arithmetic ---------------------------------


def test_n_eff_uses_the_nearest_threshold_above():
    """Binding is measured to the threshold the state sits BELOW, not the lowest.

    A Rydberg state is bound against the ion level it belongs to. Using the
    lowest threshold instead would report a far larger binding, hence a far
    smaller `n_eff`, for every state above the first channel.
    """
    thresholds = [-0.10, -0.08, -0.06]
    # binding to -0.08 is 0.02 -> n_eff = 1/sqrt(0.04) = 5
    assert n_eff(-0.10, thresholds) == pytest.approx(5.0)
    # the -0.10 threshold is below and must be ignored
    assert n_eff(-0.09, thresholds) == pytest.approx(1.0 / np.sqrt(0.02))


def test_n_eff_above_every_threshold_is_undefined():
    with pytest.raises(ValueError, match="above every threshold"):
        n_eff(0.0, [-0.10, -0.08])


def test_admissible_levels_shrink_the_index_as_the_vibrational_level_rises():
    """The closed-channel constraint: higher `v` implies LOWER Rydberg index.

    At fixed energy a higher vibrational threshold means a larger binding, so a
    smaller `n_eff`. This is what makes the admissible set finite -- and it is
    what separates a spurious pole from one whose partner was never built.
    """
    thresholds = np.array([-0.10, -0.08, -0.06, -0.04])
    got = admissible_levels(-0.105, thresholds)
    curves = [j for j, _ in got]
    vibs = [v for _, v in got]
    assert vibs == sorted(vibs)
    assert curves == sorted(curves, reverse=True), "index must fall as v rises"


def test_admissible_levels_ignore_open_channels():
    """Only thresholds ABOVE the state contribute -- an open channel is continuum."""
    thresholds = np.array([-0.10, -0.08, -0.06])
    got = admissible_levels(-0.07, thresholds)
    assert [v for _, v in got] == [2]


def test_admissible_levels_can_be_cut_at_an_accumulation_region():
    thresholds = np.array([-0.10, -0.08])
    just_below = -0.080001  # binding 1e-6 -> n_eff ~ 707
    assert admissible_levels(just_below, thresholds) != []
    assert admissible_levels(just_below, thresholds, n_eff_max=12.0) == []


def test_basis_covers_is_false_when_an_admissible_level_is_missing():
    """The distinction the overlap test cannot make on its own."""
    thresholds = np.array([-0.10, -0.08, -0.06])
    energies = np.array([[-0.5]], dtype=np.float64)
    empty = BoBasis(energies=energies, states={})
    assert not basis_covers(-0.105, thresholds, empty)

    admissible = admissible_levels(-0.105, thresholds)
    from qscat.core.bo import BoState

    full = BoBasis(
        energies=energies,
        states={k: BoState(np.zeros(1, dtype=np.complex128), 0.0, k[0], k[1]) for k in admissible},
    )
    assert basis_covers(-0.105, thresholds, full)


def test_basis_covers_accepts_a_neighbouring_curve_index():
    """`n_eff ~ j+1` is the asymptotic relation; low curves depart from it."""
    from qscat.core.bo import BoState

    thresholds = np.array([-0.10, -0.08])
    (j, v) = admissible_levels(-0.105, thresholds)[0]
    off_by_one = BoBasis(
        energies=np.array([[-0.5]], dtype=np.float64),
        states={
            (jj + 1, vv): BoState(np.zeros(1, dtype=np.complex128), 0.0, jj + 1, vv)
            for jj, vv in admissible_levels(-0.105, thresholds)
        },
    )
    assert (j, v) not in off_by_one
    assert basis_covers(-0.105, thresholds, off_by_one)
    assert not basis_covers(-0.105, thresholds, off_by_one, curve_tol=0)


# --- resonance_curve ---------------------------------------------------------


def test_resonance_curve_reproduces_the_local_complex_potential(grids):
    """The curve energy IS `V_d - i Gamma/2` -- the same object `lcp` computes.

    `resonance_curve` and `qscat.core.lcp.local_complex_potential` run the same
    two-angle pole walk; this one keeps the eigenvectors the other discards, so
    they must agree on what they both compute.
    """
    from qscat.core.lcp import local_complex_potential

    g_R = nuclear_grid(r_max=16.0, quadrature=6, n_complex=3, angle_deg=25.0)
    ea, eb = _elec(35.0), _elec(44.0)
    Vd, Gamma = local_complex_potential(N2, g_R, ea, eb)

    eps_e, _ = anion_electronic_states(ea, N2, g_R.R0, 1)
    seed = (eps_e[0] - 0.05, eps_e[0] + 0.05, -0.05, 0.05)
    cur = resonance_curve(N2, ea, eb, g_R, seed, with_states=True)

    real = np.flatnonzero(g_R.points.imag == 0.0)
    assert np.allclose(cur.energies[0, real].real, Vd[real].real, rtol=0, atol=1e-12)
    # `local_complex_potential` CLAMPS the width at zero because a negative
    # Gamma is not a width; `resonance_curve` keeps the raw pole, because a
    # tiny positive Im is honest numerical noise and clamping it would hide
    # which nuclear geometries the walk is struggling at. Compare like for like.
    assert np.allclose(
        np.maximum(0.0, -2.0 * cur.energies[0, real].imag), Gamma[real], rtol=0, atol=1e-10
    )
    assert cur.n_curves == 1 and cur.has_states


def test_resonance_curve_states_are_phase_aligned_and_frozen_on_the_tail(grids):
    """Adjacent columns must not flip sign, and the tail takes the outer state."""
    g_R = nuclear_grid(r_max=16.0, quadrature=6, n_complex=3, angle_deg=25.0)
    ea, eb = _elec(35.0), _elec(44.0)
    eps_e, _ = anion_electronic_states(ea, N2, g_R.R0, 1)
    cur = resonance_curve(N2, ea, eb, g_R, (eps_e[0] - 0.05, eps_e[0] + 0.05, -0.05, 0.05))

    pts = g_R.points
    real = np.flatnonzero(pts.imag == 0.0)
    order = real[np.argsort(pts[real].real)]
    for a, b in pairwise(order):
        ov = complex(np.vdot(cur.states[0, :, a], cur.states[0, :, b]))
        assert ov.real > 0, "adjacent columns flipped phase"

    tail = np.flatnonzero(pts.imag != 0.0)
    if tail.size:
        outer = order[-1]
        for k in tail:
            assert np.allclose(cur.states[0, :, k], cur.states[0, :, outer])
