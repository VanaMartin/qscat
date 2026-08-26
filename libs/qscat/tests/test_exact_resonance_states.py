"""Tests for `qscat.core.exact_resonance_states` -- the exact (non-BO) resonance
states of the 2-D model.

The oracle is the SEPARABLE LIMIT. Switch the electronic-nuclear coupling off and
the 2-D Hamiltonian becomes an exact Kronecker sum,

    H = (T_r + diag v_el) (+) (T_R + diag v_nuc),

so its eigenvalues are exactly pairwise sums `eps_el + eps_vib` of the 1-D
eigenvalues on the SAME grids, and its eigenvectors are exactly the products
`phi(r) chi(R)`. A 2-D resonance with a known position, a known width (that of
the electronic pole alone) and a known wavefunction -- which is what makes the
2-D search testable at all.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core import ExactResonanceStates, exact_resonance_states
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid, eigen, hamiltonian_nd, kinetic
from qscat.linalg import c_product
from qscat.model import N2

_R_FIXED = 2.02  # the electronic potential is frozen at this nuclear geometry
_MU = 12528.0  # N2 reduced mass, so the vibrational scale is realistic
_D_MORSE, _A_MORSE, _RE_MORSE = 0.30, 1.0, 2.02


def _v_el(r):
    """A genuine electronic resonance potential: the N2 surface at fixed R."""
    return np.asarray(N2.surface(r, _R_FIXED), dtype=np.complex128)


def _v_nuc(R):
    """A Morse well: real bound vibrational levels below a flat continuum."""
    x = 1.0 - np.exp(-_A_MORSE * (R - _RE_MORSE))
    return _D_MORSE * x * x - _D_MORSE


class _SeparableModel:
    """A deliberately SEPARABLE model: `V(r, R) = v_el(r) + v_nuc(R)`.

    Satisfies the `ResonanceModel` protocol as far as this solver uses it.
    """

    charge = 0
    mu = _MU

    def surface(self, r, R):
        return _v_el(r) + _v_nuc(R)

    def hamiltonian(self, tgrid):
        return hamiltonian_nd(tgrid, [1.0, self.mu], self.surface)


# Deliberately small (~11k unknowns, well under a second per 2-D solve) so the
# whole file stays in the fast suite. Converged physics needs a production deck;
# these tests check the machinery against an oracle evaluated on the SAME grid,
# which a coarse grid serves just as well.
#
# The ELECTRONIC grid is the exception and may not be shrunk freely: the
# acceptance criterion is an angle residual, so the grid must resolve the
# electronic pole well enough for that residual to fall under `rel_tol`. At
# order 5 / r_max 12 the 1-D pole is angle-stable only to 7e-4 and every
# resonance is rejected at the default 1e-4; at order 6 / n_complex 5 / r_max 14
# it is stable to 2e-7. The nuclear side has no such constraint here -- its
# factor is a bound state, stable to 1e-17.
def _elec(theta_r: float):
    return electronic_grid(r_max=14.0, order=6, n_complex=5, angle_deg=theta_r)


def _nuc(theta_R: float):
    return nuclear_grid(r_max=16.0, quadrature=6, n_complex=3, angle_deg=theta_R)


def _grids(theta_r: float = 35.0, theta_R: float = 25.0) -> TensorGrid:
    return TensorGrid([_elec(theta_r), _nuc(theta_R)])


def _three_grids() -> tuple[TensorGrid, TensorGrid, TensorGrid]:
    """Base, electronic-angle-moved, nuclear-angle-moved."""
    return _grids(35.0, 25.0), _grids(44.0, 25.0), _grids(35.0, 30.0)


def _one_d_reference() -> tuple[complex, np.ndarray, float, np.ndarray]:
    """The 1-D electronic pole and the lowest nuclear level, on the base grids.

    Computed densely, so the separable-limit expectation is independent of the
    sparse machinery under test.
    """
    from qscat.ecs import find_resonance_pole

    g_r_a, g_r_b = _elec(35.0), _elec(44.0)
    h_a = kinetic(g_r_a, 1.0) + np.diag(_v_el(g_r_a.points))
    h_b = kinetic(g_r_b, 1.0) + np.diag(_v_el(g_r_b.points))
    w_a, v_a = eigen(h_a)
    e_pole, _ = find_resonance_pole(w_a, eigen(h_b)[0], (-1.0, 0.0, -0.1, 0.0))
    phi = v_a[:, int(np.argmin(np.abs(w_a - e_pole)))]

    g_R = _nuc(25.0)
    w_R, v_R = eigen(kinetic(g_R, _MU) + np.diag(_v_nuc(g_R.points)))
    bound = np.flatnonzero((w_R.real < 0.0) & (np.abs(w_R.imag) < 1e-6))
    j = bound[int(np.argmin(w_R.real[bound]))]
    return complex(w_a[int(np.argmin(np.abs(w_a - e_pole)))]), phi, complex(w_R[j]), v_R[:, j]


_WINDOW = (-1.5, 0.5, -0.5, 0.0)
# Seed placement is not arbitrary: `k` eigenpairs nearest the shift is a LOCAL
# window, and the vibrational spacing here is omega ~ 0.0069 Ha. A seed 0.01
# above the target puts one and a half quanta between them, and v=0 is crowded
# out of the k nearest by v=1, v=2 -- the search then reports real states that
# are simply not the one asked for. Seed inside the level spacing.
_SEED_OFFSET = 0.002 - 0.001j


def _separable_search(shifts, k: int = 10):
    ga, gb, gc = _three_grids()
    return exact_resonance_states(_SeparableModel(), ga, gb, gc, shifts=shifts, k=k, window=_WINDOW)


def test_separable_limit_pole_is_the_sum_of_the_one_dimensional_eigenvalues() -> None:
    """The exact oracle: with no coupling, the 2-D pole IS eps_el + eps_vib.

    The tensor Hamiltonian is then an exact Kronecker sum of the two 1-D
    Hamiltonians ON THESE GRIDS, so this holds to solver precision and not
    merely to discretization accuracy -- a coarse grid does not weaken it.
    """
    e_el, _, e_vib, _ = _one_d_reference()
    expected = e_el + e_vib

    res = _separable_search([expected + _SEED_OFFSET])
    assert res.energies.size >= 1
    hit = res.energies[int(np.argmin(np.abs(res.energies - expected)))]
    assert abs(hit - expected) <= 1e-9 * abs(expected)


def test_separable_limit_width_is_the_electronic_width_alone() -> None:
    """The nuclear factor is bound, so it contributes no width."""
    e_el, _, e_vib, _ = _one_d_reference()
    expected = e_el + e_vib

    res = _separable_search([expected + _SEED_OFFSET])
    i = int(np.argmin(np.abs(res.energies - expected)))
    assert res.widths[i] == pytest.approx(-2.0 * e_el.imag, rel=1e-6)


def test_separable_limit_residuals_reduce_to_the_one_dimensional_ones() -> None:
    """The two residuals are separable too, and say different things.

    With no coupling the ELECTRONIC residual of the 2-D state must equal the
    1-D electronic pole's own angle residual -- the nuclear factor is identical
    on both grids and cancels -- while the NUCLEAR residual collapses to zero,
    because a bound vibrational state does not care where the nuclear contour
    turns. Anything else means the two-angle bookkeeping is wrong.
    """
    from qscat.ecs import find_resonance_pole

    e_el, _, e_vib, _ = _one_d_reference()
    g_a, g_b = _elec(35.0), _elec(44.0)
    _, residual_1d = find_resonance_pole(
        eigen(kinetic(g_a, 1.0) + np.diag(_v_el(g_a.points)))[0],
        eigen(kinetic(g_b, 1.0) + np.diag(_v_el(g_b.points)))[0],
        (-1.0, 0.0, -0.1, 0.0),
    )

    expected = e_el + e_vib
    res = _separable_search([expected + _SEED_OFFSET])
    i = int(np.argmin(np.abs(res.energies - expected)))
    assert res.residual_electronic[i] == pytest.approx(residual_1d, rel=1e-6)
    assert res.residual_nuclear[i] < 1e-12


def test_separable_limit_eigenvector_is_the_product_state() -> None:
    """phi(r) chi(R), up to the scale a non-Hermitian eigenvector carries."""
    e_el, phi, e_vib, chi = _one_d_reference()
    expected = e_el + e_vib

    res = _separable_search([expected + _SEED_OFFSET])
    i = int(np.argmin(np.abs(res.energies - expected)))
    psi = res.states[i]
    product = np.outer(phi, chi).ravel()
    product = product / np.sqrt(c_product(product, product))
    psi = psi / np.sqrt(c_product(psi, psi))
    assert abs(abs(c_product(psi, product)) - 1.0) < 1e-6


def test_states_are_row_per_state_and_load_rejects_column_caches(tmp_path) -> None:
    """lib-M9: one orientation convention (rows, like chi/phi); pre-flip
    caches must fail loudly, not load transposed."""
    e_el, _, e_vib, _ = _one_d_reference()
    expected = e_el + e_vib
    res = _separable_search([expected + _SEED_OFFSET])
    assert res.energies.size >= 1
    assert res.states.shape[0] == res.energies.size
    assert res.energies.size != res.states.shape[1]

    p = tmp_path / "states.npz"
    res.save(p)
    loaded = ExactResonanceStates.load(p)
    np.testing.assert_array_equal(loaded.states, res.states)

    # forge a pre-flip archive: transpose states
    from dataclasses import fields

    legacy = {f.name: getattr(res, f.name) for f in fields(res)}
    legacy["states"] = legacy["states"].T
    np.savez(tmp_path / "legacy.npz", **legacy)
    with pytest.raises(ValueError, match="column-per-state"):
        ExactResonanceStates.load(tmp_path / "legacy.npz")


@pytest.mark.slow
def test_electronic_continuum_is_rejected_by_the_electronic_angle_only() -> None:
    """The diagnostic the three-spectrum design buys.

    A state built on the ELECTRONIC continuum moves when theta_r changes but not
    when theta_R does, so the (A,B) pair rejects it while the (A,C) pair does
    not. Asserting this is what distinguishes 'no resonance here' from 'the
    search is broken'.
    """
    from qscat.core.resonance import _pooled_spectrum
    from qscat.ecs import match_angle_stable

    model = _SeparableModel()
    ga, gb, gc = _three_grids()
    # A shift well above the electronic threshold: continuum x bound-vibrational.
    sigma = -0.35 - 0.25j
    window = (-1.5, 0.5, -1.0, 0.0)
    va, _ = _pooled_spectrum(model, ga, [sigma], k=10)
    vb, _ = _pooled_spectrum(model, gb, [sigma], k=10)
    vc, _ = _pooled_spectrum(model, gc, [sigma], k=10)

    stable_ab, _, _ = match_angle_stable(va, vb, window)
    stable_ac, _, _ = match_angle_stable(va, vc, window)
    assert stable_ab.size == 0  # electronic angle moved it
    assert stable_ac.size > 0  # nuclear angle did not


def test_n2_yields_a_ladder_of_quasi_bound_states() -> None:
    """With the real coupling on, N2 gives a LADDER of angle-stable states.

    Not a single pole: a vibrational progression of the anion, each member
    stable under both ECS angles, with widths that grow as the state climbs --
    the boomerang picture, recovered without any Born-Oppenheimer or local
    approximation.

    Deliberately coarse grid, so the numbers here are not converged physics;
    what is asserted is the structure (a ladder at the vibrational scale, all
    members angle-stable, all widths positive), which the grid does not
    manufacture.
    """
    res = exact_resonance_states(
        N2,
        *_three_grids(),
        shifts=[-0.664 - 0.004j],
        k=12,
        window=(-0.75, -0.55, -0.05, 0.0),
    )
    assert res.energies.size >= 3
    assert np.all(res.widths > 0.0)
    assert np.all(res.residual_electronic < 1e-5)
    assert np.all(res.residual_nuclear < 1e-10)

    spacings = np.diff(res.energies.real)
    assert np.all(spacings > 0.005)  # a vibrational progression, not a cluster
    assert np.all(spacings < 0.015)


def test_electronic_continuum_states_do_not_reach_the_result() -> None:
    """The rejected states are in the raw spectrum and must not survive.

    On the same N2 search, the base spectrum contains states whose electronic
    factor rotates with theta_r (they sit at Im E ~ -0.008 and -0.020, moving
    2.5e-3 and 6.2e-3 when the electronic angle changes). None of them may
    appear in the accepted set -- that filtering IS the capability.
    """
    from qscat.core.resonance import _pooled_spectrum

    ga, gb, gc = _three_grids()
    shifts = [-0.664 - 0.004j]
    raw, _ = _pooled_spectrum(N2, ga, shifts, 12)
    res = exact_resonance_states(
        N2, ga, gb, gc, shifts=shifts, k=12, window=(-0.75, -0.55, -0.05, 0.0)
    )
    rejected = [e for e in raw if np.min(np.abs(res.energies - e)) > 1e-12]
    assert rejected  # the raw spectrum really does contain interlopers
    for e in rejected:
        assert np.min(np.abs(res.energies - e)) > 1e-6


# --- the grid-family guard ---------------------------------------------------
#
# `exact_resonance_states` accepted any three grids before this. Two failure
# modes were reachable and neither announced itself: a partner identical to the
# base makes every eigenvalue match itself with residual zero (so the search
# accepts the whole rotated continuum while reporting perfect stability), and a
# partner differing in the WRONG axis turns the residual into a discretization
# difference. Both are caught here, at the door.


def test_an_identical_partner_grid_is_rejected() -> None:
    """The dangerous one: it would 'pass' every state instead of failing loudly."""
    from qscat.exceptions import GridError

    ga, _, gc = _three_grids()
    with pytest.raises(GridError, match="identical to grid_base"):
        exact_resonance_states(
            _SeparableModel(), ga, ga, gc, shifts=[-0.6 - 0.001j], k=4, window=_WINDOW
        )
    ga, gb, _ = _three_grids()
    with pytest.raises(GridError, match="identical to grid_base"):
        exact_resonance_states(
            _SeparableModel(), ga, gb, ga, shifts=[-0.6 - 0.001j], k=4, window=_WINDOW
        )


def test_a_partner_that_moved_the_wrong_axis_is_rejected() -> None:
    """`grid_electronic` must move the ELECTRONIC angle, and only that."""
    from qscat.exceptions import GridError

    ga = _grids(35.0, 25.0)
    both_moved = _grids(44.0, 30.0)
    gc = _grids(35.0, 30.0)
    with pytest.raises(GridError, match="ELECTRONIC ECS angle only"):
        exact_resonance_states(
            _SeparableModel(), ga, both_moved, gc, shifts=[-0.6 - 0.001j], k=4, window=_WINDOW
        )


def test_a_partner_with_a_different_real_mesh_is_rejected() -> None:
    """Real-region error must cancel between the two spectra, or the test lies."""
    from qscat.core.grids import electronic_grid
    from qscat.exceptions import GridError

    ga = _grids(35.0, 25.0)
    coarser = TensorGrid(
        [electronic_grid(r_max=14.0, order=5, n_complex=5, angle_deg=44.0), _nuc(25.0)]
    )
    gc = _grids(35.0, 30.0)
    with pytest.raises(GridError, match="real nodes"):
        exact_resonance_states(
            _SeparableModel(), ga, coarser, gc, shifts=[-0.6 - 0.001j], k=4, window=_WINDOW
        )


def test_ecs_angle_family_builds_a_family_that_passes() -> None:
    """The builder exists so callers stop assembling this triple by hand."""
    from qscat.core.grids import ecs_angle_family

    base, moved_el, moved_nu = ecs_angle_family(
        _elec, _nuc, electronic_angles=(35.0, 44.0), nuclear_angles=(25.0, 30.0)
    )
    res = exact_resonance_states(
        _SeparableModel(), base, moved_el, moved_nu, shifts=[-0.664 - 0.004j], k=8, window=_WINDOW
    )
    assert res.energies.size >= 1


def test_ecs_angle_family_rejects_a_degenerate_angle_pair() -> None:
    from qscat.core.grids import ecs_angle_family
    from qscat.exceptions import GridError

    with pytest.raises(GridError, match="electronic_angles must differ"):
        ecs_angle_family(_elec, _nuc, electronic_angles=(35.0, 35.0), nuclear_angles=(25.0, 30.0))
    with pytest.raises(GridError, match="nuclear_angles must differ"):
        ecs_angle_family(_elec, _nuc, electronic_angles=(35.0, 44.0), nuclear_angles=(25.0, 25.0))


# --- persistence -------------------------------------------------------------


def test_save_load_round_trips_every_field(tmp_path) -> None:
    """Field names are the dataclass's business, not each call site's.

    Hand-rolled caches stored `res_el`/`res_nuc` where the dataclass says
    `residual_electronic`/`residual_nuclear` -- one rename from silently loading
    the wrong array into the wrong attribute.
    """
    from dataclasses import fields

    res = _separable_search([-0.664 - 0.004j], k=8)
    path = tmp_path / "poles.npz"
    res.save(path)
    back = ExactResonanceStates.load(path)
    for f in fields(ExactResonanceStates):
        assert np.array_equal(getattr(back, f.name), getattr(res, f.name)), f.name


def test_load_rejects_a_foreign_archive(tmp_path) -> None:
    path = tmp_path / "not-poles.npz"
    np.savez(path, energies=np.zeros(3, dtype=np.complex128))
    with pytest.raises(ValueError, match="missing"):
        ExactResonanceStates.load(path)
