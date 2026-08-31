"""Channel vectors in a channel-outermost layout."""

from __future__ import annotations

import numpy as np
from qscat.core.channels import channel_vector
from qscat.core.driven import ve_cross_section
from qscat.core.grids import segmented_grid
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel, DiagonalChannelModel
from projects.no_coupled_channels.scattering import (
    channel_block,
    coupled_channel_vector,
    coupled_ve_cross_section,
)


def _tensor_grid() -> TensorGrid:
    """A deliberately small 2-D grid: these tests are about layout, not physics."""
    el = segmented_grid(((4, 8.0),), ((2, 12.0),), angle_deg=35.0, quadrature=6)
    nu = segmented_grid(((3, 4.0),), ((2, 6.0),), angle_deg=30.0, quadrature=6, x_min=1.0)
    return TensorGrid([el, nu])


def test_the_vector_lands_in_its_own_block_and_nowhere_else() -> None:
    tg = _tensor_grid()
    n = tg.size
    _eps, chi = vibrational_states(tg.grids[1], NO.mu, 2, NO.v0)
    ells = (1, 2, 3)
    for c in range(3):
        vec = coupled_channel_vector(tg, 0.3, chi[0], ells, c)
        assert vec.shape == (3 * n,)
        np.testing.assert_allclose(
            channel_block(vec, c, n), channel_vector(tg, 0.3, chi[0], ells[c])
        )
        for other in range(3):
            if other != c:
                assert np.all(channel_block(vec, other, n) == 0.0)


def test_one_channel_reproduces_the_shipped_vector_exactly() -> None:
    """With a single channel the coupled vector IS the shipped one -- the
    layout must add nothing at n_channels = 1."""
    tg = _tensor_grid()
    _eps, chi = vibrational_states(tg.grids[1], NO.mu, 2, NO.v0)
    got = coupled_channel_vector(tg, 0.3, chi[0], (1,), 0)
    np.testing.assert_array_equal(got, channel_vector(tg, 0.3, chi[0], 1))


def test_the_partial_wave_selects_the_bessel_order() -> None:
    """Same k, different l -- the momentum is shared, the order is not."""
    tg = _tensor_grid()
    n = tg.size
    _eps, chi = vibrational_states(tg.grids[1], NO.mu, 2, NO.v0)
    ells = (1, 2)
    a = channel_block(coupled_channel_vector(tg, 0.3, chi[0], ells, 0), 0, n)
    b = channel_block(coupled_channel_vector(tg, 0.3, chi[0], ells, 1), 1, n)
    assert float(np.max(np.abs(a - b))) > 1e-6


E_TEST = np.array([0.02, 0.05])
VPRIMES = [0, 1]


def _basis(tg: TensorGrid):
    return vibrational_states(tg.grids[1], NO.mu, 3, NO.v0)


def test_s0_reproduces_the_shipped_solver() -> None:
    """The duplicated sweep must give the shipped answer where the models are
    the same Hamiltonian. This is what certifies the duplication."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)

    got = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    want = ve_cross_section(tg, DiagonalChannelModel(well=well, l=1), eps, chi, 0, VPRIMES, E_TEST)
    np.testing.assert_allclose(got.total, want, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(got.restricted, want, rtol=1e-12, atol=0.0)


def test_s0_many_channels_equal_one_channel() -> None:
    """The embedding identity, carried all the way to an observable.

    At s = 0 the coupled Hamiltonian is block-diagonal, so an electron
    entering channel 0 can never reach another channel: a four-channel run
    must give exactly the one-channel answer, and its total must equal its
    restricted part because no other exit is reachable. The preceding phase
    gated this at the Hamiltonian; this gates it through the solve, the
    projection and the normalisation as well.
    """
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.5)
    one = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    four = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=4), eps, chi, 0, VPRIMES, E_TEST
    )
    np.testing.assert_allclose(four.total, one.total, rtol=1e-10, atol=0.0)
    np.testing.assert_allclose(four.total, four.restricted, rtol=1e-10, atol=0.0)


def test_kappa_zero_two_channels_equal_one() -> None:
    """Parity: only even Legendre components survive a symmetric well, so
    l = 1 cannot reach l = 2 at ANY anisotropy. An identity, at s well away
    from zero where the coupling is otherwise fully on."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.0)
    one = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    two = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=2), eps, chi, 0, VPRIMES, E_TEST
    )
    np.testing.assert_allclose(two.total, one.total, rtol=1e-9, atol=0.0)
    np.testing.assert_allclose(two.restricted, one.restricted, rtol=1e-9, atol=0.0)


def test_the_coupling_moves_the_cross_section() -> None:
    """With kappa on, l = 1 reaches l = 2 and the answer must change --
    otherwise the previous test proves nothing."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.5)
    one = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=1), eps, chi, 0, VPRIMES, E_TEST
    )
    two = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=2), eps, chi, 0, VPRIMES, E_TEST
    )
    rel = np.abs(two.total - one.total) / np.maximum(one.total, 1e-30)
    assert float(np.max(rel)) > 1e-3


def test_total_is_at_least_the_restricted_part() -> None:
    """The total sums over exit channels and the restricted one is a single
    term of that sum, so the total can never be smaller."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.6, kappa=0.5)
    out = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=3), eps, chi, 0, VPRIMES, E_TEST
    )
    assert np.all(out.total >= out.restricted - 1e-15)


def test_a_closed_channel_contributes_nothing() -> None:
    """Below its threshold an exit channel is closed and must be zero, not a
    small number from an imaginary momentum."""
    tg = _tensor_grid()
    eps, chi = _basis(tg)
    well = TwoCentreWell(base=NO, s=0.3, kappa=0.5)
    tiny = np.array([1e-4])  # far below the v' = 1 threshold
    out = coupled_ve_cross_section(
        tg, CoupledModel(well=well, n_channels=2), eps, chi, 0, VPRIMES, tiny
    )
    assert out.total[0, 1] == 0.0
