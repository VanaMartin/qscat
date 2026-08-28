"""Spec gate 2: at s = 0 the coupled pole walk must reproduce the shipped
LCP curve. Both routes diagonalize the SAME 1-D Hamiltonian on the SAME
grids, so they should agree at eigenvalue round-off, not at a tolerance."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid
from qscat.core.lcp import local_complex_potential
from qscat.model import NO

from projects.no_coupled_channels.anisotropy import TwoCentreWell
from projects.no_coupled_channels.model import CoupledModel, DiagonalChannelModel
from validation.coupled.screen import ANGLES, NO_ELECTRONIC, coupled_resonance_curve
from validation.diatomic.config import CONFIGS

R_SAMPLE = np.linspace(2.0, 4.0, 9)
# The electronic Hamiltonian's diagonal CONTAINS v0(R), so the pole sits at
# v0(R) + eps_res, not at eps_res. A constant seed is adrift by up to 0.25 Ha
# -- five window half-widths -- and finds nothing.
SEED_OFFSET = 0.03 - 0.01j


def _seeds(R: np.ndarray) -> np.ndarray:
    """Window centres that track the neutral curve."""
    return np.asarray(NO.v0(R).real + SEED_OFFSET, dtype=np.complex128)


def _grids() -> tuple:
    return tuple(electronic_grid(angle_deg=a, **NO_ELECTRONIC) for a in ANGLES)


def _nuclear_nodes() -> np.ndarray:
    """A subsample of NO's own nuclear deck inside the well-resolved region.

    Evaluated AT grid nodes, never interpolated, so the comparison below has
    no interpolation error in it. Every 10th node keeps the test in the fast
    tier (~13 points instead of ~130).
    """
    R = np.asarray(CONFIGS["NO"].da_grid().grids[1].real_points, dtype=np.float64)
    inside = R[(R >= 2.0) & (R <= 4.0)]
    return inside[::10]


def test_s0_curve_reproduces_the_shipped_lcp_curve() -> None:
    """Spec gate 2. The coupled walk at s = 0 with one channel, against the
    SHIPPED `local_complex_potential` on `DiagonalChannelModel` -- a genuinely
    independent implementation of the same pole walk.

    Both diagonalize the same 1-D Hamiltonian on the same grids, so where they
    select the same state they agree at eigenvalue round-off. A looser
    agreement would be hiding a selection difference, not a tolerance.
    """
    ga, gb = _grids()
    well = TwoCentreWell(base=NO, s=0.0, kappa=0.3)
    nuclear = CONFIGS["NO"].da_grid().grids[1]
    vd_ref, gamma_ref = local_complex_potential(
        DiagonalChannelModel(well=well, l=1), nuclear, ga, gb
    )

    R_nodes = _nuclear_nodes()
    R_all = np.asarray(nuclear.real_points, dtype=np.float64)
    idx = np.searchsorted(R_all, R_nodes)

    curve = coupled_resonance_curve(
        CoupledModel(well=well, n_channels=1),
        R_nodes,
        ga,
        gb,
        seeds=np.asarray(vd_ref[idx].real - 0.5j * gamma_ref[idx], dtype=np.complex128),
    )
    np.testing.assert_allclose(curve.v_d, vd_ref[idx].real, atol=1e-9)
    np.testing.assert_allclose(curve.gamma, gamma_ref[idx], atol=1e-9)


def test_the_pole_is_insensitive_to_the_seed() -> None:
    """An angle-stable pole is determined by the grids, not by the window it
    was searched in. If this fails, the walk is selecting on the seed."""
    ga, gb = _grids()
    # s = 0.3: Gamma/eps is 0.56-0.94 across the sampled R, so the pole is
    # comfortably isolated and the residual is ~1e-9. At s = 1 it is 2.9 and
    # the state is no longer a resonance -- a poor place to test selection.
    model = CoupledModel(well=TwoCentreWell(base=NO, s=0.3, kappa=0.3), n_channels=2)
    base = coupled_resonance_curve(model, R_SAMPLE, ga, gb, seeds=_seeds(R_SAMPLE))
    shifted = coupled_resonance_curve(model, R_SAMPLE, ga, gb, seeds=base.E_res + (0.004 - 0.003j))
    np.testing.assert_allclose(shifted.E_res, base.E_res, atol=1e-9)


def test_n_poles_counts_only_the_residual_survivors() -> None:
    """`n_poles` is what a gate may use, so it needs its own check.

    `n_stable` counts every angle-stable state in the window and is 2 almost
    everywhere, because a spurious near-threshold state is always present.
    `n_poles` counts those that also pass the residual cut, so it can never
    exceed `n_stable`, and it must be at least 1 wherever a pole was actually
    recorded.
    """
    ga, gb = _grids()
    model = CoupledModel(well=TwoCentreWell(base=NO, s=0.3, kappa=0.3), n_channels=2)
    curve = coupled_resonance_curve(model, R_SAMPLE, ga, gb, seeds=_seeds(R_SAMPLE))
    assert curve.n_poles.shape == R_SAMPLE.shape
    assert np.all(curve.n_poles <= curve.n_stable)
    found = np.isfinite(curve.E_res)
    assert np.all(curve.n_poles[found] >= 1)
    assert np.all(curve.n_poles[~found] == 0)


def test_the_curve_records_how_many_stable_states_it_found() -> None:
    ga, gb = _grids()
    model = CoupledModel(well=TwoCentreWell(base=NO, s=0.3, kappa=0.3), n_channels=2)
    curve = coupled_resonance_curve(model, R_SAMPLE, ga, gb, seeds=_seeds(R_SAMPLE))
    assert curve.n_stable.shape == R_SAMPLE.shape
    assert np.all(curve.n_stable >= 1)


def test_kappa_zero_channel_parity() -> None:
    """PARITY ORACLE, at full anisotropy: the one identity gate in this
    campaign that is not evaluated at `s = 0`, where the coupling is
    switched off.

    At `kappa = 0` the two-centre well (`TwoCentreWell`) is left-right
    symmetric under `r -> -r` along the molecular axis, so only EVEN
    Legendre components of the interaction survive (see
    `anisotropy.v_lambda_closed_form`'s docstring: the `lam_A + (-1)^lambda
    lam_B` bracket vanishes for odd `lambda` when `lam_A = lam_B`). Within
    `Lambda = 1` that means `l = 1` couples only to `l = 3, 5, ...`, never to
    `l = 2` -- the `l = 1 <-> l = 2` block that a nonzero `kappa` is what
    opens (see the module docstring). Adding channel `l = 2` (`N_l = 2`)
    therefore changes NOTHING relative to the fixed-`l` model (`N_l = 1`):
    the off-diagonal `V_{12}` block is identically zero, at ANY `s`, not
    just `s = 0`.

    This is a genuine oracle rather than a restatement of the `s = 0`
    embedding gate elsewhere in this suite: it runs at `s = 0.3`, where the
    coupling itself is fully switched on, so it exercises the angular
    quadrature (`angular.theta_lm`), `_coupling_table`, and
    `assemble_coupled` under the exact conditions the campaign uses them --
    unlike every other identity check in this project, which runs where the
    coupling is off and so cannot see a bug in any of those three.
    """
    ga, gb = _grids()
    R = np.linspace(2.0, 3.6, 5)
    seeds = _seeds(R)
    one = coupled_resonance_curve(
        CoupledModel(well=TwoCentreWell(base=NO, s=0.3, kappa=0.0), n_channels=1),
        R,
        ga,
        gb,
        seeds=seeds,
    )
    two = coupled_resonance_curve(
        CoupledModel(well=TwoCentreWell(base=NO, s=0.3, kappa=0.0), n_channels=2),
        R,
        ga,
        gb,
        seeds=seeds,
    )
    np.testing.assert_allclose(two.E_res, one.E_res, atol=1e-9)


@pytest.mark.slow
def test_four_channels_run_at_the_production_deck_size() -> None:
    """The cost check the campaign is sized from: N_l = 4 on the real deck."""
    ga, gb = _grids()
    model = CoupledModel(well=TwoCentreWell(base=NO, s=0.3, kappa=0.3), n_channels=4)
    curve = coupled_resonance_curve(model, R_SAMPLE, ga, gb, seeds=_seeds(R_SAMPLE))
    assert np.all(np.isfinite(curve.E_res))
    assert np.all(curve.residual < 1e-3)
