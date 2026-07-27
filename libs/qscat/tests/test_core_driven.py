"""Gates for the promoted `qscat.core.channels`/`qscat.core.driven` (sub-project
#A, Task 4): the exact TI driven-Lippmann-Schwinger VE cross section, now
taking a `model: qscat.model.ResonanceModel` instead of a hardcoded N2
Hamiltonian.

Fixture mirrors `projects/n2_2d_cross_section/test_cross_section_2d.py`'s
small (deliberately unconverged) grid so the reference values below are on
a fast, real, already-exercised system.

Regression values (`_SIGMA_REF`/`_SIGMA_SCALAR_REF`) were captured from the
PRE-promotion `projects.n2_2d_cross_section.cross_section_2d.ve_cross_section_2d`
(checked out at the base commit, before this task rewired it into a shim over
this very module) -- comparing against the NOW-rewired shim would be
tautological, since the shim just calls straight through to
`qscat.core.driven.ve_cross_section`. Matched to round-off (`atol=0`,
`rtol=1e-12`): the promotion changes only where the Hamiltonian/interaction
diagonal/`l` come from (`model.hamiltonian`/`model.interaction_diag`/
`model.ell` instead of `build_h2d`/`interaction_diag`/`ELL`), not any
arithmetic, and `N2.hamiltonian`/`N2.interaction_diag` are themselves already
gated bit-identical to `build_h2d`/`interaction_diag`
(`libs/qscat/tests/test_model.py`).

Beyond the regression pin, this file gates the physics/contract invariants
per the task brief: sigma real & >=0, correct shape, the scalar/array +
`return_wavefunction` contract, and sigma>0 at the resonance (the cheap
"OR" alternative to a live Houfek-anchor comparison, which would need the
full converged `WORKING_GRID` -- out of scope for a fast unit test; the
Houfek gate itself lives in `validation/n2`'s Group E, run via the harness).
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.channels import channel_vector
from qscat.core.driven import ve_cross_section
from qscat.dvr import TensorGrid
from qscat.model import N2

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

TG = TensorGrid(
    [
        n2_electronic_grid(r_max=16.0, order=7, n_complex=5),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], N2.mu, 4)

# Captured from the pre-promotion `ve_cross_section_2d(TG, EPS, CHI, 0, [0, 1,
# 2], [0.1, 0.15, 0.2])` -- see module docstring. Row 2 (E=0.2, v'=1) also
# matches `projects/n2_2d_cross_section/convergence.py`'s independently
# documented WORKING_GRID value `sigma = 1.256450927036e-01`.
_E_ARR = np.array([0.1, 0.15, 0.2])
_SIGMA_REF = np.array(
    [
        [23.67447212318286, 6.122995193806422, 4.007759847511016],
        [7.238577710979102, 0.6257633940406637, 0.12794860997342516],
        [5.150658185623486, 0.12564509270361957, 0.012029795320066768],
    ]
)


def test_ve_cross_section_matches_pre_promotion_reference() -> None:
    sigma = ve_cross_section(TG, N2, EPS, CHI, 0, [0, 1, 2], _E_ARR)
    assert sigma.shape == (3, 3)
    np.testing.assert_allclose(sigma, _SIGMA_REF, rtol=1e-12, atol=0.0)


def test_scalar_energy_matches_the_corresponding_array_row() -> None:
    """Scalar `E` returns shape `(len(vprimes),)` and agrees with the array
    call's corresponding row -- the scalar/array contract `ve_cross_section_2d`
    documents, preserved verbatim by the promotion."""
    sigma_scalar = ve_cross_section(TG, N2, EPS, CHI, 0, [0, 1, 2], 0.2)
    assert sigma_scalar.shape == (3,)
    np.testing.assert_allclose(sigma_scalar, _SIGMA_REF[2], rtol=1e-12, atol=0.0)


def test_return_wavefunction_contract() -> None:
    """`return_wavefunction=True` also returns `psi_plus`: `None` below
    threshold, else a complex array of length `tgrid.size` for a scalar `E`."""
    sigma, psi = ve_cross_section(
        TG, N2, EPS, CHI, 0, [0, 1, 2], 0.2, return_wavefunction=True
    )
    np.testing.assert_allclose(sigma, _SIGMA_REF[2], rtol=1e-12, atol=0.0)
    assert isinstance(psi, np.ndarray)  # scalar E -> one array, never a list
    assert psi.shape == (TG.size,)
    assert psi.dtype == np.complex128

    sigma0, psi0 = ve_cross_section(
        TG, N2, EPS, CHI, 0, [0, 1, 2], 0.0, return_wavefunction=True
    )
    assert np.all(sigma0 == 0.0)
    assert psi0 is None


def test_sigma_is_real_and_non_negative() -> None:
    sigma = ve_cross_section(TG, N2, EPS, CHI, 0, [0, 1, 2, 3], 0.2)
    assert sigma.dtype == np.float64
    assert np.all(sigma >= 0.0)


def test_sigma_is_positive_at_the_resonance() -> None:
    """A cheap physics sanity gate (the brief's "OR" alternative to a live
    Houfek-anchor comparison, which needs the full converged grid): at E=0.2,
    v'=1 -- the resonance region this fixture is built around -- sigma must
    be strictly positive, not just non-negative."""
    sigma = ve_cross_section(TG, N2, EPS, CHI, 0, [1], 0.2)
    assert sigma[0] > 0.0


def test_closed_channels_are_zero() -> None:
    e_small = 0.005
    sigma = ve_cross_section(TG, N2, EPS, CHI, 0, [0, 1, 2, 3], e_small)
    open_ = (e_small + EPS[0] - EPS) > 0.0
    assert np.all(sigma[~open_[:4]] == 0.0)


def test_array_of_energies_matches_scalar_calls() -> None:
    energies = [0.1, 0.2]
    both = ve_cross_section(TG, N2, EPS, CHI, 0, [1], energies)
    assert both.shape == (2, 1)
    for i, e in enumerate(energies):
        assert both[i, 0] == pytest.approx(
            ve_cross_section(TG, N2, EPS, CHI, 0, [1], e)[0], rel=1e-12
        )


def test_channel_vector_matches_pre_promotion_reference() -> None:
    """`l` is now a required (positional) parameter -- `N2.ell` supplies the
    value the old default `ELL` used to. Value/shape/dtype captured from the
    pre-promotion `channel_vector(TG, 0.6, CHI[0])` (default `l=ELL=2`)."""
    cv = channel_vector(TG, 0.6, CHI[0], N2.ell)
    assert cv.shape == (TG.size,)
    assert cv.dtype == np.complex128
    assert np.abs(cv).sum() == pytest.approx(60.17329455106118, rel=1e-12)


def test_channel_vector_is_masked_to_the_unscaled_region() -> None:
    psi = channel_vector(TG, 0.6, CHI[0], N2.ell)
    assert np.all(psi[~TG.real_mask()] == 0.0)
    assert np.abs(psi[TG.real_mask()]).max() > 0.0
