"""Independent analytic oracle and metrics; never imports an ODE solver."""

import numpy as np

from projects.ukrmol_co.verify_coupled_jost import fundamental_maps, measurements


def test_zero_coupling_reduces_to_canonical_free_waves():
    y, a = fundamental_maps(1, 0)
    np.testing.assert_allclose(y, np.diag([1j, -1]), atol=1e-14, rtol=0)
    np.testing.assert_allclose(a, [[np.exp(-1j), 0]], atol=1e-14, rtol=0)


def test_coupled_flux_and_wronskian_are_independent_checks():
    y, a = fundamental_maps(1, 0.4)
    correct = measurements(y, a, (y, a))
    np.testing.assert_allclose(
        [correct["flux_defect"], correct["wronskian_residual"]], 0, atol=1e-13, rtol=0
    )
    assert measurements(y, 2 * a, (y, a))["flux_defect"] > 1
