"""Independent polynomial oracle recovery and isolation of a dense unused audit."""

import numpy as np

from projects.ukrmol_co.phase_identifiability import audit, candidate_s, fit_background, oracle_s


def test_fixed_true_pole_recovers_nonconstant_polynomial_s_factor():
    energy = np.linspace(0.04, 0.14, 61)
    case = {
        "position_hartree": 0.09,
        "full_width_hartree": 0.01,
        "background": [0.1, -0.1, 0.03],
        "exponential_background_amplitude": 0,
    }
    values = oracle_s(energy, case)
    fit = fit_background(energy, values, [0.04, 0.14], [0.09, 0.01], 3)
    dense = np.linspace(0.04, 0.14, 501)
    np.testing.assert_allclose(candidate_s(dense, fit), oracle_s(dense, case), atol=1e-13, rtol=0)
    np.testing.assert_allclose(abs(candidate_s(dense, fit)), 1, atol=3e-15, rtol=0)
    assert fit["full_rank_and_finite"]


def test_narrow_window_can_hide_known_wrong_width_on_unused_points():
    energy = np.linspace(0.13, 0.15, 61)
    case = {
        "position_hartree": 0.14,
        "full_width_hartree": 0.085,
        "background": [0.1, -0.1],
        "exponential_background_amplitude": 0.03,
    }
    fit = fit_background(energy, oracle_s(energy, case), [0.13, 0.15], [0.14, 0.0935], 7)
    unused = np.linspace(0.131, 0.149, 1001)
    difference = np.angle(candidate_s(unused, fit) * oracle_s(unused, case).conj()) / 2
    assert max(abs(difference)) < 0.001
    assert (
        abs(fit["full_width_hartree"] - case["full_width_hartree"])
        > 0.05 * case["full_width_hartree"]
    )


def test_inner_window_endpoints_are_unused_even_when_original_indices_are_development():
    import json
    from pathlib import Path

    contract = json.loads(
        Path(__file__).with_name("phase-identifiability-contract.json").read_text()
    )
    # Original 61-point grid indices 3 and 57 lie at these inner-window endpoints.
    # Both would be development points under the every-fifth-point rule alone.
    contract.update(grid_points=[61], background_degrees=[3], exact_pole_width_factors=[1])
    contract["cases"] = [contract["cases"][0] | {"windows_hartree": [[0.045, 0.135]]}]
    result = audit(contract)
    grid = result["grids"][0]
    develop = np.array(grid["development_energy_hartree"])
    for key in ["endpoint_energy_hartree", "reserved_energy_hartree", "dense_audit_energy_hartree"]:
        assert np.min(abs(develop[:, None] - np.array(grid[key]))) > 9e-14
