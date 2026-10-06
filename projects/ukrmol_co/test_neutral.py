"""A correlated neutral energy must retain its reference/density qualification."""

import json

import pytest

from projects.ukrmol_co.analyze import neutral_properties


@pytest.mark.parametrize("defect", [None, "unstable", "lambda", "density", "triples", "reference"])
def test_neutral_pilot_rejects_inapplicable_or_inconsistent_records(tmp_path, defect):
    report = {
        "basis": "cc-pVDZ",
        "bond_length_bohr": 2.1323,
        "rhf_converged": True,
        "rhf_internally_stable": True,
        "rhf_externally_stable": True,
        "ccsd_converged": True,
        "lambda_converged": True,
        "frozen_spatial_orbitals": 2,
        "correlated_electrons": 10,
        "rhf_energy_hartree": -112.75,
        "ccsd_energy_hartree": -113.0,
        "triples_correction_hartree": -0.01,
        "neutral_energy_hartree": -113.01,
        "density_electrons": 14,
        "t1_frobenius_over_sqrt_correlated_electrons": 0.02,
        "t1_largest_singular_value": 0.04,
        "ccsd_dipole_au": [0, 0, 0.04],
        "independent_reference": {"rhf_energy_hartree": -112.75, "neutral_energy_hartree": -113.01},
    }
    if defect == "unstable":
        report["rhf_externally_stable"] = False
    elif defect == "lambda":
        report["lambda_converged"] = False
    elif defect == "density":
        report["density_electrons"] = 10
    elif defect == "triples":
        report["neutral_energy_hartree"] = -113.0
    elif defect == "reference":
        report["independent_reference"]["neutral_energy_hartree"] += 1e-5
    config = {"reference_check": True, "basis": "cc-pVDZ", "bond_length": 2.1323}
    (tmp_path / "neutral.json").write_text(json.dumps(report))
    if defect:
        with pytest.raises((ValueError, AssertionError)):
            neutral_properties(tmp_path, config)
    else:
        result = neutral_properties(tmp_path, config)
        assert result["neutral_only"]
        assert (
            "basis and single-reference convergence not established" in result["validation_scope"]
        )
