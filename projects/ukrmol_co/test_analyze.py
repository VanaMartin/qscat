"""Output-format regression checks; numerical engine validation uses saved runs."""

import hashlib
import json

import numpy as np
import pytest

from projects.ukrmol_co.analyze import ground_cross_sections, target_properties


def test_ground_cross_sections_joins_columns_and_excludes_excited_initial_states(tmp_path):
    """CC output splits final-state columns across differently sized blocks."""
    path = tmp_path / "xsec"
    path.write_text(
        "# CROSS SECTIONS IN BOHR**2 FOR INITIAL STATE 1\n"
        "# I E TOTAL 1 2\n"
        "1 0.1 2.0 2.0 0.0\n2 0.2 3.0 3.0 0.0\n\n"
        "# CROSS SECTIONS IN BOHR**2 FOR INITIAL STATE 1\n"
        "# I E TOTAL 7\n"
        "1 0.1 0.5\n2 0.2 0.7\n"
        "# CROSS SECTIONS IN BOHR**2 FOR INITIAL STATE 2\n"
        "# I E TOTAL 1\n"
        "1 0.1 9.0 9.0\n2 0.2 8.0 8.0\n"
    )
    np.testing.assert_allclose(
        ground_cross_sections(path),
        [[1, 0.1, 2, 2, 0, 0.5], [2, 0.2, 3, 3, 0, 0.7]],
        atol=0,
        rtol=0,
    )


@pytest.mark.parametrize(
    ("excited_error", "native_dipole"),
    [(0, "-0.400000000000D-01"), (2e-6, "-0.400000000000D-01"), (0, "0.400000000000D-01")],
)
@pytest.mark.parametrize("nonuniform_roots", [False, True])
def test_state_average_checks_excited_roots_not_ensemble_energy(
    tmp_path, excited_error, native_dipole, nonuniform_roots
):
    """Ground agreement must not hide a broken excited-state orbital import."""
    geometry = tmp_path / "output/CO/geom1"
    geometry.mkdir(parents=True)
    molden = geometry / "co.molden"
    molden.write_text("retained orbital export")
    (geometry / "outputs").mkdir()
    (geometry / "outputs/target.denprop.out").write_text(
        f"1  1  0  1  0  3  1  0 {native_dipole}\n"
    )
    states = [
        {"spin": "singlet", "irrep": "A1", "root": 1, "energy_hartree": -113},
        {"spin": "triplet", "irrep": "B1", "root": 1, "energy_hartree": -112.75},
        {"spin": "triplet", "irrep": "B2", "root": 1, "energy_hartree": -112.75},
    ]
    (tmp_path / "target.json").write_text(
        json.dumps(
            {
                "converged": True,
                "states": states,
                "ground_energy_hartree": -113,
                "ensemble_energy_hartree": -112.8,
                "ground_dipole_au": [0, 0, 0.04],
                "molden_sha256": hashlib.sha256(molden.read_bytes()).hexdigest(),
            }
        )
    )
    (geometry.parent / "target.energies").write_text(
        "# R_bohr singlet.A1.1 singlet.B1.1 singlet.B2.1 "
        "triplet.B1.1 triplet.B2.1\n"
        f"2.1323 -113 -112.7 -112.7 {-112.75 + excited_error} {-112.75 + excited_error}\n"
    )
    config = {
        "bond_length": 2.1323,
        "model": "CAS-A",
        "orbitals": "state-averaged",
        "target_roots": 1,
    }
    if nonuniform_roots:
        config.update(
            target_roots=5,
            target_singlet_roots=[1, 1, 1, 0],
            target_triplet_roots=[0, 1, 1, 0],
        )
    if excited_error or not native_dipole.startswith("-"):
        with pytest.raises(AssertionError):
            target_properties(tmp_path, config)
    else:
        result = target_properties(tmp_path, config)
        assert result["target_backend"] == "pyscf"
        assert result["target_state_max_energy_difference_hartree"] < 1e-12
        assert result["target_energy_difference_hartree"] < 1e-12
        assert abs(result["target_dipole_difference_au"]) < 1e-12
