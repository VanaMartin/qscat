"""Known-eigenphase and malformed-native-record checks for saved scattering data."""

import numpy as np
import pytest
from scipy.io import FortranFile

from projects.ukrmol_co.saved_matrices import read_k, scattering_matrix


def test_rotated_known_eigenphases_and_flux_conservation():
    angle = 0.37
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    phases = np.array([-0.4, 1.3])
    k = (rotation * np.tan(phases)) @ rotation.T
    expected = (rotation * np.exp(2j * phases)) @ rotation.T
    result = scattering_matrix(k)
    np.testing.assert_allclose(result, expected, atol=2e-15, rtol=0)
    np.testing.assert_allclose(result.conj().T @ result, np.eye(2), atol=2e-15, rtol=0)


@pytest.mark.parametrize("truncate", [False, True])
def test_native_triangle_order_and_missing_energy_rejection(tmp_path, truncate):
    path = tmp_path / "k.bin"
    with FortranFile(path, "w", header_dtype="<u4") as stream:
        stream.write_record(np.array([11, 1, 6, 1, 5], dtype="<i8"))
        stream.write_record(np.array([b"known triangle".ljust(80)], dtype="S80"))
        stream.write_record(np.array([1, 2, 0, 0], dtype="<i8"), np.array([2.1, 1.0], dtype="<f8"))
        stream.write_record(np.array([1, 0, 0, 3, 1], dtype="<i8"))
        stream.write_record(np.array([1, 2], dtype="<i8"), np.array([0.1, 0.1], dtype="<f8"))
        for energy in [0.1] if truncate else [0.1, 0.2]:
            stream.write_record(
                np.array([3, 0, 6], dtype="<i8"), np.array([energy, 1, 2, 3, 4, 5, 6], dtype="<f8")
            )
    if truncate:
        with pytest.raises(ValueError, match="Incomplete K"):
            read_k(path)
    else:
        energies, matrices, header = read_k(path)
        np.testing.assert_allclose(energies, [0.1, 0.2], atol=1e-15, rtol=0)
        np.testing.assert_array_equal(matrices[0], [[1, 2, 4], [2, 3, 5], [4, 5, 6]])
        assert header["actual_energy_records"] == 2
