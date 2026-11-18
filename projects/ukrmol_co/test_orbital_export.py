"""Generated double/int64 records; no molecular inputs or native executable."""

import struct

import numpy as np
import pytest
from scipy.io import FortranFormattingError

from projects.ukrmol_co.orbital_export import read_orbital_export


def payloads():
    def integers(values):
        return np.array(values, "<i8").tobytes(order="F")

    def floats(values):
        return np.array(values, "<f8").tobytes(order="F")

    return [
        b"CO_ORBITAL_EXPORT_V1".ljust(32),
        integers([0, 3, 2, 0, 1, 1]),
        floats([18]),
        integers([2]),
        integers([[1, 2], [1, 1], [1, 2]]),
        integers([0, 1]),
        floats([[1, 4], [2, 5], [3, 6]]),
        floats([[1, 0], [0, 1]]),
        integers([[1, 1, 1], [1, 2, 3]]),
        integers([[1], [1], [0], [1], [3], [1]]),
        b"",
        integers([1, 1, 0]),
        floats([0, 0, 0, 1]),
        floats([0.25]),
        floats([-0.5]),
        floats([2]),
    ]


def write_records(path, rows):
    with path.open("xb") as stream:
        for body in rows:
            marker = struct.pack("<I", len(body))
            stream.write(marker + body + marker)


def test_full_column_major_export_preserves_coefficients_and_normalization(tmp_path):
    path = tmp_path / "orbitals.bin"
    write_records(path, payloads())
    original = path.read_bytes()
    result = read_orbital_export(path)
    assert result["kind"] == 0 and result["radius_bohr"] == 18
    np.testing.assert_allclose(result["coefficients"], [[1, 4], [2, 5], [3, 6]], rtol=0, atol=0)
    np.testing.assert_array_equal(result["continuum_flags"], [0, 1])
    np.testing.assert_allclose(result["shells"][0]["contractions"], [-0.5], rtol=0, atol=0)
    np.testing.assert_allclose(result["shells"][0]["primitive_norms"], [2], rtol=0, atol=0)
    assert path.read_bytes() == original


@pytest.mark.parametrize(
    "change",
    ["nan", "counts", "ao-address", "shell-angular", "negative-exponent", "extra", "truncated"],
)
def test_malformed_or_incomplete_export_rejects(tmp_path, change):
    rows = payloads()
    if change == "nan":
        rows[6] = np.full(6, np.nan, "<f8").tobytes()
    elif change == "counts":
        rows[3] = np.array([10**12], "<i8").tobytes()
    elif change == "ao-address":
        rows[8] = np.array([[1, 1, 1], [1, 1, 3]], "<i8").tobytes(order="F")
    elif change == "shell-angular":
        rows[11] = np.array([2, 1, 0], "<i8").tobytes()
    elif change == "negative-exponent":
        rows[13] = np.array([-1], "<f8").tobytes()
    elif change == "extra":
        rows.append(b"junk")
    path = tmp_path / "bad.bin"
    write_records(path, rows)
    if change == "truncated":
        path.write_bytes(path.read_bytes()[:-5])
    with pytest.raises((ValueError, OSError, FortranFormattingError)):
        read_orbital_export(path)
