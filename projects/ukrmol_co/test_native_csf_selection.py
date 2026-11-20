"""Hand-built selected-CSF format/parity controls; no molecular fixtures."""

import struct

import numpy as np
import pytest

from projects.ukrmol_co.native_csf_selection import (
    corrected_target_coefficients,
    selected_csf_determinants,
)


def fixture_file(path, *, reference=(1, 2), coefficients=None, groups=False):
    coefficients = (
        np.array([1 / np.sqrt(2), 1 / np.sqrt(2), 1.0]) if coefficients is None else coefficients
    )
    # Two CSFs: a two-determinant open-shell singlet and a closed-shell state.
    substitutions = [[1, 2, 4], [1, 1, 3], [0]]
    packed = np.concatenate(substitutions).astype("<i8")
    fields = [
        0,
        0.0,
        0.0,
        0.0,
        0.0,
        3,
        6,
        2,
        2,
        3,
        1,
        1,
        2,
        len(packed),
        *([0] * 6),
        1e-12,
        2,
        1 if groups else -1,
    ]
    header = b"generated selected CSF".ljust(120) + struct.pack("<qdddd9q6qdqq", *fields)
    phases = np.array([1, -1, 2, 1, 0, 0] if groups else [1, -1], "<i8")
    meta = np.array([3, *reference, 2, 1, 0, 3, 3, 0, 3, 0], "<i8")
    offsets = np.array([1, 3, 4, 1, 7, 8], "<i8")
    rows = [
        header,
        phases.tobytes(),
        meta.tobytes(),
        offsets.tobytes(),
        packed.tobytes(),
        np.asarray(coefficients, "<f8").tobytes(),
    ]
    with path.open("wb") as stream:
        for body in rows:
            marker = struct.pack("<I", len(body))
            stream.write(marker + body + marker)


def test_selected_signed_determinants_and_source_metadata(tmp_path):
    path = tmp_path / "csf"
    fixture_file(path)
    result = selected_csf_determinants(path, [0])
    assert result["norb"] == 3 and result["inventory"] == (3,)
    assert set(result["selected"]) == {0}
    assert result["selected"][0] == {(0, 4): 1 / np.sqrt(2), (1, 3): 1 / np.sqrt(2)}
    np.testing.assert_array_equal(result["phase_labels"], [1, -1])


def test_reference_slot_parity_is_not_a_native_sort_assumption(tmp_path):
    path = tmp_path / "csf"
    fixture_file(path, reference=(2, 1))
    value = selected_csf_determinants(path, [1])["selected"][1]
    assert value == {(0, 3): -1.0}


def test_unselected_nonfinite_coefficients_are_not_expanded(tmp_path):
    path = tmp_path / "csf"
    fixture_file(path, coefficients=np.array([np.nan, np.nan, 1.0]))
    assert selected_csf_determinants(path, [1])["selected"][1] == {(0, 3): 1.0}
    with pytest.raises(ValueError, match="nonfinite"):
        selected_csf_determinants(path, [0])


def test_explicit_target_group_metadata_and_selected_term_cap(tmp_path):
    path = tmp_path / "csf"
    fixture_file(path, groups=True)
    assert selected_csf_determinants(path, [1])["target_groups"] == ((2,), (1,), (0,), (0,))
    with pytest.raises(ValueError, match="term cap"):
        selected_csf_determinants(path, [0], max_terms=1)


@pytest.mark.parametrize("indices", [[], [0, 0], [-1], [2], list(range(33))])
def test_selected_index_guards(tmp_path, indices):
    path = tmp_path / "csf"
    fixture_file(path)
    with pytest.raises(ValueError):
        selected_csf_determinants(path, indices)


def test_corrupt_record_trailer_is_rejected(tmp_path):
    path = tmp_path / "csf"
    fixture_file(path)
    body = bytearray(path.read_bytes())
    body[-1] ^= 1
    path.write_bytes(body)
    with pytest.raises(ValueError, match="trailer"):
        selected_csf_determinants(path, [0])


def test_native_row_phase_product_is_explicit_complex_and_covariant():
    coefficients = np.array([[1 + 2j, 3j], [4.0, 5 - 1j]])
    original = coefficients.copy()
    corrected = corrected_target_coefficients(coefficients, [1, -1], [-1, -1])
    np.testing.assert_allclose(corrected, [[-1 - 2j, -3j], [4.0, 5 - 1j]], rtol=0, atol=0)
    np.testing.assert_array_equal(coefficients, original)
    np.testing.assert_array_equal(
        corrected_target_coefficients(corrected, [1, -1], [-1, -1]), original
    )
    with pytest.raises(ValueError):
        corrected_target_coefficients(coefficients, [0, 1], [1, 1])
