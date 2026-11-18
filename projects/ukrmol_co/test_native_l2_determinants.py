"""Synthetic packed-CSF framing, parity, sparse transform and rejection controls."""

import struct

import numpy as np
import pytest

from projects.ukrmol_co.native_l2_determinants import l2_determinant_transform


def write_csfs(path, *, continuum_in_tail=False, nonunit=False, mixed=False):
    # Five electrons, two closed core electrons, three bound spatial orbitals.
    # The third A1 orbital is a continuum placeholder, B1 starts at spatial3.
    reference = [1, 2, 3, 4, 5]
    terms = [
        [([5], [6], 1.0)],  # Discarded continuum CSF.
        [([5], [7], 1.0)],
        [([4, 5], [7, 8], 1.0)],
    ]
    if continuum_in_tail:
        terms[-1] = [([4, 5], [6, 7], 1.0)]
    if mixed:
        r = 1 / np.sqrt(2)
        terms[1:] = [
            [([5], [7], r), ([4, 5], [7, 8], r)],
            [([5], [7], r), ([4, 5], [7, 8], -r)],
        ]
    packed, coefficients, co, po = [], [], [1], [1]
    for csf in terms:
        for holes, particles, coefficient in csf:
            packed.extend([len(holes), *holes, *particles])
            coefficients.append(coefficient * (0.5 if nonunit else 1))
        co.append(len(coefficients) + 1)
        po.append(len(packed) + 1)
    fields = [0, 0.5, 0.5, 0.0, 0.0, 5, 10, 3, 5, len(coefficients), 1, 2, 2, len(packed)]
    fields += [0] * 6 + [1e-10, 1, 1]
    header = b"Synthetic doublet".ljust(120) + struct.pack("<qdddd9q6qdqq", *fields)
    metadata = [3, 2, *reference, *[len(csf) for csf in terms]]
    records = [
        header,
        np.array([1], "<i8").tobytes(),
        np.array(metadata, "<i8").tobytes(),
        np.array(co + po, "<i8").tobytes(),
        np.array(packed, "<i8").tobytes(),
        np.array(coefficients, "<f8").tobytes(),
    ]
    with path.open("xb") as stream:
        for record in records:
            marker = struct.pack("<I", len(record))
            stream.write(marker + record + marker)


def test_only_bound_tail_expands_with_full_occupied_parity(tmp_path):
    path = tmp_path / "csfs"
    write_csfs(path)
    transform = l2_determinant_transform(path, (2, 1), 1, (2, 1), 1)
    # The native occupied lists need three and two swaps, respectively.
    np.testing.assert_allclose(transform.toarray(), np.diag([-1, 1]), rtol=0, atol=0)
    np.testing.assert_allclose((transform.T @ transform).toarray(), np.eye(2), rtol=0, atol=0)
    q = np.array([1, 1j]) / np.sqrt(2)
    np.testing.assert_allclose(np.vdot(transform @ q, transform @ q), 1, rtol=0, atol=1e-14)


def test_multiterm_orthonormal_csfs_preserve_complex_combinations(tmp_path):
    path = tmp_path / "csfs"
    write_csfs(path, mixed=True)
    transform = l2_determinant_transform(path, (2, 1), 1, (2, 1), 1)
    r = 1 / np.sqrt(2)
    np.testing.assert_allclose(transform.toarray(), [[-r, -r], [r, -r]], rtol=0, atol=1e-14)
    np.testing.assert_allclose((transform.T @ transform).toarray(), np.eye(2), rtol=0, atol=1e-14)


@pytest.mark.parametrize("kwargs", [dict(continuum_in_tail=True), dict(nonunit=True)])
def test_continuum_leakage_and_nonunit_coefficients_reject(tmp_path, kwargs):
    path = tmp_path / "csfs"
    write_csfs(path, **kwargs)
    with pytest.raises(ValueError):
        l2_determinant_transform(path, (2, 1), 1, (2, 1), 1)


def test_caps_spin_tail_and_record_integrity_reject(tmp_path):
    path = tmp_path / "csfs"
    write_csfs(path)
    with pytest.raises(ValueError, match="term cap"):
        l2_determinant_transform(path, (2, 1), 1, (2, 1), 1, max_terms=1)
    for counts, core, electrons, last in [
        ((2, 1), 1, (1, 2), 1),
        ((2, 1), 0, (2, 1), 1),
        ((2, 1), 1, (2, 1), 3),
        ((2, 1), 1, (2, 1), -1),
    ]:
        with pytest.raises(ValueError):
            l2_determinant_transform(path, counts, core, electrons, last)
    malformed = tmp_path / "malformed"
    malformed.write_bytes(path.read_bytes()[:-1])
    with pytest.raises(ValueError):
        l2_determinant_transform(malformed, (2, 1), 1, (2, 1), 1)
