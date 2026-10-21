"""Analytic electron-slot/sign controls for the finite native CSF correspondence."""

from dataclasses import replace

import numpy as np
import pytest

from projects.ukrmol_co.native_csf_map import (
    TargetCSFs,
    address_map,
    canonical_packed,
    read_target_csfs,
    signed_correspondence,
    sorting_phase,
)


def toy():
    # Exciting alpha electron 13->15 passes beta electron 14: sorting sign -1.
    return TargetCSFs(
        np.array([7, 3, 3]),
        np.arange(1, 15),
        np.array([2]),
        np.array([1]),
        np.array([1, 3]),
        np.array([1, 5]),
        np.array([0, 1, 13, 15]),
        np.array([1, -1]) / np.sqrt(2),
        0,
        0,
    )


def test_fermionic_slot_parity_and_native_alpha_beta_embedding():
    assert sorting_phase([1, 3, 2]) == -1
    assert sorting_phase([3, 1, 2]) == 1
    assert address_map([58, 37, 37, 22])[117:123].tolist() == list(range(15, 21))
    assert address_map([58, 37, 37, 22])[191:197].tolist() == list(range(21, 27))
    assert address_map([58, 37, 37, 22])[15] == 0
    ref, packed = canonical_packed(toy())
    np.testing.assert_array_equal(ref, np.arange(1, 15))
    np.testing.assert_array_equal(packed, toy().packed)


def test_negative_csf_orientation_and_involutory_vector_map():
    source = toy()
    dest = replace(
        source,
        inventory=np.array([58, 37, 37, 22]),
        phases=-source.phases,
        coefficients=-source.coefficients,
        packed=np.array([0, 1, 13, 117]),
    )
    signs = signed_correspondence(source, dest)
    np.testing.assert_array_equal(signs, [-1])
    vectors = np.array([[0.3, -0.4]])
    assert (signs[:, None] * signs[:, None] * vectors).tobytes() == vectors.tobytes()


@pytest.mark.parametrize("defect", ["phase", "nonleading", "continuum", "missing", "duplicate"])
def test_rejects_false_correspondence(defect):
    source = toy()
    if defect == "phase":
        dest = replace(source, phases=-source.phases)
    elif defect == "nonleading":
        dest = replace(source, coefficients=np.array([1, 1]) / np.sqrt(2))
    elif defect == "continuum":
        dest = replace(
            source, inventory=np.array([58, 37, 37, 22]), packed=np.array([0, 2, 13, 14, 16, 15])
        )
    elif defect == "missing":
        dest = replace(source, counts=np.array([3]))
    else:
        dest = replace(source, packed=np.array([0, 2, 13, 14, 13, 13]))
    with pytest.raises(ValueError):
        signed_correspondence(source, dest)


def write_toy(path):
    from scipy.io import FortranFile

    csf = toy()
    header = np.zeros(23, dtype="<i8")
    header[[5, 6, 7, 8, 9, 10, 11, 12, 13, 21, 22]] = [13, 26, 1, 14, 2, 1, 3, 2, 4, 1, -1]
    metadata = np.r_[
        csf.inventory,
        csf.reference,
        csf.counts,
        0,
        csf.inventory,
        csf.inventory,
        [0] * 3,
        csf.inventory,
        [0] * 3,
    ].astype("<i8")
    with FortranFile(path, "w") as stream:
        stream.write_record(
            np.frombuffer(b"CO - singlet A1".ljust(120, b" ") + header.tobytes(), "u1")
        )
        for array in (
            csf.phases,
            metadata,
            np.r_[csf.coefficient_offsets, csf.packed_offsets],
            csf.packed,
            csf.coefficients,
        ):
            stream.write_record(array)


def test_six_record_reader_requires_complete_data_and_eof(tmp_path):
    path = tmp_path / "fort.710"
    write_toy(path)
    np.testing.assert_array_equal(read_target_csfs(path).packed, toy().packed)
    with path.open("ab") as stream:
        stream.write(b"extra")
    with pytest.raises(ValueError, match="Trailing"):
        read_target_csfs(path)
    write_toy(path)
    path.write_bytes(path.read_bytes()[:-1])
    from scipy.io import FortranFormattingError

    with pytest.raises(FortranFormattingError):
        read_target_csfs(path)


def test_separate_canonical_decoder_recovers_analytic_superposition(tmp_path):
    from projects.ukrmol_co.native_csf_reconstruction import expansion

    path = tmp_path / "fort.710"
    write_toy(path)
    csfs = list(expansion(path))
    assert len(csfs) == 1 and len(csfs[0]) == 2
    # Stored second coefficient is negative; the electron-slot sorting is odd.
    np.testing.assert_allclose(list(csfs[0].values()), [1 / np.sqrt(2)] * 2, atol=0, rtol=0)
    excited = tuple(
        sorted([(0, k, s) for k in range(1, 7) for s in (0, 1)] + [(0, 7, 1), (1, 1, 0)])
    )
    assert excited in csfs[0]
