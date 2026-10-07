"""Guard native binary parsing and the immutable selected-root contraction deck."""

import struct

import numpy as np
import pytest
from scipy.io import FortranFile

from projects.ukrmol_co.sparse_scattering import (
    boundary_spectrum,
    ci_spectrum,
    dense_root_indices,
    interface_digest,
    mpi_boundary_deck,
    selected_deck,
    truncate_boundary,
)


def write_ci(path, last_index=2):
    # Independent fixture following the pinned SolutionHandler binary layout.
    header = (
        struct.pack("<qq", 1, 4)
        + b"synthetic".ljust(120)
        + struct.pack("<qqqqddqdq", 1, 4, 2, 1, 0.5, 0.5, 15, -100.0, 0)
    )
    with FortranFile(path, "w") as stream:
        stream.write_record(np.frombuffer(b"********                CIDATA  ", dtype="u1"))
        stream.write_record(np.frombuffer(header, dtype="u1"))
        stream.write_record(np.zeros(40, dtype="u1"))
        stream.write_record(np.array([1.0, 2.0, 3.0, 4.0, -1.0, -0.5, 9.0, 8.0, 7.0, 6.0]))
        for index, coeff in [(1, [0.25, -0.5]), (last_index, [0.75, 0.125])]:
            record = struct.pack("<qdd", index, *coeff)
            stream.write_record(np.frombuffer(record, dtype="u1"))


def test_selected_ci_has_dimension_sized_phase_and_core_shift(tmp_path):
    path = tmp_path / "fort.25"
    write_ci(path)
    energies, vectors = ci_spectrum(path)
    np.testing.assert_allclose(energies, [-101.0, -100.5], atol=1e-14, rtol=0)
    np.testing.assert_allclose(vectors, [[0.25, -0.5], [0.75, 0.125]], atol=1e-14, rtol=0)


def test_ci_rejects_incomplete_root_sequence(tmp_path):
    path = tmp_path / "fort.25"
    write_ci(path, last_index=3)
    with pytest.raises(ValueError, match="Incomplete CIDATA"):
        ci_spectrum(path)


def test_dense_root_matching_handles_endpoints_and_selected_order():
    np.testing.assert_array_equal(
        dense_root_indices(
            np.array([-99.1, -112.9, -109.9]), np.array([-113.0, -111.0, -110.0, -99.0])
        ),
        [3, 0, 2],
    )
    with pytest.raises(ValueError, match="distinct"):
        dense_root_indices(np.array([-113.0, -112.9]), np.array([-113.0, -111.0]))
    with pytest.raises(ValueError, match="ordered"):
        dense_root_indices(np.array([-113.0]), np.array([-111.0, -113.0]))


def test_scattering_selection_preserves_target_and_export_contract():
    text = "&input\n nftg = 26, numtgt = 5,5,\n/\n&cinorn\n vecstore = 1,\n!  igh = 1,\n/\n"
    selected = selected_deck(text, 128, 1e-12, 1000)
    assert selected.split("&cinorn")[0] == text.split("&cinorn")[0]
    assert "vecstore = 1," in selected
    assert "igh = -1, nstat = 128," in selected
    with pytest.raises(ValueError, match="already selects"):
        selected_deck(selected, 512, 1e-12, 1000)


@pytest.mark.parametrize("slots", [0, 2])
def test_scattering_selection_rejects_ambiguous_native_dispatch(slots):
    with pytest.raises(ValueError, match="unique native"):
        selected_deck("&cinorn\n" + "!  igh = 1,\n" * slots + "/\n", 128, 1e-12, 1000)


def test_native_boundary_export_reuses_target_order_and_enables_memory_vectors():
    text = "&input\n nftg = 26, numtgt = 5,5,\n/\n&cinorn\n vecstore = 1,\n/\n"
    selected = mpi_boundary_deck(text, "ntarg = 3, idtarg = 2,1,3, rmatr = 18,")
    assert selected.split("&cinorn")[0] == text.split("&cinorn")[0]
    assert "vecstore = 3," in selected
    assert "ntarg = 3, idtarg = 2,1,3," in selected
    assert "write_amp = .true., write_dip = .false., write_rmt = .false." in selected
    with pytest.raises(ValueError, match="Incomplete retained-state"):
        mpi_boundary_deck(text, "ntarg = 3, idtarg = 2,1, rmatr = 18,")


def write_boundary(path, poles=2, buttle=0, amplitude_count=None, multipoles=2):
    with FortranFile(path, "w") as stream:
        stream.write_record(np.array([11, 1, 7, 1, 1], dtype="<i8"))
        stream.write_record(np.frombuffer(b"native boundary".ljust(80), dtype="u1"))
        stream.write_record(np.array([40, 0, 0, 3], dtype="<i8"))
        stream.write_record(np.array([1, 2, 2, 0, 0, 0], dtype="<i8"))
        data = struct.pack("<qqqqd", multipoles, poles, 0, buttle, 18.0)
        stream.write_record(np.frombuffer(data, dtype="u1"))
        if multipoles:
            stream.write_record(np.zeros(12, dtype="<f8"))
        stream.write_record(-112.0 + np.arange(poles) * 0.1)
        stream.write_record(
            np.arange(3 * poles if amplitude_count is None else amplitude_count, dtype="<f8")
        )


def test_boundary_parser_matches_channel_by_pole_native_order(tmp_path):
    path = tmp_path / "fort.21"
    write_boundary(path)
    energies, amplitudes = boundary_spectrum(path)
    np.testing.assert_allclose(energies, [-112.0, -111.9], atol=1e-14, rtol=0)
    np.testing.assert_array_equal(amplitudes, [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    other = tmp_path / "dense.21"
    write_boundary(other, poles=4)
    assert interface_digest(path, boundary=True) == interface_digest(other, boundary=True)


@pytest.mark.parametrize("multipoles", [0, 2])
def test_dense_oracle_truncation_preserves_full_bytes_and_static_interface(tmp_path, multipoles):
    source = tmp_path / "dense.21"
    write_boundary(source, poles=4, multipoles=multipoles)
    full = tmp_path / "full.21"
    truncate_boundary(source, full, 4)
    assert full.read_bytes() == source.read_bytes()
    selected = tmp_path / "subset.21"
    truncate_boundary(source, selected, 2)
    energy, amplitudes = boundary_spectrum(selected)
    np.testing.assert_allclose(energy, [-112.0, -111.9], atol=1e-14, rtol=0)
    np.testing.assert_array_equal(amplitudes, [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    assert interface_digest(selected, boundary=True) == interface_digest(source, boundary=True)
    with pytest.raises(FileExistsError):
        truncate_boundary(source, selected, 2)
    with pytest.raises(ValueError, match="valid root count"):
        truncate_boundary(source, tmp_path / "invalid.21", 5)


@pytest.mark.parametrize("buttle,amplitudes,message", [(1, 6, "Partitioned"), (0, 5, "Incomplete")])
def test_boundary_parser_rejects_partitioned_and_incomplete_records(
    tmp_path, buttle, amplitudes, message
):
    path = tmp_path / "fort.21"
    write_boundary(path, buttle=buttle, amplitude_count=amplitudes)
    with pytest.raises(ValueError, match=message):
        boundary_spectrum(path)
