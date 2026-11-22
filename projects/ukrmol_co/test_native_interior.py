"""Generated native row inventories, coherent assembly and spectral-frame controls."""

import struct

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from projects.ukrmol_co.core_sector_overlap import CoreSectorCI
from projects.ukrmol_co.native_interior import (
    assemble_native_interior,
    continuum_row_layout,
    read_retained_spectrum,
    recover_retained_state,
)


@pytest.mark.parametrize("payload", [2, 5])
def test_spectrum_explicit_full_or_partial_inventory_and_exact_eof(tmp_path, payload):
    path = tmp_path / "spectrum"
    header = bytearray(208)
    struct.pack_into("<4q", header, 136, 0, 5, 3, 1)
    struct.pack_into("<d", header, 192, -2.0)
    vectors = np.arange(3 * payload).reshape(3, payload) / 10
    records = [
        b"********                CIDATA  ",
        header,
        np.r_[np.zeros(5), [0.1, 0.2, 0.3], np.zeros(5)].astype("<f8").tobytes(),
    ]
    records += [struct.pack("<q", i + 1) + v.astype("<f8").tobytes() for i, v in enumerate(vectors)]
    with path.open("xb") as stream:
        for body in records:
            marker = struct.pack("<I", len(body))
            stream.write(marker + body + marker)
    energies, actual, metadata = read_retained_spectrum(
        path, dimension=5, roots=3, coefficients=payload
    )
    np.testing.assert_allclose(energies, [-1.9, -1.8, -1.7], rtol=0, atol=1e-14)
    np.testing.assert_allclose(actual, vectors, rtol=0, atol=0)
    assert metadata["coefficients"] == payload
    with pytest.raises(ValueError):
        read_retained_spectrum(path, dimension=5, roots=3, coefficients=5 if payload == 2 else 2)
    with path.open("ab") as stream:
        stream.write(b"x")
    with pytest.raises(ValueError):
        read_retained_spectrum(path, dimension=5, roots=3, coefficients=payload)


def test_explicit_native_flags_and_sector_state_continuum_order():
    layout = continuum_row_layout([1, 1, 1, 2, 2], [0, 1, 1, 0, 1], (2, 1), (1, 0), (1, 2))
    assert layout == ((2,), (0, 1))


@pytest.mark.parametrize("case", ["counts", "flags", "label", "roots"])
def test_inconsistent_row_metadata_rejects(case):
    kwargs = dict(
        orbital_irreps=[1, 1, 2],
        flags=[0, 1, 1],
        target_counts=(1,),
        continuum_irreps=(0,),
        continuum_counts=(1,),
    )
    if case == "counts":
        kwargs["continuum_counts"] = (2,)
    elif case == "flags":
        kwargs["flags"] = [0, 2, 1]
    elif case == "label":
        kwargs["continuum_irreps"] = (8,)
    else:
        kwargs["target_counts"] = (65,)
    with pytest.raises(ValueError):
        continuum_row_layout(**kwargs)


def test_native_P_Q_assembly_preserves_complex_coefficients_and_input():
    target = np.zeros((1, 2, 2))
    target[0, 0, 0] = 1
    state = CoreSectorCI(0, 2, (), (), (1, 1), target)
    coefficients = np.array([2 + 3j, -4j, 0.2 + 0.1j, -0.3j])
    transform = csr_matrix([[0, -1], [1, 0]])
    result = assemble_native_interior((state,), ((1, 0),), coefficients, transform, 2)
    np.testing.assert_allclose(result.products[0].amplitudes, [[-4j, 2 + 3j]], rtol=0, atol=0)
    np.testing.assert_allclose(
        result.bound.coefficients.reshape(-1), [0.3j, 0.2 + 0.1j], rtol=0, atol=0
    )
    assert result.bound.active_electrons == (2, 1)
    np.testing.assert_array_equal(coefficients, [2 + 3j, -4j, 0.2 + 0.1j, -0.3j])
    with pytest.raises(ValueError):
        assemble_native_interior((state,), ((1, 0),), coefficients[:-1], transform, 2)


def spectral_case():
    rng = np.random.default_rng(731)
    vectors = np.linalg.qr(rng.normal(size=(9, 4)))[0].T
    energies = np.array([-0.3, -0.1, -0.1 + 1e-12, 0.2])
    rotation = np.eye(4)
    angle = 0.73
    rotation[0, 0] = -1
    rotation[1:3, 1:3] = [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    rotation[3, 3] = -1
    recovered = rotation @ vectors
    c = rng.normal(size=4) + 1j * rng.normal(size=4)
    return energies, vectors, recovered, c


def test_recovery_transports_signs_degenerate_rotations_and_complex_matching_c():
    energies, old, new, c = spectral_case()
    state, report = recover_retained_state(energies, old[:, :5], c, energies, new)
    np.testing.assert_allclose(state, old.T @ c, rtol=0, atol=1e-12)
    assert report["clusters"] == [1, 2, 1]
    assert report["maximum_component_error"] < 1e-12
    assert np.max(abs(new.T @ c - state)) > 0.1


@pytest.mark.parametrize("case", ["energy", "invisible", "components", "cluster_cap"])
def test_unqualified_recovered_frames_reject_instead_of_guessing(case):
    energies, old, new, c = spectral_case()
    old_p = old[:, :5].copy()
    other = energies.copy()
    kwargs = {}
    if case == "energy":
        other[-1] += 1e-5
    elif case == "invisible":
        old_p[0] = 0
        new[0] = 0
    elif case == "components":
        old_p[0, 0] += 0.1
    else:
        kwargs["max_cluster"] = 1
    with pytest.raises(ValueError):
        recover_retained_state(energies, old_p, c, other, new, **kwargs)
