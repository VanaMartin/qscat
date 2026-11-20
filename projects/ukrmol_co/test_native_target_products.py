"""Generated raw-CIDATA/CSF phase, sparse contraction and signed row controls."""

import struct

import numpy as np
import pytest

from projects.ukrmol_co.core_sector_overlap import CoreSectorCI
from projects.ukrmol_co.native_target_products import (
    expand_primary_target,
    native_target_products,
    read_primary_cidata,
)


def framed(stream, body):
    marker = struct.pack("<I", len(body))
    stream.write(marker + body + marker)


def write_cidata(path, *, phases=(1, -1), charge=8.0, bad_order=False, trailing=False):
    with path.open("xb") as stream:
        for irrep in range(4):
            for spin in range(2):
                framed(stream, b"********                CIDATA  ")
                header = bytearray(208)
                struct.pack_into("<2q", header, 0, 2 * irrep + spin + 1, 5)
                struct.pack_into("<4q", header, 136, 2, 2, 2, irrep)
                struct.pack_into("<2d", header, 168, spin, spin)
                struct.pack_into("<q", header, 184, 14)
                framed(stream, header)
                for q in [6.0, charge]:
                    framed(stream, struct.pack("<5d", 0.0, 0.0, 0.0, 0.0, q))
                framed(
                    stream,
                    np.array(phases, "<i8").tobytes()
                    + np.array([0.1, 0.2], "<f8").tobytes()
                    + bytes(16),
                )
                for root in range(2):
                    framed(
                        stream,
                        struct.pack("<q", 2 if bad_order else root + 1)
                        + np.eye(2, dtype="<f8")[:, root].tobytes(),
                    )
        if trailing:
            stream.write(b"x")


def write_target(path, *, spin=0, phase=(1, -1)):
    reference = list(range(1, 15)) if spin == 0 else [*range(1, 13), 13, 15]
    substitutions = [((), ()), ((13, 14), (15, 16))] if spin == 0 else [((), ()), ((15,), (17,))]
    packed = []
    po = [1]
    for holes, particles in substitutions:
        packed.extend([len(holes), *holes, *particles])
        po.append(len(packed) + 1)
    fields = [0, float(spin), float(spin), 0.0, 0.0, 13, 26, 2, 14, 2, 1, 3, 2, len(packed)]
    fields += [0] * 6 + [1e-10, 2, -1]
    bodies = [
        b"Generated target".ljust(120) + struct.pack("<qdddd9q6qdqq", *fields),
        np.array(phase, "<i8").tobytes(),
        np.array([7, 3, 3, *reference, 1, 1], "<i8").tobytes(),
        np.array([1, 2, 3, *po], "<i8").tobytes(),
        np.array(packed, "<i8").tobytes(),
        np.array([1.0, 1.0], "<f8").tobytes(),
    ]
    with path.open("xb") as stream:
        for body in bodies:
            framed(stream, body)


def test_cidata_preserves_raw_configuration_state_axes_and_phases(tmp_path):
    path = tmp_path / "cidata"
    write_cidata(path)
    sectors = read_primary_cidata(path)
    assert len(sectors) == 8
    assert [(s["irrep"], s["spin"]) for s in sectors] == [
        (i, s) for i in range(4) for s in range(2)
    ]
    for sector in sectors:
        np.testing.assert_array_equal(sector["phases"], [1, -1])
        np.testing.assert_allclose(sector["coefficients"], np.eye(2), rtol=0, atol=0)


@pytest.mark.parametrize(
    "kwargs", [dict(phases=(0, 1)), dict(charge=7.0), dict(bad_order=True), dict(trailing=True)]
)
def test_cidata_rejects_mismatched_source_or_framing(tmp_path, kwargs):
    path = tmp_path / "cidata"
    write_cidata(path, **kwargs)
    with pytest.raises(ValueError):
        read_primary_cidata(path)


@pytest.mark.parametrize("spin", [0, 1])
def test_full_raw_target_transform_preserves_complex_root_and_no_extra_phase(tmp_path, spin):
    path = tmp_path / "csf"
    write_target(path, spin=spin)
    c = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
    state = expand_primary_target(path, c, [1, -1], spin)
    assert state.active_electrons == (5 + spin, 5 - spin)
    flat = state.coefficients.reshape(2, -1)
    np.testing.assert_allclose(flat.conj() @ flat.T, np.eye(2), rtol=0, atol=1e-13)
    identity = expand_primary_target(path, np.eye(2), [1, -1], spin).coefficients.reshape(2, -1)
    np.testing.assert_allclose(flat, c.T @ identity, rtol=0, atol=1e-13)
    # Applying delta again would change the second CSF's sign; raw C is untouched.
    np.testing.assert_array_equal(c, np.array([[1, 1j], [1j, 1]]) / np.sqrt(2))
    with pytest.raises(ValueError, match="phases"):
        expand_primary_target(path, c, [-1, 1], spin)


def tiny_targets():
    singlet = np.zeros((1, 2, 2))
    singlet[0, 0, 0] = 1
    triplet = np.ones((1, 1, 1))
    return (
        CoreSectorCI(0, 2, (), (), (1, 1), singlet),
        CoreSectorCI(0, 2, (), (), (2, 0), triplet),
    )


def test_declared_native_row_order_and_frozen_signed_factors():
    singlet, triplet = tiny_targets()
    rows = np.array([1 + 2j, 3 - 4j, 5 + 6j])
    alpha, beta = native_target_products((singlet, triplet), ((2, 0), (1,)), rows, 3)
    np.testing.assert_allclose(
        alpha.amplitudes, [[rows[1], 0, rows[0]], [0, -rows[2] / np.sqrt(3), 0]], rtol=0, atol=1e-14
    )
    np.testing.assert_allclose(
        beta.amplitudes, [[0, rows[2] * np.sqrt(2 / 3), 0]], rtol=0, atol=1e-14
    )
    # S-minus on |alpha0 alpha1> is |alpha0 beta1>-|alpha1 beta0>.
    np.testing.assert_allclose(
        alpha.targets.coefficients[1],
        [[0, 1 / np.sqrt(2)], [-1 / np.sqrt(2), 0]],
        rtol=0,
        atol=1e-14,
    )
    assert alpha.continuum_spin == 0 and beta.continuum_spin == 1
    np.testing.assert_array_equal(rows, [1 + 2j, 3 - 4j, 5 + 6j])
    # Arbitrary root rephasing must be compensated in channel rows, not silently fixed.
    phase = np.exp(0.71j)
    rephased = CoreSectorCI(0, 2, (), (), (2, 0), phase * triplet.coefficients)
    aa, bb = native_target_products((singlet, rephased), ((2, 0), (1,)), rows / [1, 1, phase], 3)
    for before, after in [(alpha, aa), (beta, bb)]:
        np.testing.assert_allclose(
            np.einsum("iab,ic->abc", before.targets.coefficients, before.amplitudes),
            np.einsum("iab,ic->abc", after.targets.coefficients, after.amplitudes),
            rtol=0,
            atol=1e-13,
        )


@pytest.mark.parametrize("case", ["duplicate", "outside", "rows", "cap", "nonfinite"])
def test_declared_row_layout_resource_and_finiteness_guards(case):
    sectors = tiny_targets()
    indices = ((0,), (1,))
    rows = np.ones(2)
    count = 2
    if case == "duplicate":
        indices = ((0, 0), (1,))
    elif case == "outside":
        indices = ((2,), (1,))
    elif case == "rows":
        rows = np.ones(3)
    elif case == "cap":
        count = 161
    else:
        rows[0] = np.nan
    with pytest.raises(ValueError):
        native_target_products(sectors, indices, rows, count)
