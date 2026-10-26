"""Known-truth binary and raw-log controls without a molecular calculation."""

import struct

import numpy as np
import pytest

from projects.ukrmol_co.stretched_basis_feasibility import (
    SECTORS,
    native_tokens,
    read_cidata,
    records,
)


def fixture(defect=None):
    chunks = []

    def put(data):
        chunks.append(struct.pack("<I", len(data)) + data + struct.pack("<I", len(data)))

    for index, sector in enumerate(SECTORS, 1):
        put(b"********                CIDATA  ")
        header = bytearray(208)
        struct.pack_into("<qq", header, 0, index, 11)
        header[16:136] = ("CO CI calculation - " + sector.replace(".", " ")).encode().ljust(120)
        spin = float(index % 2 == 0)
        struct.pack_into(
            "<4q2dqdq",
            header,
            136,
            2,
            8,
            7 if defect == "count" else 8,
            (index - 1) // 2,
            spin,
            spin,
            14,
            -100.0,
            0,
        )
        put(header)
        for charge in (6, 8):
            put(b"atom    " + struct.pack("<4d", 0, 0, 0, charge))
        phases = np.ones(8, dtype="<i8")
        if defect == "phase":
            phases[0] = 0
        put(phases.tobytes() + np.arange(8, dtype="<f8").tobytes() + bytes(64))
        for root in range(8):
            vector = np.eye(8, dtype="<f8")[root]
            if defect == "norm":
                vector *= 2
            put(struct.pack("<q", root + (2 if defect == "root" else 1)) + vector.tobytes())
    return b"".join(chunks) + (b"bad" if defect == "trailing" else b"")


def test_all64_binary_roots_are_read_without_rounding():
    answer = read_cidata(fixture())
    assert set(answer) == set(SECTORS)
    for row in answer.values():
        np.testing.assert_allclose(row["energies"], np.arange(8) - 100, atol=0, rtol=0)
        assert row["orthogonality_error"] <= 1e-14


@pytest.mark.parametrize("defect", ["count", "phase", "norm", "root", "trailing"])
def test_bad_native_records_are_not_qualified(defect):
    with pytest.raises(ValueError):
        read_cidata(fixture(defect))


def test_bad_fortran_tail_is_rejected():
    with pytest.raises(ValueError):
        list(records(struct.pack("<I", 2) + b"hi" + struct.pack("<I", 3)))


def test_inherited_energy_echo_is_not_a_new_sector():
    echo = "Old CIDATA energies -88.0\nCI data will be stored as set number 3\nEIGEN-ENERGIES\n"
    number, energies = native_tokens(echo + " ".join(str(-100.0 + i) for i in range(8)) + "\nTimer")
    assert number == 3 and len(energies) == 8
    with pytest.raises(ValueError):
        native_tokens(echo + "-100.0\nTimer")
