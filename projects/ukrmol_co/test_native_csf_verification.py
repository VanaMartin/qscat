"""Mixed integer-phase/real-energy CIDATA controls without a target solver."""

import struct
from dataclasses import replace

import numpy as np
import pytest
from scipy.io import FortranFile

from projects.ukrmol_co.native_csf_verification import verify_cidata
from projects.ukrmol_co.test_native_csf_map import toy


def cidata_fixture(path, defect=None):
    pairs = {}
    with FortranFile(path, "w") as stream:
        for symmetry, irrep in enumerate(("A1", "B1", "B2", "A2")):
            for spin, label in ((0, "singlet"), (1, "triplet")):
                sector = f"{label}.{irrep}"
                source = replace(
                    toy(), symmetry=symmetry, spin=spin, phases=np.ones(16, dtype="<i8")
                )
                signs = np.ones(16, dtype="<i8")
                signs[0] = -1
                pairs[sector] = (source, source, signs)
                roots = 13 if symmetry == 0 else 8
                header = bytearray(208)
                struct.pack_into("<qq", header, 0, len(pairs), roots + 3)
                header[16:136] = f"CO CI calculation - {label} {irrep}".encode().ljust(120, b" ")
                centers = 3 if defect == "centres" else 2
                struct.pack_into(
                    "<qqqqddqdq",
                    header,
                    136,
                    centers,
                    16,
                    roots,
                    symmetry,
                    spin,
                    spin,
                    14,
                    -100.0,
                    0,
                )
                stream.write_record(np.frombuffer(b"********                CIDATA  ", "u1"))
                stream.write_record(np.frombuffer(header, "u1"))
                for charge in (6, 8):
                    stream.write_record(
                        np.frombuffer(b"atom    " + struct.pack("<4d", 0, 0, 0, charge), "u1")
                    )
                phases = source.phases.copy()
                if defect == "phase":
                    phases[0] = -1
                energies = np.arange(roots, dtype="<f8")
                if defect == "energies":
                    energies[0] = np.nan
                stream.write_record(
                    np.frombuffer(
                        phases.tobytes() + energies.tobytes() + np.zeros(16, dtype="<f8").tobytes(),
                        "u1",
                    )
                )
                for root in range(roots):
                    index = root + 2 if defect == "root" else root + 1
                    vector = np.eye(roots, 16, dtype="<f8")[root]
                    stream.write_record(
                        np.frombuffer(struct.pack("<q", index) + vector.tobytes(), "u1")
                    )
        if defect == "trailing":
            stream.write_record(np.zeros(1))
    return pairs


def test_all_74_root_records_and_integer_phases_roundtrip(tmp_path):
    path = tmp_path / "fort.26"
    pairs = cidata_fixture(path)
    result = verify_cidata(path, pairs)
    assert sum(s["roots"] for s in result) == 74
    assert all(s["exact_coefficient_roundtrip"] for s in result)


@pytest.mark.parametrize("defect", ["phase", "energies", "root", "centres", "trailing"])
def test_saved_cidata_defects_do_not_pass(tmp_path, defect):
    path = tmp_path / "fort.26"
    pairs = cidata_fixture(path, defect)
    with pytest.raises(ValueError):
        verify_cidata(path, pairs)
