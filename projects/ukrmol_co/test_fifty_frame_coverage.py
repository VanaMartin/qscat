"""Complete native-vector reading rejects changed data, not just changed counters."""

import struct

import numpy as np
import pytest
from scipy.io import FortranFile

from projects.ukrmol_co.fifty_import_verification import IRREPS, native_binary_sectors


def native_fixture(path, *, corrupt_energy=False, corrupt_vector=False):
    references = {}
    with FortranFile(path, "w") as stream:
        for irrep in IRREPS:
            for spin, name in ((0, "singlet"), (2, "triplet")):
                roots = 13 if irrep == "A1" else 8
                dimension = 20
                energies = -100 + np.arange(roots) * 0.1
                references[(spin, irrep)] = energies.copy()
                header = bytearray(208)
                title = f"CO CI calculation - {name} {irrep}"
                header[16:136] = title.ljust(120).encode()
                struct.pack_into(
                    "<qqqqddqd",
                    header,
                    136,
                    2,
                    dimension,
                    roots,
                    IRREPS.index(irrep),
                    spin / 2,
                    spin / 2,
                    14,
                    3.0,
                )
                stream.write_record(np.frombuffer(b"********                CIDATA  ", dtype="u1"))
                stream.write_record(np.frombuffer(header, dtype="u1"))
                for charge in (6, 8):
                    stream.write_record(
                        np.frombuffer(struct.pack("<qdddd", 0, 0, 0, 0, charge), dtype="u1")
                    )
                if corrupt_energy:
                    energies[-1] += 2e-7
                stream.write_record(np.r_[np.zeros(dimension), energies - 3, np.zeros(dimension)])
                for index in range(roots):
                    vector = np.eye(dimension)[index]
                    if corrupt_vector and index == 0:
                        vector *= 0.9
                    record = struct.pack("<q", index + 1) + vector.astype("<f8").tobytes()
                    stream.write_record(np.frombuffer(record, dtype="u1"))
    return references


def test_all74_native_vectors_and_headers_read_to_eof(tmp_path):
    path = tmp_path / "fort.26"
    references = native_fixture(path)
    rows = native_binary_sectors(path, references)
    assert len(rows) == 8 and sum(r["roots"] for r in rows) == 74
    assert max(r["vector_gram_error"] for r in rows) < 1e-15


@pytest.mark.parametrize("corruption", ["energy", "vector", "trailing"])
def test_native_corruptions_reject(tmp_path, corruption):
    path = tmp_path / "fort.26"
    refs = native_fixture(
        path, corrupt_energy=corruption == "energy", corrupt_vector=corruption == "vector"
    )
    if corruption == "trailing":
        with path.open("ab") as stream:
            stream.write(b"unqualified trailing bytes")
    with pytest.raises(AssertionError):
        native_binary_sectors(path, refs)


def test_reused_ci_reference_against_hubbard_analytic_and_perturbed_controls():
    pytest.importorskip("pyscf")
    from projects.ukrmol_co.coverage_reference import benchmark

    result = benchmark()
    assert result["exact_eigenvector"]["physical_residual_norm_hartree"] <= 1e-12
    assert result["perturbed_vector"]["physical_residual_norm_hartree"] >= 1e-3
