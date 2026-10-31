"""Independent raw-log/header checks on synthetic known truths."""

import io
import struct
import tarfile

import pytest

from projects.ukrmol_co import verify_stretched_preparation as reader
from projects.ukrmol_co.stretched_basis_preparation import cards
from projects.ukrmol_co.test_compressed_cas12_preparation import inventory, log
from projects.ukrmol_co.test_stretched_basis_preparation import SOURCE


def completed_log():
    return (
        log()
        + """
Wavefn projection (projec()) - completed
Data on CSFs has been written to file (MEGUL) = 721
Number of CSFs (nocsf) = 990990
"""
    )


def test_independent_contraction_counts_and_corruptions():
    _, scatci = cards(inventory(), SOURCE)
    result = reader.contraction(completed_log(), scatci, inventory())
    assert (
        result["contracted_dimension"] == 87198 and result["contracted_continuum_dimension"] == 2256
    )
    for broken in [
        completed_log().replace("84942 CSFS", "84941 CSFS"),
        completed_log().replace("triplet A2", "triplet B1"),
        log(),
    ]:
        with pytest.raises(ValueError):
            reader.contraction(broken, scatci, inventory())


def test_independent_native_csf_header():
    data = bytearray(304)
    data[:120] = b"e + CO - doublet B1".ljust(120)
    fields = [0] * 23
    fields[0] = 1
    for i, v in [(5, 155), (6, 310), (7, 990990), (8, 15)]:
        fields[i] = v
    struct.pack_into("<23q", data, 120, *fields)
    struct.pack_into("<2d", data, 128, 0.5, 0.5)
    record = struct.pack("<I", 304) + data + struct.pack("<I", 304)
    reader.csf_header(record, [59, 37, 37, 22], 990990)
    with pytest.raises(ValueError):
        reader.csf_header(record, [59, 37, 37, 22], 990989)
    with pytest.raises(ValueError):
        reader.csf_header(record[:-1], [59, 37, 37, 22], 990990)


def test_independent_orbital_inventory_rejects_incomplete_rank_or_sum():
    text = "\n".join(
        f"Final number of orbitals for {name}\n" + " ".join(map(str, row))
        for name, row in inventory().items()
    )
    text += "\nTwo_electron_integrals took [s]: 1\n" + "NUMBER OF ERROR MESSAGES: 0\n" * 4
    assert reader.orbital_inventory(text) == inventory()
    for broken in [
        text.replace("59 37 37 22", "58 37 37 22"),
        text.replace("MESSAGES: 0", "MESSAGES: 1", 1),
    ]:
        with pytest.raises(ValueError):
            reader.orbital_inventory(broken)


def test_stream_authenticates_uncaptured_large_binary_and_rejects_changed_pin(tmp_path):
    data = b"native binary" * 100000
    path = tmp_path / "evidence.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo("fort.721")
        member.size = len(data)
        archive.addfile(member, io.BytesIO(data))
    collection = dict(
        archive=dict(bytes=path.stat().st_size, sha256=reader.sha(path.read_bytes())),
        files={"fort.721": dict(bytes=len(data), sha256=reader.sha(data))},
    )
    kept, digests, headers = reader.stream_files(path, collection)
    assert (
        not kept and digests["fort.721"] == reader.sha(data) and headers["fort.721"] == data[:312]
    )
    collection["files"]["fort.721"]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        reader.stream_files(path, collection)
    collection["archive"]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        reader.stream_files(path, collection)
