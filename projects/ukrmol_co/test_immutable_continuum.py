"""Continuum preparation changes only native addresses and rejects incomplete files."""

import struct

import pytest

from projects.ukrmol_co.immutable_continuum import congen_deck, csf_records, native_inventory


def inventory_text():
    rows = {
        "target": [7, 3, 3, 0],
        "PCOs": [0, 0, 0, 0],
        "TGT+PCOs": [7, 3, 3, 0],
        "continuum": [51, 34, 34, 22],
        "total": [58, 37, 37, 22],
    }
    return (
        "\n".join(
            f"Final number of orbitals for {kind}  \n" + " ".join(map(str, counts + [0] * 4))
            for kind, counts in rows.items()
        )
        + "\nTwo_electron_integrals took [s]: 5\n"
        + "NUMBER OF ERROR MESSAGES: 0\n" * 4
    )


def test_native_inventory_reconstructs_full_continuum():
    assert native_inventory(inventory_text())["total"] == [58, 37, 37, 22, 0, 0, 0, 0]


@pytest.mark.parametrize("defect", ["duplicate", "missing", "error", "asymmetric", "total"])
def test_bad_native_inventory_rejects(defect):
    text = inventory_text()
    if defect == "duplicate":
        text += "\nFinal number of orbitals for target\n7 3 3 0 0 0 0 0\n"
    elif defect == "missing":
        text = text.replace("Final number of orbitals for PCOs", "Missing PCOs")
    elif defect == "error":
        text = text.replace("ERROR MESSAGES: 0", "ERROR MESSAGES: 1", 1)
    elif defect == "asymmetric":
        text = text.replace("51 34 34 22", "51 33 34 22")
    else:
        text = text.replace("58 37 37 22", "57 37 37 22")
    with pytest.raises(ValueError):
        native_inventory(text)


def test_only_congen_full_orbital_address_counts_change():
    text = "&state\n  nob  = 58,37,37,22, ! total\n  nob0 = 7,3,3,0,\n/\n&wfngrp\n untouched\n/\n"
    result = congen_deck(text, [59, 38, 38, 22, 0, 0, 0, 0])
    assert result.replace("59,38,38,22,", "58,37,37,22,") == text
    with pytest.raises(ValueError):
        congen_deck(text.replace("58,37,37,22", "59,37,37,22"), [58, 37, 37, 22, 0, 0, 0, 0])


def csf_fixture(path):
    header = [0] * 23
    for index, value in {5: 154, 6: 308, 7: 15548, 8: 14, 10: 1}.items():
        header[index] = value
    rows = [b"CO - singlet A1".ljust(120, b" ") + struct.pack("<23q", *header)] + [
        b"saved CSF record"
    ] * 5
    with path.open("wb") as stream:
        for data in rows:
            marker = struct.pack("<I", len(data))
            stream.write(marker + data + marker)


def test_csf_header_records_and_complete_eof(tmp_path):
    path = tmp_path / "fort.710"
    csf_fixture(path)
    assert len(csf_records(path, [58, 37, 37, 22, 0, 0, 0, 0], "singlet.A1", 15548)) == 6
    with pytest.raises(ValueError):
        csf_records(path, [58, 37, 37, 22, 0, 0, 0, 0], "triplet.A1", 15548)
    with path.open("ab") as stream:
        stream.write(b"junk")
    with pytest.raises(ValueError):
        csf_records(path, [58, 37, 37, 22, 0, 0, 0, 0], "singlet.A1", 15548)
