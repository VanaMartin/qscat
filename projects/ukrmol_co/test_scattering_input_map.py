"""Exact fermionic signs and malformed saved-interface controls."""

import io
import itertools
import struct

import pytest

from projects.ukrmol_co.scattering_input_map import (
    assignment,
    continuum_phase,
    orbital_irrep,
    property_map,
    record_spans,
)


def test_continuum_phase_matches_direct_target_then_continuum_order():
    # Independent oracle: adjacent swaps into the target-sorted/continuum-last order.
    for occupied in itertools.permutations((1, 2, 5, 8)):
        for continuum in occupied:
            destination = [*sorted(a for a in occupied if a != continuum), continuum]
            work = list(occupied)
            sign = 1
            for i, want in enumerate(destination):
                j = work.index(want)
                while j > i:
                    work[j], work[j - 1] = work[j - 1], work[j]
                    sign = -sign
                    j -= 1
            for coefficient in (-1, 1):
                assert continuum_phase(list(occupied), coefficient, continuum) == sign * coefficient


@pytest.mark.parametrize(
    "occupied,coefficient,continuum", [([1, 1], 1, 1), ([1, 2], 0, 2), ([1, 2], 1, 3)]
)
def test_invalid_phase_inputs_reject(occupied, coefficient, continuum):
    with pytest.raises(ValueError):
        continuum_phase(occupied, coefficient, continuum)


def test_sequential_records_require_full_framing_and_eof():
    data = struct.pack("<I4sI", 4, b"test", 4)
    assert record_spans(data) == [(4, 4)]
    for bad in (data[:-1], data + b"x", data[:8] + struct.pack("<I", 3)):
        with pytest.raises(ValueError):
            record_spans(bad)


def test_card_parser_ignores_comments_but_not_duplicate_assignments():
    assert assignment("&input\nntgts=1,2,3, ! not 99\n ntgtf=4,4,4, /", "ntgts") == [1, 2, 3]
    for text in ("/", "ntgts=1,\nntgts=2,/"):
        with pytest.raises(ValueError):
            assignment(text, "ntgts")


def test_continuum_address_boundaries():
    assert orbital_irrep(14) == (0, 7, True)
    assert orbital_irrep(15) == (0, 8, False)
    assert orbital_irrep(117) == (1, 1, True)
    assert orbital_irrep(123) == (1, 4, False)
    assert orbital_irrep(191) == (2, 1, True)
    with pytest.raises(ValueError):
        orbital_irrep(309)


def test_property_ids_are_derived_from_energy_order_not_sector_order():
    from projects.ukrmol_co.test_fifty_property_reading import fixture

    text, reference = fixture()
    text = text.replace("6 1 0 2 40 0", "6 1 0 2 40 80")
    result = property_map(text, reference, expected_moments=80)
    assert result["idtarg"][:5] == [1, 9, 17, 25, 33]
    assert len(set(result["idtarg"])) == 40


@pytest.mark.parametrize("defect", ["duplicate", "symmetry", "energy", "missing", "count"])
def test_property_corruption_rejects(defect):
    from projects.ukrmol_co.test_fifty_property_reading import fixture

    text, reference = fixture()
    text = text.replace("6 1 0 2 40 0", "6 1 0 2 40 80")
    if defect == "duplicate":
        text += "\n1 1 0 1 0 3 1 0 -0.125D+00"
    elif defect == "symmetry":
        text = text.replace("1 1 0 1 0 3 1 0", "1 1 1 1 0 3 1 0")
    elif defect == "energy":
        reference["singlet.A2"][-1] += 2e-7
    elif defect == "missing":
        text = "\n".join(line for line in text.splitlines() if not line.startswith("1 40 "))
    else:
        text = text.replace("40 80", "40 79")
    with pytest.raises(ValueError):
        property_map(text, reference, expected_moments=80)


def congen_fixture(path, monkeypatch, defect=None):
    from projects.ukrmol_co import scattering_input_map as module

    monkeypatch.setattr(module, "DIMS", (1,) * 8)
    reference = [*range(1, 13), 117, 118, 123]
    cont = [123, 15, 265, 197]
    dimensions, phases, packed, offsets = [1] * 8, [], [], [1]
    for i in range(17):
        continuum = cont[min(i // 2, 7) // 2]
        # The first two sectors have B1 continuum, the next two A1, etc.
        if i < 16:
            continuum = {0: 15, 1: 123, 2: 197, 3: 265}[1 ^ (i // 4)]
        occupied = [*reference[:-1], continuum]
        if i % 2 == 0 and i < 16:
            phases.append(continuum_phase(occupied, 1, continuum))
        packed.extend([0] if continuum == 123 else [1, 123, continuum])
        offsets.append(len(packed) + 1)
    header = bytearray(304)
    values = [1, 0, 0, 0, 0, 154, 308, 17, 15, 17, 1, 4, 2, len(packed), 0, 0, 0, 0, 0, 0, 0, 8, 8]
    struct.pack_into("<23q", header, 120, *values)
    struct.pack_into("<2d", header, 128, 0.5, 0.5)
    if defect == "phase":
        phases[3] *= -1
    phase_record = [*phases, *dimensions, *([2] * 8), 1, 1, 0, 0, 3, 3, 2, 2, *([0] * 8)]
    metadata = [
        58,
        37,
        37,
        22,
        *reference,
        *([1] * 17),
        0,
        7,
        3,
        3,
        0,
        58,
        37,
        37,
        22,
        0,
        0,
        0,
        0,
        7,
        3,
        3,
        0,
        0,
        0,
        0,
        0,
    ]
    indices = [*range(1, 19), *offsets]
    if defect == "offset":
        indices[2] = indices[1]

    def integer_record(values):
        return struct.pack(f"<{len(values)}q", *values)

    records = [
        header,
        integer_record(phase_record),
        integer_record(metadata),
        integer_record(indices),
        integer_record(packed),
        struct.pack("<17d", *([1] * 17)),
    ]
    data = b"".join(struct.pack("<I", len(r)) + r + struct.pack("<I", len(r)) for r in records)
    if defect == "trailing":
        data += struct.pack("<II", 0, 0)
    path.write_bytes(data)


def test_complete_congen_first_term_reader(tmp_path, monkeypatch):
    from projects.ukrmol_co.scattering_input_map import audit_congen

    path = tmp_path / "congen"
    congen_fixture(path, monkeypatch)
    witness = io.BytesIO()
    result = audit_congen(path, 1, [(1,)] * 8, witness)
    assert result["phase_defining_terms"] == 8 and result["l2_csfs"] == 1
    assert len(witness.getvalue()) == 8 * 32


@pytest.mark.parametrize("defect", ["phase", "offset", "trailing"])
def test_congen_phase_offset_and_eof_corruptions_reject(tmp_path, monkeypatch, defect):
    from projects.ukrmol_co.scattering_input_map import audit_congen

    path = tmp_path / "congen"
    congen_fixture(path, monkeypatch, defect)
    with pytest.raises(ValueError):
        audit_congen(path, 1, [(1,)] * 8, io.BytesIO())
