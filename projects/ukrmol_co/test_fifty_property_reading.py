"""The independent property reader checks symmetry IDs as well as field coverage."""

import pytest

from projects.ukrmol_co.fifty_property_reading import IRREPS, SECTORS, read_selected


def fixture():
    reference = {
        s: [-112 + i / 1000 + root / 10 for root in range(5)] for i, s in enumerate(SECTORS)
    }
    ordered = sorted(
        (energy, sector) for sector, energies in reference.items() for energy in energies
    )
    text = ["6 1 0 2 40 0 0 0 0 0"]
    for state_id, (energy, sector) in enumerate(ordered, 1):
        spin, irrep = sector.split(".")
        symmetry = IRREPS.index(irrep)
        text.append(f"5 {state_id} 0 0 0 1 0 0 {energy:.12f} State No. {state_id} {spin} {irrep}")
        for rank in (1, 2):
            text.append(f"1 {state_id} {symmetry} {state_id} {symmetry} 3 {rank} 0 -0.125D+00")
    return "\n".join(text), reference


def test_all_forty_states_are_explicitly_mapped():
    text, reference = fixture()
    result = read_selected(text, reference, 0.125)
    assert len(result["target_map"]) == 40
    assert len(set(result["target_map"])) == 40
    assert result["multipole_record_count"] == 80


@pytest.mark.parametrize("defect", ["symmetry", "nan", "duplicate", "missing", "energy", "dipole"])
def test_corrupt_native_records_reject(defect):
    text, reference = fixture()
    dipole = 0.125
    if defect == "symmetry":
        text = text.replace("1 1 0 1 0 3 1 0", "1 1 1 1 0 3 1 0")
    elif defect == "nan":
        text = text.replace("-0.125D+00", "NaN", 1)
    elif defect == "duplicate":
        text += "\n1 1 0 1 0 3 1 0 -0.125D+00"
    elif defect == "missing":
        text = "\n".join(line for line in text.splitlines() if not line.startswith("1 40 "))
    elif defect == "energy":
        reference["singlet.A2"][-1] += 2e-7
    else:
        dipole += 2e-5
    with pytest.raises(ValueError):
        read_selected(text, reference, dipole)
