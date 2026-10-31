"""Saved-only independent reader controls; no molecular/native calculation."""

import pytest

from projects.ukrmol_co.test_stretched_basis_feasibility import fixture
from projects.ukrmol_co.test_stretched_basis_properties import table
from projects.ukrmol_co.verify_stretched_properties import inventory, native_energies


def test_independent_native_reader_known_truth():
    energies, norm = native_energies(fixture())
    assert len(energies) == 8 and norm == 0
    assert all(e == list(range(-100, -92)) for e in energies.values())


@pytest.mark.parametrize("defect", ["count", "phase", "norm", "root", "trailing"])
def test_independent_native_reader_rejects_corruption(defect):
    with pytest.raises(ValueError):
        native_energies(fixture(defect))


def test_independent_inventory_known_truth_and_rejections():
    text, expected = table()
    checked = inventory(text, expected, 0)
    assert checked["target_count"] == 64 and checked["roots_per_sector"] == [8] * 8
    for broken in [
        text.replace("6 0 0 0 64", "6 0 0 0 40"),
        text.replace("-100 ", "-99 ", 1),
        "\n".join(line for line in text.splitlines() if not line.startswith("1 64 ")),
        text + "\n1 1 0 1 0 3 1 0 0.0",
    ]:
        with pytest.raises(ValueError):
            inventory(broken, expected, 0)
    with pytest.raises(ValueError):
        inventory(text, expected, 0.1)
