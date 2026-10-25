"""Matched-property inventory controls; no native call or quantum calculation."""

import unittest
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co.a2_channel_preflight import property_inventory
from projects.ukrmol_co.compressed_cas12_properties import matched_deck

DECK = """&input
 ntgt=8,
 nftsor=710,730,711,731,712,732,713,733,
 ntgtf=1,2,3,4,5,6,7,8,
 ntgtl=6,4,5,7,5,7,3,3, ! Original lowest40 is not firstfive per sector.
 nftg=26, nftmt=24, nftint=16,
/
"""


def test_only_the_original_property_selection_changes():
    actual = matched_deck(DECK)
    assert actual.replace("5,5,5,5,5,5,5,5, ", "6,4,5,7,5,7,3,3, ") == DECK
    for text in [
        DECK.replace("ntgtf=1,2", "ntgtf=2,1"),
        DECK.replace("ntgtl=6,4", "ntgtl=5,5"),
        DECK.replace("ntgt=8", "ntgt=7"),
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            matched_deck(text)


def test_first_five_inventory_rejects_original_uneven_counts_and_energy_corruption():
    sectors = [s + "." + i for i in ["A1", "B1", "B2", "A2"] for s in ["singlet", "triplet"]]
    expected = {s: [-100 + 0.1 * (i * 5 + j) for j in range(5)] for i, s in enumerate(sectors)}
    rows = ["6 0 0 0 40"]
    for i in range(40):
        spin, irrep = sectors[i // 5].split(".")
        rows.append(f"5 {i + 1} 0 0 0 0 0 0 {-100 + 0.1 * i:.17g} {spin} {irrep}")
        for rank in [1, 2]:
            rows.append(f"1 {i + 1} 0 {i + 1} 0 3 {rank} 0 0.0")
    text = "\n".join(rows)
    result = property_inventory(text, expected)
    assert result["target_count"] == 40 and result["target_map"] == list(range(1, 41))
    np.testing.assert_allclose(result["maximum_energy_error_hartree"], 0, atol=1e-14, rtol=0)
    with unittest.TestCase().assertRaises(ValueError):
        property_inventory(text.replace("singlet A2", "singlet A1", 1), expected)
    with unittest.TestCase().assertRaises(ValueError):
        property_inventory(text.replace("-100 ", "-99.999 ", 1), expected)


def test_deck_translation_never_launches_a_solver_or_native_program():
    with patch("subprocess.Popen") as process:
        matched_deck(DECK)
    process.assert_not_called()


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
