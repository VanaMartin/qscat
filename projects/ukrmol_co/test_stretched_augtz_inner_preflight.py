"""Catch a newly open excited target before allocating the expensive inner export."""

from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co.stretched_augtz_inner_preflight import DECK, threshold_screen


def energies():
    sectors = [f"{s}.{i}" for i in ["A1", "B1", "B2", "A2"] for s in ["singlet", "triplet"]]
    e = {s: [-112.0 + i for i in range(8)] for s in sectors}
    e["singlet.A1"][0] = -113.0
    return e


def test_known_closed_grid_pass_and_aug_like_crossing_reject_without_retuning():
    e = energies()
    deck = Path(DECK).read_text()
    assert threshold_screen(e, deck)["unchanged_six_open_grid_prerequisite_verdict"] is True
    for sector in ["triplet.B1", "triplet.B2"]:
        e[sector][0] = -113.0 + 0.18159458918989912
    r = threshold_screen(e, deck)
    assert r["unchanged_six_open_grid_prerequisite_verdict"] is False
    assert r["affected_grid_indices_zero_based"] == [97, 98]
    np.testing.assert_allclose(r["affected_requested_energies_ev"], [4.95, 5.0], atol=1e-12, rtol=0)
    assert len(r["crossed_retained_target_thresholds"]) == 2
    assert r["expensive_inner_execution_released"] is False
    assert r["new_native_calls"] == r["new_ci_calls"] == r["new_maps"] == 0


def test_narrowed_grid_is_not_silently_presented_as_the_same_contract():
    deck = Path(DECK).read_text().replace("nescat = 99", "nescat = 97")
    with pytest.raises(ValueError, match="unchanged legacy"):
        threshold_screen(energies(), deck)


@pytest.mark.parametrize(
    "change",
    [
        ("nerang = 1", "nerang = 2"),
        ("ieunit = 2", "ieunit = 1"),
        ("einc = 0.1, 0.05", "einc = 0.1, 0.04"),
    ],
)
def test_other_energy_prescription_changes_are_rejected(change):
    deck = Path(DECK).read_text().replace(*change)
    with pytest.raises(ValueError, match="unchanged legacy"):
        threshold_screen(energies(), deck)


def test_invalid_inventory_or_nonfinite_native_energy_stops():
    e = energies()
    deck = Path(DECK).read_text()
    del e["triplet.B1"]
    with pytest.raises(ValueError, match="inventory"):
        threshold_screen(e, deck)
    e = energies()
    e["triplet.B1"][0] = float("nan")
    with pytest.raises(ValueError, match="Finite"):
        threshold_screen(e, deck)


def test_threshold_on_grid_endpoint_fails_the_strict_closed_channel_condition():
    e = {s: [float(i + 1) for i in range(8)] for s in energies()}
    e["singlet.A1"][0] = 0.0
    e["triplet.B1"][0] = 0.18375
    r = threshold_screen(e, Path(DECK).read_text())
    assert r["unchanged_six_open_grid_prerequisite_verdict"] is False
    assert r["affected_grid_indices_zero_based"] == [98]
