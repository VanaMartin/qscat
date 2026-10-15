"""Explicit-channel maps must cover the requested sectors rather than the lowest N states."""

import json

import numpy as np
import pytest

from projects.ukrmol_co.a2_channel_preflight import (
    SECTORS,
    expected_states,
    property_inventory,
    replace_assignment,
    target_solver_deck,
)
from projects.ukrmol_co.sparse_comparison_owner import sha


def property_table(a2_roots=9):
    expected = {}
    states = []
    for index, (spin, irrep) in enumerate(SECTORS):
        label = f"{spin}.{irrep}"
        count = a2_roots if irrep == "A2" else 5
        # Interleave sectors in energy order; A2 root9 need not be in the lowest40.
        expected[label] = [-112 + root * 0.1 + index * 0.001 for root in range(count)]
        states.extend((e, spin, irrep) for e in expected[label])
    states.sort()
    rows = [f"6 1 699 2 {len(states)} 656 0 0 0 0"]
    for state_id, (energy, spin, irrep) in enumerate(states, 1):
        rows.append(f"5 {state_id} 0 0 0 1 0 0 {energy:.12f} State No. {state_id} {spin} {irrep}")
        for rank in (1, 2):
            rows.append(f"1 {state_id} 0 {state_id} 0 3 {rank} 0 0.125D+00")
    return "\n".join(rows), expected


def test_nonuniform_explicit_a2_map_contains_both_ninth_roots():
    text, expected = property_table()
    result = property_inventory(text, expected)
    assert result["target_count"] == 48
    labels = result["state_labels_by_id"]
    assert [labels[i] for i in result["target_map"]] == [
        f"{spin}.{irrep}.{root}"
        for spin, irrep in SECTORS
        for root in range(1, len(expected[f"{spin}.{irrep}"]) + 1)
    ]
    assert labels[result["target_map"][-1]] == "triplet.A2.9"
    assert result["ground_dipole_z_au"] == -0.125
    assert result["multipole_record_count"] == 96


def test_lowest40_does_not_qualify_explicit48_target_selection():
    text, _ = property_table(a2_roots=5)
    _, expected = property_table(a2_roots=9)
    with pytest.raises(ValueError, match="header/state IDs"):
        property_inventory(text, expected)


def test_wrong_sector_counts_cannot_pass_total_target_count():
    text, expected = property_table()
    expected["singlet.A2"] = expected["singlet.A2"][:-1]
    expected["singlet.A1"] = [*expected["singlet.A1"], -111.9]
    with pytest.raises(ValueError, match="Incomplete retained target sector"):
        property_inventory(text, expected)


@pytest.mark.parametrize(
    "replacement,match",
    [
        ("1 999 0 1 0 3 1 0 0.125D+00", "absent state"),
        ("1 1 0 1 0 3 3 0 0.125D+00", "unexpected rank"),
        ("1 1 0 1 0 3 1 0 NaN", "nonfinite"),
    ],
)
def test_invalid_multipoles_fail_before_scattering(replacement, match):
    text, expected = property_table()
    text = text.replace("1 1 0 1 0 3 1 0 0.125D+00", replacement)
    with pytest.raises(ValueError, match=match):
        property_inventory(text, expected)


def test_missing_added_state_properties_fail():
    text, expected = property_table()
    text = "\n".join(line for line in text.splitlines() if not line.startswith("1 48 "))
    with pytest.raises(ValueError, match="coverage"):
        property_inventory(text, expected)


def test_changed_native_energy_fails_without_relaxed_budget():
    text, expected = property_table()
    expected["triplet.A2"][-1] += 1e-5
    with pytest.raises(ValueError, match="energy gate failed"):
        property_inventory(text, expected)


def test_duplicate_native_assignment_is_not_silently_patched():
    with pytest.raises(ValueError, match="one native assignment"):
        replace_assignment("nstat=5,\nnstat=8,\n", "nstat", "9")
    with pytest.raises(ValueError, match="one native assignment"):
        replace_assignment("ntgt=8,\n", "ntgtl", "5,5,5,5,5,5,9,9")


def test_target_deck_change_keeps_solver_tolerance_and_comment():
    text = "  igh=-1, forse=0, crite=1e-12, maxiter=500, ! solver\n  nstat=5, ! roots\n"
    changed = replace_assignment(replace_assignment(text, "nstat", "9"), "maxiter", "1000")
    assert "crite=1e-12" in changed and "maxiter=1000, ! solver" in changed
    assert "nstat=9, ! roots" in changed
    np.testing.assert_equal(changed.count("nstat="), 1)


def test_equilibrium_dense_target_template_selects_solver_without_losing_full_vectors():
    original = (
        "&cinorn\n nstat=5, ! roots\n nciset=1, ! index\n vecstore=5,\n memp=6.0,\n!igh = 1,\n/\n"
    )
    changed = target_solver_deck(original, 9, 8)
    assert "nstat=9," in changed and "nciset=8," in changed
    assert "igh = -1, crite = 1e-12, maxiter = 1000," in changed
    assert "vecstore=5," in changed and "memp=6.0," in changed
    assert target_solver_deck(changed, 5, 1).count("igh = -1") == 1


@pytest.mark.parametrize("slot", ["igh=1,", "!igh=2,", "!igh=1,\n!igh=1,", "crite=1e-6,\n!igh=1,"])
def test_unexpected_dense_target_solver_is_rejected(slot):
    with pytest.raises(ValueError, match=r"solver settings|diagonalizer slot"):
        target_solver_deck("nstat=5,\nnciset=1,\n" + slot + "\n", 9, 8)


@pytest.mark.parametrize("defect", [None, "frame", "residual", "spin", "flags", "settings"])
def test_native_import_requires_same_frame_independent_actions(tmp_path, defect):
    run = tmp_path / "reference"
    geometry = run / "output/CO/geom1"
    geometry.mkdir(parents=True)
    independent = tmp_path / "independent"
    independent.mkdir()
    (geometry / "co.casscf.chk").write_bytes(b"immutable checkpoint")
    (geometry / "co.molden").write_bytes(b"immutable orbital export")
    _, energies = property_table()
    states = [
        {"spin": spin, "irrep": irrep, "energy_hartree": energy}
        for spin, irrep in SECTORS
        for energy in energies[f"{spin}.{irrep}"][:5]
    ]
    (run / "target.json").write_text(
        json.dumps(
            {"converged": True, "molden_sha256": sha(geometry / "co.molden"), "states": states}
        )
    )
    (run / "config.json").write_text(json.dumps({"bond_length": 2.5}))
    for spin, irrep in SECTORS:
        data = {
            "checkpoint_sha256": sha(geometry / "co.casscf.chk"),
            "gate_passed": True,
            "ci_converged": [True] * 9,
            "energy_hartree": energies[f"{spin}.{irrep}"],
            "physical_residual_norms_hartree": [1e-10] * 9,
            "penalized_residual_norms_hartree": [1e-10] * 9,
            "spin_square": [0 if spin == "singlet" else 2] * 9,
        }
        if (spin, irrep) == ("triplet", "A2"):
            if defect == "frame":
                data["checkpoint_sha256"] = "other frame"
            if defect == "residual":
                data["physical_residual_norms_hartree"][-1] = 1e-5
            if defect == "spin":
                data["spin_square"][-1] = 0
            if defect == "flags":
                data["ci_converged"][-1] = False
        (independent / f"{spin}.{irrep}.json").write_text(json.dumps(data))
    contract = {
        "target_reference_directory": "reference",
        "independent_reference_directory": "independent",
        "fixed_target_settings": {"bond_length": 2.4 if defect == "settings" else 2.5},
    }
    if defect is None:
        assert expected_states(tmp_path, contract) == energies
    else:
        with pytest.raises(ValueError):
            expected_states(tmp_path, contract)
