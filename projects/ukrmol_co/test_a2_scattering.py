"""Noncontiguous channel selection and enlarged interfaces must survive native export."""

import json
import re
import sys

import numpy as np
import pytest

from projects.ukrmol_co.a2_scattering import channel_deck, channel_map, verify_model_inputs
from projects.ukrmol_co.sparse_pilot import audit_sector


def inventory():
    # Actual energy ordering of the independently checked stretched48-state import.
    maps = [
        [1, 10, 15, 28, 33],
        [4, 5, 35, 36, 38],
        [11, 19, 23, 29, 31],
        [2, 13, 16, 21, 26],
        [12, 18, 24, 30, 32],
        [3, 14, 17, 22, 25],
        [8, 9, 27, 41, 43, 44, 46, 47, 48],
        [6, 7, 20, 34, 37, 39, 40, 42, 45],
    ]
    names = [
        "singlet.A1",
        "triplet.A1",
        "singlet.B1",
        "triplet.B1",
        "singlet.B2",
        "triplet.B2",
        "singlet.A2",
        "triplet.A2",
    ]
    return {
        "state_labels_by_id": {
            str(i): f"{s}.{root}"
            for s, ids in zip(names, maps, strict=True)
            for root, i in enumerate(ids, 1)
        }
    }


def values(text, key):
    return list(map(int, re.search(rf"\b{key}\s*=\s*([\d,]+)", text)[1].rstrip(",").split(",")))


def test_baseline_uses_first_five_per_sector_not_first40_energy_ids():
    mapping = channel_map(inventory(), 5)
    assert len(set(mapping)) == 40
    assert 41 in mapping and 43 in mapping and 37 in mapping
    assert 39 not in mapping
    assert 48 not in mapping and 45 not in mapping
    assert set(mapping) != set(range(1, 41))
    expanded = channel_map(inventory(), 9)
    assert set(expanded) == set(range(1, 49))
    assert expanded[-1] == 45 and expanded[-10] == 48


def test_contract_deck_changes_only_selection_and_keeps_root_count():
    original = (
        "&input\n megul=721,\n numtgt=5,5,5,5,5,5,5,5,\n"
        " ntgtf=1,\n ntgts=1,\n/\n&cinorn\n memp=6,\n vecstore=1,\n!  igh = 1,\n/\n"
    )
    five, nine = [channel_deck(original, inventory(), n, 2048) for n in (5, 9)]
    assert values(five, "numtgt") == [5] * 8
    assert values(nine, "numtgt") == [5] * 6 + [9] * 2
    assert values(five, "ntgtf") == [i for i in range(1, 9) for _ in range(5)]
    assert values(nine, "ntgtf")[-9:] == [8] * 9
    assert values(nine, "ntgts")[-18:] == list(range(1, 10)) * 2
    assert values(nine, "idtarg") == channel_map(inventory(), 9)
    assert values(five, "idtarg") == channel_map(inventory(), 5)
    assert values(five, "nstat") == values(nine, "nstat") == [2048]
    for text in (five, nine):
        assert "memp=6" in text and "vecstore = 3" in text
        assert "write_amp = .true." in text and "auto_nvo = .false." in text

    def normalized(text):
        for field in ("numtgt", "ntgtf", "ntgts", "ntarg", "idtarg"):
            text = re.sub(rf"(\b{field}\s*=\s*)[\d,]+", r"\1CHANNELS", text)
        return text

    assert normalized(five) == normalized(nine)


@pytest.mark.parametrize("count", [6, 8, 12, 0])
def test_undeclared_channel_ladders_are_rejected(count):
    with pytest.raises(ValueError, match="five and nine"):
        channel_map(inventory(), count)


def test_missing_root9_is_not_silently_replaced_by_lowest48():
    data = inventory()
    data["state_labels_by_id"].pop("48")
    with pytest.raises(KeyError):
        channel_map(data, 9)


def test_expanded_boundary_dimension_is_checked_without_changing_default(tmp_path, monkeypatch):
    import projects.ukrmol_co.sparse_pilot as pilot

    telemetry = {
        "requested": 2,
        "converged": 2,
        "reason": 1,
        "dimension": 27818,
        "matrix_type": "mpisbaij",
        "solver_type": "krylovschur",
        "eigenpairs": [{"absolute_residual": 1e-10}] * 2,
    }
    (tmp_path / "eps-telemetry.jsonl").write_text(json.dumps(telemetry) + "\n")
    (tmp_path / "eigenp.inp").write_text("luphso=1021,")
    (tmp_path / "fort.1021").write_text(
        "\n".join(f"{e:.6f} 0.0" for e in 0.1 + 0.05 * np.arange(99))
    )
    (tmp_path / "reson.out").write_text("")
    energies = np.array([-110.0, -109.0])
    monkeypatch.setattr(pilot, "ci_spectrum", lambda path: (energies, np.ones((2, 1682))))
    monkeypatch.setattr(pilot, "boundary_spectrum", lambda path: (energies, np.ones((2, 298))))
    assert audit_sector(tmp_path, 2, 27818, 1682, 298)["export_checks_passed"]
    with pytest.raises(ValueError, match="boundary export"):
        audit_sector(tmp_path, 2, 27818, 1682)


@pytest.mark.parametrize("change", [None, "memp", "target", "congen"])
def test_actual_model_inputs_are_checked_against_frozen_frame(tmp_path, change):
    geometry = tmp_path / "reference/output/CO/geom1"
    imported = tmp_path / "import"
    work = tmp_path / "work"
    (geometry / "inputs").mkdir(parents=True)
    (geometry / "doublet.B1").mkdir()
    imported.mkdir()
    directory = work / "nine/B1"
    directory.mkdir(parents=True)
    original = (
        "&input\n megul=721,\n numtgt=5,5,5,5,5,5,5,5,\n"
        " ntgtf=1,\n ntgts=1,\n/\n&cinorn\n memp=6,\n vecstore=1,\n!  igh = 1,\n/\n"
    )
    (geometry / "inputs/scattering.scatci.doublet.B1.inp").write_text(original)
    deck = channel_deck(original, inventory(), 9, 2048)
    (directory / "scatci.inp").write_text(
        deck.replace("memp=6", "memp=16") if change == "memp" else deck
    )
    for name, parent in [
        ("fort.16", geometry / "doublet.B1"),
        ("fort.721", geometry / "doublet.B1"),
        ("fort.24", imported),
        ("fort.26", imported),
    ]:
        (parent / name).write_bytes(name.encode())
        (directory / name).write_bytes(name.encode())
    if change == "target":
        (directory / "fort.26").write_bytes(b"different vector frame")
    if change == "congen":
        (directory / "fort.721").write_bytes(b"different CSFs")
    contract = {
        "target_reference_directory": "reference",
        "target_preflight_directory": "import",
        "sectors": ["B1"],
        "scattering_roots": 2048,
    }
    if change is None:
        verify_model_inputs(tmp_path, work, contract, inventory(), "nine")
    else:
        with pytest.raises(ValueError, match="Actual scattering"):
            verify_model_inputs(tmp_path, work, contract, inventory(), "nine")


@pytest.mark.parametrize("verdict", [False, None])
def test_rejected_or_unset_baseline_never_creates_expanded_native_run(
    tmp_path, monkeypatch, verdict
):
    import projects.ukrmol_co.a2_scattering as scattering

    contract = tmp_path / "contract.json"
    contract.write_text("{}")
    (tmp_path / "baseline-qualification.json").write_text(
        json.dumps({"dense_baseline_verdict": verdict})
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "a2",
            "model",
            "--model",
            "nine",
            "--root",
            str(tmp_path),
            "--work",
            str(tmp_path),
            "--contract",
            str(contract),
        ],
    )
    with pytest.raises(ValueError, match="passing dense baseline"):
        scattering.main()
    assert not (tmp_path / "nine").exists()


@pytest.mark.parametrize("verdict", [False, None])
def test_baseline_reader_preserves_rejection_and_exits_nonzero(tmp_path, monkeypatch, verdict):
    import projects.ukrmol_co.a2_scattering as scattering

    contract = tmp_path / "contract.json"
    contract.write_text("{}")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "a2",
            "baseline-gate",
            "--root",
            str(tmp_path),
            "--work",
            str(tmp_path),
            "--contract",
            str(contract),
        ],
    )
    monkeypatch.setattr(scattering, "compare", lambda *args: {"dense_baseline_verdict": verdict})
    with pytest.raises(SystemExit) as error:
        scattering.main()
    assert error.value.code == 1
    assert (
        json.loads((tmp_path / "baseline-qualification.json").read_text())["dense_baseline_verdict"]
        is verdict
    )
