"""Single-energy deck fidelity and independent vetoes in the fixed attribution."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co import native_k_attribution as runner
from projects.ukrmol_co.native_k_attribution import evaluate, single_energy_deck


@pytest.fixture
def contract():
    return json.loads(Path(__file__).with_name("native-k-attribution-contract.json").read_text())


@pytest.mark.parametrize("sector,unit", [("B1", 921), ("B2", 922)])
def test_slice_mutates_only_declared_fields(sector, unit):
    deck = (
        f"&rslvin ! {sector}\n  nescat = 99, ! count\n  einc = 0.1,0.05, ! grid\n"
        f"  lukmt = {unit},\n/\n&bprop ! comment\n  raf = 100.,\n/\n"
    ).encode()
    sliced = single_energy_deck(deck, "0.99999999999999978", 32)
    restored = (
        sliced.replace(b"nescat = 1,", b"nescat = 99,")
        .replace(b"einc = 0.99999999999999978,0.05,", b"einc = 0.1,0.05,")
        .replace(b"  nrange = 32,\n", b"")
    )
    assert restored == deck
    assert single_energy_deck(deck, "0.99999999999999978", 64) == sliced.replace(
        b"nrange = 32", b"nrange = 64"
    )


@pytest.mark.parametrize(
    "deck",
    [
        b"&rslvin\nnescat=98,\neinc=0.1,0.05,\n/\n&bprop\n/\n",
        b"nescat=99,\nnescat=99,\neinc=0.1,0.05,\n&bprop\n/\n",
    ],
)
def test_unfrozen_or_duplicate_energy_fields_are_rejected(deck):
    with pytest.raises(ValueError, match="Exactly one"):
        single_energy_deck(deck, "0.99999999999999978", 32)


def test_incident_energy_literal_reproduces_native_double_bits():
    assert float("0.99999999999999978") * 0.0735 == 0.07349999999999998


def rows(contract):
    keys = [(0.025, 400.0, 8), (0.0125, 400.0, 8)] + contract["independent_configurations"]
    return [
        {"step_bohr": h, "endpoint_bohr": r, "tail_order": p, "k": np.zeros((6, 6)).tolist()}
        for h, r, p in keys
    ]


@pytest.mark.parametrize(
    "selected,failed",
    [
        ((0.00625, 400.0, 8), "step_k"),
        ((0.0125, 800.0, 8), "endpoint_k"),
        ((0.0125, 400.0, 6), "order_k"),
    ],
)
def test_each_control_keeps_its_separate_k_veto(contract, selected, failed):
    family = rows(contract)
    # Near a tangent singularity the K veto can fail while S and phase pass.
    for row in family:
        k = np.diag([-185.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        if (row["step_bohr"], row["endpoint_bohr"], row["tail_order"]) == selected:
            k[0, 0] += 1e-5
        row["k"] = k.tolist()
    report = evaluate(family[:2], family[2:], np.diag([-185.0, 0.0, 0.0, 0.0, 0.0, 0.0]), contract)
    assert report["verdict"] is False
    assert not report["checks"][failed]
    assert all(
        value
        for key, value in report["checks"].items()
        if key.endswith("_s") or key.endswith("_phase_rad_mod_pi")
    )


def test_missing_control_cannot_be_reported_as_a_pass(contract):
    family = rows(contract)
    with pytest.raises(ValueError, match="Missing"):
        evaluate(family[:2], family[2:-1], np.zeros((6, 6)), contract)


@pytest.mark.parametrize("failure", ["slice_budget", "energy_identity", "native64_budget"])
def test_native_gates_preserve_exits_and_stop_refinement(tmp_path, monkeypatch, contract, failure):
    case = "stretched/nine/B2"
    raw, prior, output = (tmp_path / name for name in ["raw", "prior", "output"])
    source = raw / case
    source.mkdir(parents=True)
    for name in ["fort.10", "fort.21"]:
        (source / name).write_bytes(b"frozen input")
        contract["input_sha256"][case + "/" + name] = runner.digest(source / name)
    (source / "rsolve.inp").write_bytes(b"&rslvin\n nescat=99,\n einc=0.1,0.05,\n/\n&bprop\n/\n")
    engine = tmp_path / "engine"
    engine.write_text(
        f"#!{sys.executable}\nimport re,sys\ns=sys.stdin.read()\n"
        "n=re.search(r'nrange = (\\d+)',s)[1]\n"
        "print('PROPAGATED ACROSS '+n+' SUBRANGES')\n"
    )
    engine.chmod(0o755)
    header = {
        name: 0
        for name in ["target_components", "boundary_channels", "symmetry_charge", "geometry_mass"]
    }
    energy = contract["energy_rydberg"]
    monkeypatch.setattr(
        runner,
        "load_model",
        lambda path: {"energies": np.full(99, energy), "unit": "fort.922", "header": header},
    )

    def read(path):
        if path.is_relative_to(prior):
            return np.full(99, energy), np.zeros((99, 6, 6)), header
        changed = energy + 1e-15 if failure == "energy_identity" else energy
        k = np.zeros((1, 6, 6))
        if failure == "slice_budget" or (failure == "native64_budget" and path.parent.name == "64"):
            k[0, 0, 0] = 1e-5
        return np.array([changed]), k, header

    monkeypatch.setattr(runner, "read_k", read)
    result = runner.native_case((raw, prior, output, engine, case, contract))
    assert result["verdict"] is (None if failure == "energy_identity" else False)
    assert all(stage["original_exit"] == 0 for stage in result["stages"])
    assert len(result["stages"]) == (2 if failure == "native64_budget" else 1)
    assert (output / "native" / case / "64").exists() is (failure == "native64_budget")
