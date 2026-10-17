"""Decision-changing failures in the prospective reconstruction control plane."""

import copy
import hashlib
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co import native_full_reconstruction as runner
from projects.ukrmol_co.native_reconstruction_checks import (
    complete_coverage,
    configurations,
    energy_checks,
)


@pytest.fixture
def contract():
    return json.loads(
        Path(__file__).with_name("native-full-reconstruction-contract.json").read_text()
    )


def family(contract, perturbation=None):
    rows = []
    for key in configurations(contract, True):
        k = np.zeros((6, 6)) if perturbation is None else perturbation(key)
        rows.append(
            dict(
                completed=True,
                step_bohr=key[0],
                endpoint_bohr=key[1],
                tail_order=key[2],
                k=k.tolist(),
            )
        )
    return rows


def test_phase_sum_veto_even_when_every_k_and_s_entry_passes(contract):
    # Six independently small eigenphases accumulate above the determinant budget.
    rows = family(contract, lambda key: np.eye(6) * np.tan(0.2e-6))
    checks, metrics = energy_checks(rows, np.zeros((6, 6)), contract, True)
    assert checks["native_k"] and checks["native_s"]
    assert not checks["native_phase_rad_mod_pi"]
    assert metrics["native_comparison"]["phase_rad_mod_pi"] == pytest.approx(1.2e-6)


@pytest.mark.parametrize("knob", ["step", "endpoint", "order"])
def test_each_final_refinement_is_an_independent_veto(contract, knob):
    rows = family(contract)
    selected = {"step": (0.025, 400, 8), "endpoint": (0.0125, 800, 8), "order": (0.0125, 400, 6)}[
        knob
    ]
    for row in rows:
        if (row["step_bohr"], row["endpoint_bohr"], row["tail_order"]) == selected:
            row["k"] = (np.eye(6) * 1e-4).tolist()
    checks, _ = energy_checks(rows, np.zeros((6, 6)), contract, True)
    assert checks["native_s"]
    assert not checks[knob + "_s"]


def test_coarse_unsymmetrized_matrix_can_veto_invariants(contract):
    rows = family(contract)
    rows[0]["k"][0][1] = 1e-3
    checks, _ = energy_checks(rows, np.zeros((6, 6)), contract, True)
    assert checks["native_s"]
    assert not checks["reciprocity"] and not checks["unitarity"]


def test_fourth_order_family_and_full_grid_scope(contract):
    def fourth_order(key):
        k = np.zeros((6, 6))
        k[0, 0] = key[0] ** 4
        return k

    checks, _ = energy_checks(family(contract, fourth_order), np.zeros((6, 6)), contract, True)
    assert all(checks.values())
    assert len(configurations(contract, True)) == 7
    assert len(configurations(contract, False)) == 2
    assert 18 * 7 + (594 - 18) * 2 == contract["independent_solves_if_complete"]


def test_duplicate_or_missing_configuration_cannot_pass(contract):
    rows = family(contract)
    rows[-1] = copy.deepcopy(rows[0])
    with pytest.raises(ValueError, match="duplicated"):
        energy_checks(rows, np.zeros((6, 6)), contract, True)


def test_duplicate_energy_count_cannot_establish_full_coverage(contract):
    rows = [
        dict(
            case=case,
            energy_index=index,
            verdict=True,
            rows=[{}] * (7 if index in [0, 48, 98] else 2),
        )
        for case in contract["cases"]
        for index in range(99)
    ]
    report = dict(probes=rows, remainder=[])
    assert complete_coverage(report, contract)
    rows[-1] = copy.deepcopy(rows[0])
    assert not complete_coverage(report, contract)


@pytest.mark.parametrize("verdict", [False, None])
def test_false_or_unset_stops_new_submissions(verdict):
    executed, received = [], []

    def worker(index):
        executed.append(index)
        return dict(verdict=verdict, index=index)

    with ThreadPoolExecutor(max_workers=1) as executor:
        result = runner.gated_pool(executor, worker, range(10), 1, received.append)
    assert result is verdict
    assert executed == [0] and len(received) == 1


@pytest.mark.parametrize("timeout", [True, False])
def test_native_failure_preserves_actual_exit(tmp_path, monkeypatch, contract, timeout):
    source = tmp_path / "input" / "compressed/nine/B1"
    source.mkdir(parents=True)
    for name in ["fort.10", "fort.21"]:
        (source / name).write_bytes(b"pinned input")
        contract["input_sha256"]["compressed/nine/B1/" + name] = hashlib.sha256(
            b"pinned input"
        ).hexdigest()
    (source / "rsolve.inp").write_bytes(b"&bprop\n&end\n")
    engine = tmp_path / "engine"
    body = "import time\ntime.sleep(30)" if timeout else "print('PROPAGATED ACROSS 8 SUBRANGES')"
    engine.write_text(f"#!{sys.executable}\n{body}\n")
    engine.chmod(0o755)
    contract["resources"]["native_stage_timeout_seconds"] = 0.05 if timeout else 10
    monkeypatch.setattr(runner, "load_model", lambda source: {"unit": "fort.921"})
    result = runner.native_case(
        (tmp_path / "input", tmp_path / "output", engine, "compressed/nine/B1", contract)
    )
    assert result["verdict"] is None
    assert len(result["stages"]) == 1
    assert result["stages"][0]["timed_out"] is timeout
    if timeout:
        assert result["stages"][0]["original_exit"] < 0
    else:
        assert result["stages"][0]["original_exit"] == 0
    assert not result["stages"][0]["completed"]


@pytest.mark.parametrize("failed_stage", ["native", "probes"])
def test_stage_failure_cannot_authorize_a_later_stage(
    tmp_path, monkeypatch, contract, failed_stage
):
    from argparse import Namespace

    prerequisite = {
        "attribution_control_verdict": True,
        "original_pilot_verdict": False,
        "public_reconstruction_passed": True,
    }
    (tmp_path / "prerequisite.json").write_text(json.dumps(prerequisite))
    contract["prerequisite_sha256"] = {
        "prerequisite.json": runner.digest(tmp_path / "prerequisite.json")
    }
    engine = tmp_path / "engine"
    engine.write_bytes(b"frozen")
    contract["engine_rsolve_sha256"] = runner.digest(engine)
    monkeypatch.setattr(runner, "verify_inputs", lambda *args: None)
    monkeypatch.setattr(runner, "ProcessPoolExecutor", ThreadPoolExecutor)
    submissions = []

    def native(task):
        return dict(case=task[3], verdict=failed_stage != "native")

    def independent(task):
        submissions.append(task[-1])
        return dict(case=task[2], energy_index=task[3], verdict=False)

    monkeypatch.setattr(runner, "native_case", native)
    monkeypatch.setattr(runner, "independent_energy", independent)
    args = Namespace(directory=tmp_path, output=tmp_path, engine=engine, prerequisites=tmp_path)
    report = dict(native=[], probes=[], remainder=[], stage_verdicts={})
    runner.run(args, contract, report)
    assert report["reconstruction_verdict"] is False
    assert report["remainder"] == []
    assert not submissions if failed_stage == "native" else all(submissions)
    assert len(submissions) <= contract["resources"]["workers"]
