"""Test expanded interfaces and decisions that cannot erase the original phase rejection."""

import json
import struct

import numpy as np
import pytest
from scipy.io import FortranFile

from projects.ukrmol_co.a2_root_control import assess_decision, refined_deck
from projects.ukrmol_co.sparse_comparison import read_run
from projects.ukrmol_co.test_sparse_scattering import write_ci


def test_root_replay_changes_exactly_the_selected_count():
    text = "&cinorn\n memp=6, igh=-1, nstat = 2048, crite=1e-12,\n/\nntarg=48,"
    candidate = refined_deck(text)
    assert candidate == text.replace("nstat = 2048", "nstat = 4096")


@pytest.mark.parametrize("text", ["nstat=4096,", "nstat=2048,nstat=2048,", "igh=-1,"])
def test_unmatched_or_ambiguous_reference_roots_stop(text):
    with pytest.raises(ValueError, match="exactly one"):
        refined_deck(text)


@pytest.mark.parametrize(
    "pair,cross,phase_change,expected",
    [
        (True, False, 0.001, True),
        (True, True, 0.001, False),
        (False, False, 0.001, False),
        (True, False, 0.01, False),
        (None, False, 0.001, None),
        (True, None, 0.001, None),
    ],
)
def test_attribution_requires_tight_refinement_and_persistent_rejection(
    pair, cross, phase_change, expected
):
    result = assess_decision(
        {"root_pair_stability_verdict": pair},
        {
            s: {"observable_budget_verdict": cross, "phase_budget_passed": cross}
            for s in ("B1", "B2")
        },
        {s: phase_change for s in ("B1", "B2")},
        0.005,
    )
    assert result["resolution_decision_verdict"] is expected
    assert result["original_a2_channel_sensitivity_verdict"] is False
    assert result["spectral_truncation_verdict"] is None


def test_both_pi_sectors_must_preserve_phase_rejection():
    cross = {
        "B1": {"observable_budget_verdict": False, "phase_budget_passed": False},
        "B2": {"observable_budget_verdict": True, "phase_budget_passed": True},
    }
    result = assess_decision(
        {"root_pair_stability_verdict": True}, cross, {"B1": 0.001, "B2": 0.001}, 0.005
    )
    assert result["resolution_decision_verdict"] is False


def test_expanded_raw_run_reader_uses_declared298_channel_binary_shape(tmp_path):
    saved = {"all_stages_complete": True, "target_map": list(range(1, 49)), "sectors": {}}
    for sector in ("B1", "B2"):
        p = tmp_path / sector
        p.mkdir()
        write_ci(p / "fort.25")
        with FortranFile(p / "fort.21", "w") as stream:
            stream.write_record(np.array([11, 1, 7, 1, 1], dtype="<i8"))
            stream.write_record(np.frombuffer(b"expanded fixture".ljust(80), dtype="u1"))
            stream.write_record(np.array([48, 0, 0, 298], dtype="<i8"))
            stream.write_record(np.array([1, 2, 2, 0, 0, 0], dtype="<i8"))
            stream.write_record(np.frombuffer(struct.pack("<qqqqd", 0, 2, 0, 0, 18.0), dtype="u1"))
            stream.write_record(np.array([-101.0, -100.5], dtype="<f8"))
            stream.write_record(np.zeros(2 * 298, dtype="<f8"))
        telemetry = {
            "requested": 2,
            "converged": 2,
            "reason": 1,
            "dimension": 4,
            "matrix_type": "mpisbaij",
            "solver_type": "krylovschur",
            "which": 1,
            "eigenpairs": [{"absolute_residual": 1e-10}] * 2,
        }
        (p / "eps-telemetry.jsonl").write_text(json.dumps(telemetry) + "\n")
        (p / "eigenp.inp").write_text("luphso=1021,")
        np.savetxt(p / "fort.1021", np.column_stack((0.1 + 0.05 * np.arange(99), np.zeros(99))))
        (p / "reson.out").write_text("")
        saved["sectors"][sector] = {
            "stages": [
                {"name": n, "original_exit": 0, "timed_out": False}
                for n in ("scatci", "rsolve", "eigenp", "reson")
            ]
        }
    (tmp_path / "pilot.json").write_text(json.dumps(saved))
    (tmp_path / "resources.json").write_text('{"exit_code":0}')
    contract = {
        "sectors": ["B1", "B2"],
        "dimension": 4,
        "continuum_coefficients": 2,
        "boundary_channels": 298,
        "raw_checks": {"native_selection_which": 1, "pi_phase_rad_mod_pi": 0.002},
    }
    result, _ = read_run(tmp_path, 2, contract)
    assert result["sectors"]["B1"]["export_checks_passed"]
    del contract["boundary_channels"]
    with pytest.raises(ValueError, match="boundary export"):
        read_run(tmp_path, 2, contract)
