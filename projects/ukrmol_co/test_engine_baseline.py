"""A same-root control must reject drift without claiming spectral convergence."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co.engine_baseline import (
    audit_baseline,
    audit_native_inputs,
    audit_root_sequence,
    label_engine_reproduction,
    main,
)
from projects.ukrmol_co.engine_baseline_owner import verify_intermediate
from projects.ukrmol_co.sparse_comparison import compare_phases
from projects.ukrmol_co.sparse_comparison_owner import preflight


def settings():
    return json.loads(
        Path(__file__).with_name("cas12-audited-engine-baseline-contract.json").read_text()
    )


def phases(position_ev=2.457, branch=0):
    contract = settings()
    requested = 0.1 + 0.05 * np.arange(99)
    energy = requested * contract["native_rydberg_per_requested_ev"] / 2
    center = position_ev / contract["hartree_to_ev"]
    width = 0.5 / contract["hartree_to_ev"]
    s = (center - energy + 0.5j * width) / (center - energy - 0.5j * width)
    np.testing.assert_allclose(abs(s), 1, atol=3e-15, rtol=0)
    return np.column_stack((requested, np.angle(s) / 2 + branch * np.pi))


@pytest.mark.parametrize("shift,expected", [(0, True), (0.001, False)])
def test_engine_reproduction_is_tighter_than_scientific_pair_budget(shift, expected):
    row = compare_phases(phases(), phases(position_ev=2.457 + shift, branch=2), settings())
    result = label_engine_reproduction(
        {
            "root_pair_stability_verdict": row["root_pair_stability_verdict"],
            "spectral_truncation_verdict": None,
            "sectors": {name: copy.deepcopy(row) for name in ("B1", "B2")},
        }
    )
    assert result["matched_engine_baseline_verdict"] is expected
    assert result["root_pair_stability_verdict"] is None
    assert result["spectral_truncation_verdict"] is None
    assert all(r["root_pair_stability_verdict"] is None for r in result["sectors"].values())


def test_invalid_fit_stays_unqualified_and_changed_root_count_is_rejected(tmp_path):
    flat = phases()
    flat[:, 1] = 0
    row = compare_phases(flat, flat, settings())
    result = label_engine_reproduction(
        {
            "root_pair_stability_verdict": row["root_pair_stability_verdict"],
            "sectors": {"B1": row},
        }
    )
    assert result["matched_engine_baseline_verdict"] is None
    with pytest.raises(ValueError, match="2048 roots on both sides"):
        audit_baseline(tmp_path, tmp_path, settings() | {"candidate_roots": 4096})


def test_root_pair_rejects_changed_integrals_and_non_root_deck_changes(tmp_path):
    left, right = tmp_path / "reference", tmp_path / "candidate"
    contract = settings() | {"candidate_roots": 4096, "sectors": ["B1"]}
    for directory, roots in ((left, 2048), (right, 4096)):
        sector = directory / "B1"
        sector.mkdir(parents=True)
        (sector / "scatci.inp").write_text(f"nstat = {roots}, megul = 721, crite = 1e-12")
        for name in ("fort.16", "fort.24", "fort.26", "fort.721"):
            (sector / name).write_bytes(b"same saved native input")
    audit_native_inputs(left, right, contract, same_roots=False)
    (right / "B1/fort.721").write_bytes(b"altered integral")
    with pytest.raises(ValueError, match=r"fort\.721"):
        audit_native_inputs(left, right, contract, same_roots=False)
    (right / "B1/scatci.inp").write_text("nstat = 4096, megul = 721, crite = 1e-8")
    with pytest.raises(ValueError, match="beyond root count"):
        audit_native_inputs(left, right, contract, same_roots=False)
    (right / "B1/scatci.inp").write_text("nstat = 2048, megul = 721, crite = 1e-12")
    with pytest.raises(ValueError, match="requested root count"):
        audit_native_inputs(left, right, contract, same_roots=False)


@pytest.mark.parametrize("verdict", [False, None])
@pytest.mark.parametrize("action", ["root-compare", "sequence-compare"])
def test_failed_root_reader_preserves_result_and_exits_nonzero(
    tmp_path, monkeypatch, verdict, action
):
    from projects.ukrmol_co import engine_baseline

    contract = tmp_path / "contract.json"
    contract.write_text(json.dumps(settings() | {"candidate_roots": 4096}))
    output = tmp_path / "result.json"
    key = (
        "root_sequence_stability_verdict"
        if action == "sequence-compare"
        else "root_pair_stability_verdict"
    )
    reader = "audit_root_sequence" if action == "sequence-compare" else "audit_root_pair"
    monkeypatch.setattr(engine_baseline, reader, lambda *args: {key: verdict})
    monkeypatch.setattr(
        "sys.argv",
        [
            "reader",
            action,
            "--reference",
            str(tmp_path),
            "--candidate",
            str(tmp_path),
            "--contract",
            str(contract),
            "--output",
            str(output),
            "--intermediate",
            str(tmp_path),
        ],
    )
    with pytest.raises(SystemExit) as stopped:
        main()
    assert stopped.value.code == 1
    assert json.loads(output.read_text())[key] is verdict


@pytest.mark.parametrize("owner_exit,verdict", [(0, False), (0, None), (1, True)])
def test_root_owner_requires_completed_passing_engine_baseline(tmp_path, owner_exit, verdict):
    baseline = tmp_path / "baseline/native/run"
    baseline.mkdir(parents=True)
    (baseline.parent.parent / "execution.json").write_text(
        json.dumps({"original_controller_exit": owner_exit})
    )
    (baseline / "independent-engine-baseline.json").write_text(
        json.dumps({"matched_engine_baseline_verdict": verdict})
    )
    with pytest.raises(ValueError, match="Baseline owner/reader gate failed"):
        preflight(
            tmp_path,
            tmp_path,
            {
                "baseline_directory": "baseline/native/run",
                "baseline_reader_filename": "independent-engine-baseline.json",
                "baseline_reader_verdict": "matched_engine_baseline_verdict",
            },
            {},
        )


def test_sequence_requires_endpoints_even_when_both_successive_pairs_pass(tmp_path, monkeypatch):
    from projects.ukrmol_co import engine_baseline

    contract = settings() | {
        "baseline_roots": 2048,
        "intermediate_roots": 4096,
        "candidate_roots": 8192,
    }
    contract["observable_budgets"] = contract["observable_budgets"] | {
        "phase_rad_mod_pi": 0.05,
        "position_ev": 0.05,
        "full_width_relative": 0.05,
        "full_width_floor_ev": 0.001,
    }
    spectra = {}
    for roots, shift in ((2048, 0), (4096, 0.03), (8192, 0.06)):
        spectra[roots] = phases()
        spectra[roots][:, 1] += shift

    def analytic_pair(reference, candidate, settings):
        return compare_phases(
            spectra[settings["baseline_roots"]], spectra[settings["candidate_roots"]], settings
        ) | {"spectral_truncation_verdict": None}

    monkeypatch.setattr(engine_baseline, "audit_root_pair", analytic_pair)
    result = audit_root_sequence(tmp_path, tmp_path, tmp_path, contract)
    assert all(row["root_pair_stability_verdict"] for row in result["successive_pairs"].values())
    assert result["root_pair_stability_verdict"] is False
    assert result["root_sequence_stability_verdict"] is False
    assert result["spectral_truncation_verdict"] is None


@pytest.mark.parametrize("verdict", [True, False, None])
def test_sequence_preserves_unset_pairs_and_omission_scope(tmp_path, monkeypatch, verdict):
    from projects.ukrmol_co import engine_baseline

    def pair(reference, candidate, contract):
        return {
            "root_pair_stability_verdict": verdict if contract["baseline_roots"] == 4096 else True,
            "spectral_truncation_verdict": None,
        }

    monkeypatch.setattr(engine_baseline, "audit_root_pair", pair)
    result = audit_root_sequence(
        tmp_path,
        tmp_path,
        tmp_path,
        {"baseline_roots": 2048, "intermediate_roots": 4096, "candidate_roots": 8192},
    )
    assert result["root_sequence_stability_verdict"] is verdict
    assert result["spectral_truncation_verdict"] is None


@pytest.mark.parametrize(
    "owner_exit,verdict,expected",
    [(0, True, True), (0, False, False), (0, None, False), (1, True, False)],
)
def test_third_resolution_requires_completed_intermediate(tmp_path, owner_exit, verdict, expected):
    directory = tmp_path / "control/native/run"
    directory.mkdir(parents=True)
    (directory.parent.parent / "execution.json").write_text(
        json.dumps({"original_controller_exit": owner_exit, "finished_utc": "completed"})
    )
    (directory / "independent-comparison.json").write_text(
        json.dumps({"root_pair_stability_verdict": verdict, "independent_raw_reading_passed": True})
    )
    contract = {"intermediate_directory": "control/native/run"}
    if expected:
        assert "comparison_sha256" in verify_intermediate(tmp_path, contract)
    else:
        with pytest.raises(ValueError, match="Intermediate owner/root-pair prerequisite failed"):
            verify_intermediate(tmp_path, contract)
