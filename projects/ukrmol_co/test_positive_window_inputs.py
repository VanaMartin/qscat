"""Explicit narrow six-open windows; generated inputs only, no CO or S fixture."""

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from projects.ukrmol_co import compressed_a2_pole_contrast as dispatch
from projects.ukrmol_co import compressed_pole_reference, local_pole_inputs
from projects.ukrmol_co.test_closed_pole_inputs import mocked_load, shallow_fixture


def contract():
    return dict(
        case="stretched/nine/B1",
        model_view="local-cas11-six-open-window-v1",
        positive_energy_interval_hartree=[0.001, 0.04],
        search_box_incident_hartree=[0.002, 0.03, -0.01, -1e-8],
    )


def test_positive_window_accepts_only_declared_six_open_domain():
    info, raw = shallow_fixture()
    model = mocked_load(
        info, raw, lambda p: local_pole_inputs.load_positive_boundary_model(p, [0.001, 0.04])
    )
    closed = mocked_load(info, raw, local_pole_inputs.load_closed_boundary_model)
    for key in ["poles", "amplitudes", "c2", "c3", "thresholds"]:
        np.testing.assert_allclose(model[key], closed[key], atol=0, rtol=0)
    with pytest.raises(ValueError, match="six-open-channel"):
        mocked_load(info, raw, local_pole_inputs.load_boundary_model)


@pytest.mark.parametrize(
    "interval",
    [
        [0, 0.04],
        [-0.01, 0.04],
        [0.04, 0.001],
        [0.01, np.inf],
        [np.nan, 0.04],
        [0.01j, 0.04],
        [0.01, 0.02, 0.03],
        [0.01, 0.05],
        [0.01, 0.06],
    ],
)
def test_invalid_or_threshold_touching_positive_windows_are_rejected(interval):
    info, raw = shallow_fixture()
    with pytest.raises(ValueError):
        mocked_load(
            info, raw, lambda p: local_pole_inputs.load_positive_boundary_model(p, interval)
        )


def test_windowed_dispatch_is_explicit_and_bounds_complex_search_real_projection():
    c = contract()
    info, raw = shallow_fixture()
    with (
        patch.object(local_pole_inputs, "read_interface", return_value=info),
        patch.object(local_pole_inputs, "records", return_value=raw),
    ):
        model = dispatch.load_pole_model(Path("raw/stretched/nine/B1"), c)
        assert model["amplitudes"].shape == (2048, 298)
        for box in [[0, 0.03, -0.01, -1e-8], [0.002, 0.045, -0.01, -1e-8]]:
            with pytest.raises(ValueError, match="window"):
                dispatch.load_pole_model(
                    Path("raw/stretched/nine/B1"), c | {"search_box_incident_hartree": box}
                )
        with pytest.raises(ValueError, match="six-open-channel"):
            dispatch.load_pole_model(
                Path("raw/stretched/nine/B1"), c | {"model_view": "local-cas11-nine-2048-v1"}
            )


def test_single_sector_windowed_inventory_still_requires_real_guard(tmp_path):
    raw = tmp_path / "raw"
    case = raw / "stretched/nine/B1"
    case.mkdir(parents=True)
    manifest = {}
    for name in ["fort.10", "fort.21"]:
        b = b"synthetic-inventory-only"
        (case / name).write_bytes(b)
        manifest["stretched/nine/B1/" + name] = hashlib.sha256(b).hexdigest()
    (raw / "input-manifest.json").write_text(json.dumps(manifest))
    (tmp_path / "canary.json").write_text(
        json.dumps({"canary_verdict": True, "rows": [{"checks": {"synthetic": True}}] * 11})
    )
    c = contract() | {
        "raw_directory": "raw",
        "input_sha256": {},
        "canary_result": "canary.json",
        "claim_boundary": "synthetic guard-order test",
    }
    with (
        patch.object(
            compressed_pole_reference.subprocess, "run", return_value=SimpleNamespace(returncode=0)
        ),
        patch.object(
            compressed_pole_reference,
            "require_real_guard",
            side_effect=ValueError("mandatory real guard"),
        ) as guard,
        patch.object(compressed_pole_reference, "load_pole_model") as load,
    ):
        assert compressed_pole_reference.run(tmp_path, tmp_path / "out", c) == 2
        guard.assert_called_once_with(tmp_path, c)
        load.assert_not_called()
    report = json.loads((tmp_path / "out/result.json").read_text())
    assert report["inputs_verified_before"] == 2
    assert "mandatory real guard" in report["error"]


def test_windowed_real_guard_rejects_holdout_index_before_use(tmp_path):
    c = contract() | {"real_development_indices": [4, 20]}
    payload = b"{}"
    (tmp_path / "guard.json").write_bytes(payload)
    c.update(
        real_guard_result="guard.json",
        input_sha256={"guard.json": hashlib.sha256(payload).hexdigest()},
    )
    with pytest.raises(ValueError, match="non-holdout"):
        dispatch.require_real_guard(tmp_path, c)
