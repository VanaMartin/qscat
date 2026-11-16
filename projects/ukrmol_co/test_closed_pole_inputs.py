"""Generated below-threshold loader controls; no molecular or K/S fixture."""

from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from scipy.io import FortranFile

from projects.ukrmol_co import local_pole_inputs
from projects.ukrmol_co.test_local_pole_inputs import fixture


def shallow_fixture():
    info, raw = fixture()
    info["target_energies_hartree"][1] = -111.95
    info["channel_thresholds_rydberg"][6:] = [0.1] * 292
    return info, raw


def mocked_load(info, raw, loader):
    with (
        patch.object(local_pole_inputs, "read_interface", return_value=info),
        patch.object(local_pole_inputs, "records", return_value=raw) as read,
    ):
        model = loader(Path("no-k-or-s-file"))
    read.assert_called_once_with(Path("no-k-or-s-file/fort.21"))
    return model


def test_closed_loader_accepts_shallow_threshold_without_relaxing_legacy_grid():
    info, raw = shallow_fixture()
    before = deepcopy(info)
    model = mocked_load(info, raw, local_pole_inputs.load_closed_boundary_model)
    np.testing.assert_allclose(model["thresholds"][6:], 0.1, atol=1e-14, rtol=0)
    np.testing.assert_allclose(model["c2"][:2, :2], [[4, 3], [3, 7]], atol=0, rtol=0)
    np.testing.assert_allclose(model["c3"][:2, :2], [[7, 11], [11, 13]], atol=0, rtol=0)
    assert info == before
    with pytest.raises(ValueError, match="six-open-channel"):
        mocked_load(info, raw, local_pole_inputs.load_boundary_model)


def test_closed_and_legacy_decoders_preserve_identical_native_arrays():
    info, raw = fixture()
    amplitudes = np.zeros((2048, 298), dtype="<f8")
    amplitudes[0, :2] = [0.25, -0.5]
    amplitudes[-1, -1] = 0.75
    raw[-1] = amplitudes.tobytes()
    closed = mocked_load(info, raw, local_pole_inputs.load_closed_boundary_model)
    legacy = mocked_load(info, raw, local_pole_inputs.load_boundary_model)
    assert closed.keys() == legacy.keys()
    for key in ["thresholds", "angular", "poles", "amplitudes", "c2", "c3"]:
        np.testing.assert_allclose(closed[key], legacy[key], atol=0, rtol=0)
    np.testing.assert_allclose(closed["amplitudes"], amplitudes, atol=0, rtol=0)
    np.testing.assert_allclose(closed["poles"], np.arange(2048), atol=0, rtol=0)


@pytest.mark.parametrize("threshold", [0.0, -0.1, np.nan, np.inf])
def test_closed_loader_rejects_nonpositive_or_nonfinite_excited_threshold(threshold):
    info, raw = fixture()
    info["target_energies_hartree"][1] = -112 + threshold / 2
    info["channel_thresholds_rydberg"][6:] = [threshold] * 292
    with pytest.raises((ValueError, AssertionError)):
        mocked_load(info, raw, local_pole_inputs.load_closed_boundary_model)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("target_components", 47),
        ("boundary_channels", 297),
        ("boundary_poles", 4096),
        ("multipole_orders", 1),
        ("radius_bohr", 17),
        ("symmetry_charge", [1, 2, 0, -1]),
        ("channel_thresholds_rydberg", [0.01] + [0] * 5 + [0.5] * 292),
        ("channel_thresholds_rydberg", [0] * 6 + [0.4] * 292),
    ],
)
def test_closed_loader_retains_layout_and_native_threshold_identity_guards(key, value):
    info, raw = fixture()
    info[key] = value
    with pytest.raises((ValueError, AssertionError)):
        mocked_load(info, raw, local_pole_inputs.load_closed_boundary_model)


def test_closed_loader_reads_complete_generated_native_records(tmp_path):
    info, raw = shallow_fixture()

    def integers(values):
        return np.asarray(values, dtype="<i8").tobytes()

    def floats(values):
        return np.asarray(values, dtype="<f8").tobytes()

    dimensions = integers([48, 0, 0, 298])
    geometry = integers([1, 2, 0, 0]) + floats([2.8, 6.0])
    channels = [integers([10, 1, 0, 0, 0]), b"synthetic".ljust(80), dimensions, geometry]
    channels.extend(
        integers([i + 1, 0, 0, 0]) + floats([energy])
        for i, energy in enumerate(info["target_energies_hartree"])
    )
    channels.extend(
        integers(label) + floats([threshold])
        for label, threshold in zip(
            info["channel_labels"], info["channel_thresholds_rydberg"], strict=True
        )
    )
    boundary = [
        integers([11, 1, 0, 0, 0]),
        b"synthetic".ljust(80),
        dimensions,
        geometry,
        integers([2, 2048, 0, 0]) + floats([18]),
        *raw[5:],
    ]
    for name, rows in [("fort.10", channels), ("fort.21", boundary)]:
        with FortranFile(tmp_path / name, "w", header_dtype="<u4") as stream:
            for row in rows:
                stream.write_record(np.frombuffer(row, dtype="u1"))
    model = local_pole_inputs.load_closed_boundary_model(tmp_path)
    np.testing.assert_allclose(model["thresholds"][6:], 0.1, atol=1e-14, rtol=0)
    np.testing.assert_allclose(model["c2"][:2, :2], [[4, 3], [3, 7]], atol=0, rtol=0)
    assert model["amplitudes"].shape == (2048, 298)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["fort.10", "fort.21"]
