"""Regression for the actual model-return contract rather than compressed mocks."""

import subprocess

from projects.ukrmol_co.verify_stretched_dz_pole_stop import SOURCE_COMMIT, returned_dict_keys


def test_original_own_model_has_no_nested_interface_metadata():
    raw = subprocess.check_output(
        ["git", "show", SOURCE_COMMIT + ":projects/ukrmol_co/native_full_reconstruction.py"]
    )
    keys = returned_dict_keys(raw, "_model_from_interface")
    assert "interface" not in keys and {"poles", "amplitudes", "thresholds", "angular"} <= set(keys)
    adapter = subprocess.check_output(
        ["git", "show", SOURCE_COMMIT + ":projects/ukrmol_co/stretched_dz_pole.py"]
    )
    assert b'model["interface"]["channel_labels"]' in adapter


def test_positive_schema_fixture_cannot_fabricate_missing_keys():
    assert returned_dict_keys(
        "def loader():\n return {'poles': [], 'angular': []}\n", "loader"
    ) == ["poles", "angular"]
