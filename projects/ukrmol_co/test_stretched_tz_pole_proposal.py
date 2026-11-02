"""Prospective own-TZ seed, sheet and inherited-budget checks; no CO maps."""

import hashlib
import json
from pathlib import Path

import numpy as np


def test_tz_seed_is_own_native_not_the_dz_pole_and_exposes_the_whole_box():
    source = Path(__file__).resolve().parents[2]
    tz = json.loads(
        Path(__file__).with_name("stretched-tz-pole-proposal-20261102.json").read_text()
    )
    raw = (source / tz["base_dz_proposal"]).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == tz["base_dz_proposal_sha256"]
    dz = json.loads(raw)
    box = tz["search_box_incident_hartree"]
    assert (
        tz["near_feature_real_index"] == 18
        and tz["seeds_incident_hartree"] != dz["seeds_incident_hartree"]
    )
    assert np.arctan(abs(box[2]) / box[0]) / 2 < np.deg2rad(min(s[1] for s in dz["settings"]))
    assert 2 * box[1] < 0.370734394
    assert tz["execution_authorized"] is tz["pole_search_released"] is False
    assert dz["resources"]["total_map_cap"] == 79 and tz["native_calls"] == 0
