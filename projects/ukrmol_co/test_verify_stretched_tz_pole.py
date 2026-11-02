"""The shared witness algebra keeps distinct own-TZ source and native identity."""

from projects.ukrmol_co.verify_stretched_tz_pole import ALLOCATION


def test_tz_provenance_and_control_count_are_not_the_dz_release():
    assert ALLOCATION["case"] == "tz" and ALLOCATION["source_count"] == 52
    assert ALLOCATION["runtime_controls"] == 31
    assert ALLOCATION["execution_path"].endswith("stretched-tz-pole-execution-20261102.json")
