"""Fixed identities and independent-reader API never imply complex release."""

from projects.ukrmol_co.verify_stretched_real_completion import MAPS


def test_real_completion_reuses_dz80_and_has_exact_tz_endpoint_triple():
    assert len(MAPS) == 5 and all(not (case == "dz" and index == 80) for case, index, _ in MAPS)
    assert [endpoint for case, index, endpoint in MAPS if case == "tz" and index == 80] == [
        60,
        90,
        120,
    ]
