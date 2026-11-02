"""Running launch placeholders cannot become terminal qualification verdicts."""

from copy import deepcopy

import pytest

from projects.ukrmol_co.verify_stretched_augtz_properties_launch import launch_disposition


def test_active_placeholders_stay_null_and_a_terminal_or_expanded_scope_rejects():
    owner = dict(state="running", original_owner_exit=2)
    science = dict(
        scientific_exit=2,
        property_inventory_verdict=None,
        native_denprop_calls=1,
        stages=[],
        qc_calls=0,
        target_diagonalizations=0,
        scattering_calls=0,
        new_maps=0,
        route_holdout_s_matrices_analyzed=0,
        full_sweep_released=False,
    )
    container = dict(State=dict(Running=True, OOMKilled=False))
    r = launch_disposition(owner, science, container)
    assert r["original_terminal_exits"] is r["property_inventory_verdict"] is None
    for changes in [
        dict(native_denprop_calls=2),
        dict(qc_calls=1),
        dict(new_maps=1),
        dict(route_holdout_s_matrices_analyzed=1),
        dict(property_inventory_verdict=True),
    ]:
        with pytest.raises(ValueError):
            launch_disposition(owner, science | changes, container)
    bad = deepcopy(container)
    bad["State"]["Running"] = False
    with pytest.raises(ValueError):
        launch_disposition(owner, science, bad)
    with pytest.raises(ValueError):
        launch_disposition(owner | dict(state="finished"), science, container)
