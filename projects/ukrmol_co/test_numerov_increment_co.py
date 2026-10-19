"""Separate candidate K vetoes and exact family coverage before CO submission."""

import numpy as np
import pytest

from projects.ukrmol_co.numerov_increment_co import evaluate


def rows():
    return [
        dict(step_bohr=h, endpoint_bohr=r, tail_order=8, k=np.zeros((6, 6)).tolist())
        for h, r in [(0.025, 400), (0.0125, 400), (0.00625, 400), (0.0125, 200), (0.0125, 800)]
    ]


@pytest.mark.parametrize("family,index", [("step", 2), ("endpoint", 4)])
def test_large_k_change_cannot_be_waived_by_small_s_and_phase(family, index):
    controls = rows()
    reference = np.diag([-185.0, 0, 0, 0, 0, 0])
    for row in controls:
        row["k"] = reference.tolist()
    controls[index]["k"][0][0] += 2e-6
    answer = evaluate(controls, reference, 1e-6)
    assert answer["verdict"] is False
    assert answer["checks"][family + "_k"] is False
    assert answer["checks"][family + "_s"] is True
    assert answer["checks"][family + "_phase_rad_mod_pi"] is True


def test_native_reference_k_is_a_separate_veto():
    controls = rows()
    reference = np.diag([-185.0, 0, 0, 0, 0, 0])
    changed = reference.copy()
    changed[0, 0] += 2e-6
    for row in controls:
        row["k"] = changed.tolist()
    answer = evaluate(controls, reference, 1e-6)
    assert not answer["checks"]["native_k"]
    assert answer["checks"]["step_k"] and answer["checks"]["endpoint_k"]
    assert answer["checks"]["native_s"] and answer["checks"]["native_phase_rad_mod_pi"]


@pytest.mark.parametrize("mutate", ["missing", "duplicate", "undeclared"])
def test_missing_duplicate_and_undeclared_grids_reject(mutate):
    controls = rows()
    if mutate == "missing":
        controls.pop()
    elif mutate == "duplicate":
        controls[-1] = controls[0]
    else:
        controls[-1]["endpoint_bohr"] = 1600
    with pytest.raises(ValueError, match="Missing, duplicate or undeclared"):
        evaluate(controls, np.zeros((6, 6)), 1e-6)
