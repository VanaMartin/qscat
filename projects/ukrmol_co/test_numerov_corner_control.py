"""All-endpoint/native veto and complete finite-grid contract checks."""

import numpy as np
import pytest

from projects.ukrmol_co.numerov_corner_control import configurations, evaluate


def controls():
    reference = np.diag([-185.0, 0, 0, 0, 0, 0])
    rows = [
        dict(step_bohr=h, endpoint_bohr=r, tail_order=p, k=reference.tolist())
        for h, r, p in configurations([0.025, 0.0125, 0.00625])
    ]
    return rows, reference


def test_endpoint_native_veto_survives_passing_pairwise_endpoint_change():
    rows, reference = controls()
    for row in rows:
        if row["step_bohr"] == 0.00625:
            row["k"][0][0] += {200: 0.9e-6, 400: 0.9e-6, 800: 1.5e-6}[row["endpoint_bohr"]]
    result = evaluate(rows, reference, [0.025, 0.0125, 0.00625], 1e-6)
    assert result["checks"]["endpoint_k"]
    assert result["checks"]["native_400_k"]
    assert not result["checks"]["native_800_k"]
    assert result["checks"]["native_800_s"]
    assert not result["verdict"]


@pytest.mark.parametrize("family,index", [("step", 2), ("order", 5)])
def test_k_veto_is_separate_from_s_phase_and_invariants(family, index):
    rows, reference = controls()
    rows[index]["k"][0][0] += 2e-6
    result = evaluate(rows, reference, [0.025, 0.0125, 0.00625], 1e-6)
    assert not result["checks"][family + "_k"]
    assert result["checks"][family + "_s"]
    assert result["checks"][family + "_phase_rad_mod_pi"]
    assert result["checks"]["unitarity"] and result["checks"]["reciprocity"]
    assert not result["verdict"]


@pytest.mark.parametrize("change", ["missing", "duplicate", "new_grid"])
def test_exact_configuration_coverage_is_required(change):
    rows, reference = controls()
    if change == "missing":
        rows.pop()
    elif change == "duplicate":
        rows[-1] = rows[0]
    else:
        rows[-1]["endpoint_bohr"] = 1600
    with pytest.raises(ValueError, match="Missing, duplicate or undeclared"):
        evaluate(rows, reference, [0.025, 0.0125, 0.00625], 1e-6)
