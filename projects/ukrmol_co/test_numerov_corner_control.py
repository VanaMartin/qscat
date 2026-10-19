"""All-endpoint/native veto and complete finite-grid contract checks."""

import mpmath as mp
import numpy as np
import pytest

from projects.ukrmol_co.native_reconstruction_replay import scattering
from projects.ukrmol_co.numerov_corner_control import configurations, evaluate
from projects.ukrmol_co.numerov_corner_replay import decisions


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


def test_independent_saved_expressions_preserve_the_new_endpoint_veto():
    rows, reference = controls()
    for row in rows:
        if row["endpoint_bohr"] == 800:
            row["k"][0][0] += 1.5e-6
        s = scattering(np.asarray(row["k"]))
        row.update(s_real=s.real.tolist(), s_imag=s.imag.tolist())
    replay = decisions(rows, reference, [0.025, 0.0125, 0.00625], 1e-6)
    producer = evaluate(rows, reference, [0.025, 0.0125, 0.00625], 1e-6)
    assert replay["checks"] == producer["checks"]
    assert not replay["checks"]["native_800_k"]
    assert not replay["verdict"]
    rows[-1]["s_real"][0][0] += 1e-8
    with pytest.raises(ValueError, match="Stored S differs"):
        decisions(rows, reference, [0.025, 0.0125, 0.00625], 1e-6)


def test_free_numerov_phase_defect_has_the_derived_fifth_order_coefficient():
    with mp.workdps(80):
        errors = []
        for text in ["0.01", "0.005", "0.0025"]:
            t = mp.mpf(text)
            theta = mp.acos((1 - 5 * t**2 / 12) / (1 + t**2 / 12))
            coefficient = (theta - t) / t**5
            errors.append(abs(coefficient * 480 - 1))
        assert errors[2] < mp.mpf("1e-6")
        assert errors[1] < errors[0] / 3
        assert errors[2] < errors[1] / 3
