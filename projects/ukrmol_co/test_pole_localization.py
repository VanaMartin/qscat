"""Known-truth globalization controls; no CO records or exterior propagation."""

import unittest
from itertools import pairwise

import numpy as np

from projects.ukrmol_co.pole_localization import configured_search, safeguarded_search
from projects.ukrmol_co.pole_reference import bordered_search
from projects.ukrmol_co.test_pole_reference import circle, exact_pole, level


def scalar(function):
    return lambda z: dict(q=np.eye(1), y=np.array([[1 - function(z)]]))


def test_historical_default_and_explicit_dispatch_do_not_silently_change_search():
    assert configured_search({}) == (bordered_search, {})
    controls = dict(max_step=0.025, max_backtracks=8, armijo=1e-4)
    contract = dict(root_search="safeguarded-newton-v1", safeguard_controls=controls)
    search, actual = configured_search(contract)
    assert search is safeguarded_search and actual == controls and actual is not controls
    for bad in [
        dict(root_search="unknown"),
        dict(safeguard_controls=controls),
        dict(root_search="safeguarded-newton-v1"),
        dict(root_search="safeguarded-newton-v1", safeguard_controls=dict(controls, observer=None)),
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            configured_search(bad)


def test_one_level_nonlocal_seed_and_probe_refinement_keep_the_exact_pole():
    roots = []
    for h in [2e-3, 2e-4, 2e-5]:
        z, report = safeguarded_search(
            level, [0.11 - 0.085j, 0.11 + h - 0.085j], [0.1, 0.25, -0.1, -0.01], 64
        )
        roots.append(z)
        np.testing.assert_allclose(z, exact_pole(), atol=1e-10, rtol=0)
        assert report["evaluations"] <= 64
        witness, _ = circle(z)
        assert all(witness["checks"].values())
    np.testing.assert_allclose(roots, [exact_pole()] * 3, atol=1e-10, rtol=0)


def test_backtracking_rejects_a_known_overshooting_newton_step():
    root = 1 - 0.04j

    def function(z):
        return (z - root) ** 2 + 2 * (z - root)

    seed = root - 0.8
    z, report = safeguarded_search(
        scalar(function), [seed, seed + 1e-5], [-0.2, 3, -0.1, -0.001], 80, max_step=2
    )
    np.testing.assert_allclose(z, root, atol=1e-10, rtol=0)
    assert any(row["reason"] == "insufficient decrease" for row in report["trials"])
    assert any(row["accepted"] and row["backtrack"] > 0 for row in report["trials"])
    centers = report["centers"]
    assert all(b["alpha_abs"] < a["alpha_abs"] for a, b in pairwise(centers))


def test_outside_proposals_never_call_the_model_and_cannot_manufacture_a_root():
    calls = []

    def model(z):
        assert 0 <= z.real <= 1 and -0.1 <= z.imag <= -0.001
        calls.append(z)
        return scalar(lambda x: x - (2 - 0.04j))(z)

    with unittest.TestCase().assertRaisesRegex(ValueError, "decreasing|cap"):
        safeguarded_search(
            model, [0.5 - 0.04j, 0.50001 - 0.04j], [0, 1, -0.1, -0.001], 40, max_step=2
        )
    assert len(calls) <= 40


def test_backward_probe_and_exact_start_do_not_need_an_extra_map():
    root = 0.99 - 0.04j
    model = scalar(lambda z: z - root)
    z, report = safeguarded_search(model, [0.995 - 0.04j, 1.005 - 0.04j], [0, 1, -0.1, -0.001], 8)
    np.testing.assert_allclose(z, root, atol=1e-10, rtol=0)
    assert report["probes"][0]["direction"] == "backward"
    z, report = safeguarded_search(model, [root, root + 1e-5], [0, 1, -0.1, -0.001], 1)
    assert report["evaluations"] == 1


def test_multiple_null_and_nonroot_constant_are_not_successes():
    case = unittest.TestCase()
    with case.assertRaises(ValueError):
        safeguarded_search(
            lambda z: dict(q=np.eye(2), y=np.eye(2)),
            [0.2 - 0.04j, 0.20001 - 0.04j],
            [0.1, 0.3, -0.1, -0.01],
            8,
        )
    with case.assertRaisesRegex(ValueError, "derivative"):
        safeguarded_search(
            scalar(lambda z: 1), [0.2 - 0.04j, 0.20001 - 0.04j], [0.1, 0.3, -0.1, -0.01], 8
        )


def test_hard_model_failure_is_not_hidden_by_backtracking():
    calls = []
    observed = []

    def broken(z):
        calls.append(z)
        if len(calls) == 3:
            raise TimeoutError("ODE control")
        return scalar(lambda x: x - (0.24 - 0.04j))(z)

    with unittest.TestCase().assertRaisesRegex(TimeoutError, "ODE control"):
        safeguarded_search(
            broken,
            [0.15 - 0.04j, 0.15001 - 0.04j],
            [0.1, 0.3, -0.1, -0.01],
            8,
            observer=lambda r: observed.append(r.copy()),
        )
    assert len(calls) == 3
    assert observed[-1]["evaluations"] == 3 and not observed[-1]["trace"][-1]["completed"]


def test_cancelling_features_and_two_simple_roots_do_not_imply_global_identity():
    root = 0.2 - 0.04j
    model = scalar(lambda z: (z - root) / (z - root))
    with unittest.TestCase().assertRaisesRegex(ValueError, "derivative|decreasing|cap"):
        safeguarded_search(model, [0.15 - 0.07j, 0.15001 - 0.07j], [0.1, 0.3, -0.1, -0.01], 16)
    other = 0.27 - 0.07j
    model = scalar(lambda z: (z - root) * (z - other))
    z, _ = safeguarded_search(model, [0.14 - 0.05j, 0.14001 - 0.05j], [0.1, 0.3, -0.1, -0.01], 64)
    np.testing.assert_allclose(z, root, atol=1e-9, rtol=0)
    # A local candidate is not a count of all roots: the second zero still exists.
    np.testing.assert_allclose(model(other)["y"], [[1]], atol=1e-15, rtol=0)


def test_invalid_controls_and_caps_stop_before_unbudgeted_calls():
    case = unittest.TestCase()
    model = scalar(lambda z: z - (0.24 - 0.04j))
    for options in [dict(limit=0), dict(max_step=0), dict(max_backtracks=-1), dict(armijo=1)]:
        with case.assertRaises(ValueError):
            kw = dict(limit=8) | options
            safeguarded_search(
                model, [0.15 - 0.04j, 0.15001 - 0.04j], [0.1, 0.3, -0.1, -0.01], **kw
            )
    with case.assertRaisesRegex(ValueError, "cap"):
        safeguarded_search(model, [0.15 - 0.04j, 0.15001 - 0.04j], [0.1, 0.3, -0.1, -0.01], 1)
    with case.assertRaises(ValueError):
        safeguarded_search(model, [0.15 - 0.04j, 0.15 - 0.03j], [0.1, 0.3, -0.1, -0.01], 8)
    with case.assertRaises(ValueError):
        safeguarded_search(
            lambda z: dict(q=np.eye(1), y=np.array([[np.nan]])),
            [0.15 - 0.04j, 0.15001 - 0.04j],
            [0.1, 0.3, -0.1, -0.01],
            8,
        )


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
