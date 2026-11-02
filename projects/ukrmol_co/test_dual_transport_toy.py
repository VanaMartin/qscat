"""Analytic/control-plane checks only; no ODE or CO propagation before the frozen toy."""

import ast
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import dual_transport_toy as toy
from projects.ukrmol_co.dual_transport_toy_owner import container_command


def test_zero_injection_exact_row_and_identity_flux():
    a = toy.exact_row(1, 1, 40, 20, 0)
    np.testing.assert_allclose(a, [1, 0], atol=1e-15, rtol=0)
    r = toy.observables(a, 1, 1)
    np.testing.assert_allclose([r["s_identity_error"], r["flux_defect"]], 0, atol=1e-15, rtol=0)


def test_critical_angle_keeps_relative_error_magnitude():
    a = toy.exact_row(1, 1, 45, 20, 1e-6)
    np.testing.assert_allclose(abs(a[1]), 1e-6, atol=1e-20, rtol=1e-12)
    assert abs(toy.exact_row(1, 1, 20, 20, 1e-6)[1]) < 1e-6
    assert abs(toy.exact_row(1, 1, 50, 20, 1e-6)[1]) > 1e-6


def test_injected_error_can_pass_s_while_failing_flux_without_any_integration():
    r = toy.observables(toy.exact_row(1, 1, 50, 40, 1e-6), 1, 1)
    assert r["s_identity_error"] < 1e-6 and r["flux_defect"] > 1e-7


def test_failed_first_integration_is_preserved_and_not_retried():
    c = json.loads(
        Path(__file__).with_name("dual-transport-toy-contract-20261102.json").read_text()
    )
    with TemporaryDirectory() as tmp:
        with patch.object(toy, "integrate", side_effect=ValueError("synthetic failure")) as call:
            assert toy.run(c, Path(tmp) / "output") == 2
        assert call.call_count == 1
        r = json.loads((Path(tmp) / "output/result.json").read_text())
        assert r["toy_verdict"] is None and r["co_maps"] == r["native_calls"] == 0


def test_owner_has_no_campaign_or_native_input_mount():
    command = container_command(Path("/frozen"), Path("/outputs"), dict(image="sha256:pinned"))
    assert (
        "--cpuset-cpus=0-3" in command and "--memory=2g" in command and "--network=none" in command
    )
    assert [command[i + 1] for i, v in enumerate(command) if v == "-v"] == [
        "/frozen:/source:ro",
        "/outputs:/out",
    ]
    assert "/work" not in command[-1] and "stretched_real_guard" not in command[-1]
    compile(command[-1], "<toy-owner>", "exec")
    tree = ast.parse(command[-1])
    write = next(
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and n.func.attr == "write_text"
    )
    assert isinstance(write.args[0], ast.BinOp) and ast.literal_eval(write.args[0].right) == "\n"


def test_frozen_case_count_and_all_resolutions_without_an_ode():
    c = json.loads(
        Path(__file__).with_name("dual-transport-toy-contract-20261102.json").read_text()
    )

    def analytic_instead(k, kappa, angle, length, epsilon, rtol, contract):
        return toy.exact_row(k, kappa, angle, length, epsilon), 0

    with TemporaryDirectory() as tmp:
        with patch.object(toy, "integrate", side_effect=analytic_instead) as call:
            assert toy.run(c, Path(tmp) / "output") == 0
        r = json.loads((Path(tmp) / "output/result.json").read_text())
        assert call.call_count == len(r["rows"]) == 72
        assert r["negative_control_cases"] > 0
        assert {row["rtol"] for row in r["rows"]} == set(c["rtols"])


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
