"""Analytic/synthetic controls only; no ODE before the frozen allocation."""

import ast
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np

from projects.ukrmol_co import coupled_jost as control
from projects.ukrmol_co.coupled_jost_owner import command


def test_independent_fundamental_seed_agrees_with_rank_one_solution():
    k = np.diag([1, np.sqrt(3)])
    c = k + 0.4 * np.ones((2, 2))
    for z in [1, 3 + 2j, 5 + 4j]:
        exponents = np.diag(k) * z
        u = np.diag(np.cosh(exponents)) + np.diag(np.sinh(exponents) / np.diag(k)) @ c
        du = np.diag(np.diag(k) * np.sinh(exponents)) + np.diag(np.cosh(exponents)) @ c
        w = np.linalg.solve(u.T, du.T).T
        np.testing.assert_allclose(control.superpotential(z, 0.4)[0], w, atol=1e-13, rtol=1e-13)


def test_exact_real_flux_and_wronskian_without_projection():
    y, a = control.exact_maps(1, 0.4)
    result = control.metrics(y, a, (y, a))
    np.testing.assert_allclose(
        [result["flux_defect"], result["wronskian_residual"]], 0, atol=1e-13, rtol=0
    )
    assert abs(y[0, 1]) > 1e-3 and abs(a[0, 1]) > 1e-3


def test_exact_coupled_maps_satisfy_transport_equations_by_finite_derivatives():
    for z in [1, 2 + 1j]:
        h = 1e-5
        y, a = control.exact_maps(z, 0.4)
        plus, minus = control.exact_maps(z + h, 0.4), control.exact_maps(z - h, 0.4)
        np.testing.assert_allclose(
            (plus[0] - minus[0]) / (2 * h),
            control.potential(z, 0.4) - np.eye(2) - y @ y,
            atol=1e-9,
            rtol=1e-7,
        )
        np.testing.assert_allclose((plus[1] - minus[1]) / (2 * h), -a @ y, atol=1e-9, rtol=1e-7)


def test_genuine_real_symmetric_noncommuting_potential_and_derivative():
    z, h = 1.0, 1e-5
    v = control.potential(z, 0.4)
    np.testing.assert_allclose(v, v.T, atol=1e-14, rtol=0)
    assert np.isrealobj(v) and abs(v[0, 1]) > 1e-3
    finite = (control.superpotential(z + h, 0.4)[0] - control.superpotential(z - h, 0.4)[0]) / (
        2 * h
    )
    np.testing.assert_allclose(control.superpotential(z, 0.4)[1], finite, atol=1e-9, rtol=1e-7)
    assert abs((v @ np.diag([0, 2]) - np.diag([0, 2]) @ v)[0, 1]) > 1e-3


def test_first_scientific_failure_is_preserved_without_retry():
    contract = json.loads(
        Path(__file__).with_name("coupled-jost-contract-20261102.json").read_text()
    )
    with TemporaryDirectory() as tmp:
        with patch.object(control, "integrate", side_effect=ValueError("synthetic stop")) as call:
            assert control.run(contract, Path(tmp) / "outputs") == 2
        assert call.call_count == 1
        result = json.loads((Path(tmp) / "outputs/result.json").read_text())
        assert result["coupled_control_verdict"] is None and result["co_maps"] == 0


def test_owner_mounts_only_minimal_source_and_output_and_has_valid_controls_json():
    args = command(Path("/frozen"), Path("/output"), dict(image="sha256:pinned"))
    assert [args[i + 1] for i, item in enumerate(args) if item == "-v"] == [
        "/frozen:/source:ro",
        "/output:/out",
    ]
    assert "--network=none" in args and "--cpuset-cpus=0-3" in args and "--memory=2g" in args
    compile(args[-1], "<coupled-owner>", "exec")
    tree = ast.parse(args[-1])
    write = next(
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and n.func.attr == "write_text"
    )
    assert ast.literal_eval(write.args[0].right) == "\n"


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(fn)
        for name, fn in list(globals().items())
        if name.startswith("test_") and callable(fn)
    )
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
