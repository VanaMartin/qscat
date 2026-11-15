"""Offline analytic controls and explicitly opted-in saved-data regression."""

import hashlib
import io
import json
import os
import tarfile
from pathlib import Path

import numpy as np
import pytest
from qscat.units import hartree_to_ev
from scipy.interpolate import CubicSpline

from projects.ukrmol_co import local_representation as local


@pytest.fixture
def analytic():
    knots = np.array([1.9, 2.0, 2.15, 2.3, 2.5])
    # Global cubics make not-a-knot interpolation exact, without fitting physics.
    er = 0.1 + 0.02 * knots - 0.003 * knots**2 + 0.004 * knots**3
    log_gamma = -3 + 0.1 * knots + 0.02 * knots**2 - 0.03 * knots**3
    return local.LocalRepresentation(
        neutral_coefficients=[0.2, 0.03, -0.01, 0.004, 0.002, -0.001, 0.0003],
        knots_bohr=knots,
        er_coefficients=CubicSpline(knots, er).c,
        log_gamma_coefficients=CubicSpline(knots, log_gamma).c,
    )


def test_analytic_curves_shape_units_and_width(analytic):
    radius = np.linspace(1.9, 2.5, 42).reshape(6, 7)
    x = (radius - 2.2) / 0.3
    # Independent explicit Chebyshev polynomials, not chebval.
    expected_v0 = (
        0.2
        + 0.03 * x
        - 0.01 * (2 * x**2 - 1)
        + 0.004 * (4 * x**3 - 3 * x)
        + 0.002 * (8 * x**4 - 8 * x**2 + 1)
        - 0.001 * (16 * x**5 - 20 * x**3 + 5 * x)
        + 0.0003 * (32 * x**6 - 48 * x**4 + 18 * x**2 - 1)
    )
    er = 0.1 + 0.02 * radius - 0.003 * radius**2 + 0.004 * radius**3
    gamma = np.exp(-3 + 0.1 * radius + 0.02 * radius**2 - 0.03 * radius**3)
    for actual, expected in [
        (analytic.v0(radius), expected_v0),
        (analytic.er(radius), er),
        (analytic.full_gamma(radius), gamma),
        (analytic.incident_pole(radius), er - 0.5j * gamma),
        (analytic.vres_hybrid(radius), expected_v0 + er),
    ]:
        assert actual.shape == radius.shape
        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-15)
    assert np.all(analytic.full_gamma(radius) > 0)
    assert np.isscalar(analytic.v0(2.2))
    assert np.isscalar(analytic.incident_pole(2.2))
    assert analytic.er([]).shape == (0,)


@pytest.mark.parametrize("method", ["v0", "er", "full_gamma", "incident_pole", "vres_hybrid"])
@pytest.mark.parametrize(
    "radius",
    [
        1.8,
        2.6,
        np.nextafter(1.9, -np.inf),
        np.nextafter(2.5, np.inf),
        np.nan,
        np.inf,
        -np.inf,
        2.2 + 0j,
        [2.2, 2.6],
        "2.2",
        True,
    ],
)
def test_all_entry_points_reject_invalid_radii(analytic, method, radius):
    with pytest.raises(ValueError):
        getattr(analytic, method)(radius)


def test_continuity_at_knots(analytic):
    # Independent, non-polynomial samples exercise actual piecewise coefficients.
    knots = np.linspace(1.9, 2.5, 8)
    spline = CubicSpline(knots, np.sin(8 * knots))
    curve = local.LocalRepresentation(
        [0] * 7, knots, spline.c, CubicSpline(knots, np.cos(4 * knots)).c
    )
    probe = np.r_[knots, np.linspace(1.9, 2.5, 83)]
    np.testing.assert_allclose(curve.er(probe), spline(probe), rtol=1e-13, atol=1e-15)
    for coeff in [curve.er_coefficients, curve.log_gamma_coefficients]:
        h = np.diff(knots)[:-1]
        left_value = ((coeff[0, :-1] * h + coeff[1, :-1]) * h + coeff[2, :-1]) * h
        left_value += coeff[3, :-1]
        left_d1 = (3 * coeff[0, :-1] * h + 2 * coeff[1, :-1]) * h + coeff[2, :-1]
        left_d2 = 6 * coeff[0, :-1] * h + 2 * coeff[1, :-1]
        for left, right in [
            (left_value, coeff[3, 1:]),
            (left_d1, coeff[2, 1:]),
            (left_d2, 2 * coeff[1, 1:]),
        ]:
            np.testing.assert_allclose(left, right, rtol=1e-12, atol=1e-12)


def test_constructor_guards_and_copies(analytic):
    source = np.zeros((4, 4))
    curve = local.LocalRepresentation([0] * 7, analytic.knots_bohr, source, source)
    source[:] = 9
    np.testing.assert_allclose(curve.er(2.2), 0, rtol=0, atol=0)
    with pytest.raises(ValueError):
        curve.er_coefficients[0, 0] = 3
    for coeff, knots, er in [
        ([0] * 6, analytic.knots_bohr, source),
        ([0] * 7, [1.8, 2, 2.2, 2.4, 2.5], source),
        ([0] * 7, [1.9, 2.1, 2.0, 2.3, 2.5], source),
        ([0] * 7, analytic.knots_bohr, np.zeros((3, 4))),
        ([np.nan] * 7, analytic.knots_bohr, source),
    ]:
        with pytest.raises(ValueError):
            local.LocalRepresentation(coeff, knots, er, source)


def test_width_underflow_and_overflow_are_not_physical_zero(analytic):
    for log_width in [-1000, 1000]:
        coeff = np.zeros_like(analytic.er_coefficients)
        coeff[-1] = log_width
        curve = local.LocalRepresentation(
            [0] * 7, analytic.knots_bohr, analytic.er_coefficients, coeff
        )
        with pytest.raises(ValueError, match="width"):
            curve.full_gamma(2.2)


def test_capsule_reader_authenticates_before_parsing(tmp_path):
    path = tmp_path / "not-a-tar.gz"
    path.write_bytes(b"untrusted bytes")
    with pytest.raises(ValueError, match="SHA256"):
        local.load_representation(path)


def test_reader_ignores_archive_paths_and_never_extracts(tmp_path):
    path = tmp_path / "synthetic.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, value in {"../escape.json": {"ignored": True}, "safe.json": {"x": 1}}.items():
            data = json.dumps(value).encode()
            entry = tarfile.TarInfo(name)
            entry.size = len(data)
            archive.addfile(entry, io.BytesIO(data))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert local._read_json_members(path, digest, ["safe.json"]) == {"safe.json": {"x": 1}}
    assert not (tmp_path.parent / "escape.json").exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["synthetic.tar.gz"]


@pytest.mark.parametrize("kind", ["missing", "duplicate", "symlink"])
def test_reader_rejects_ambiguous_or_link_members(tmp_path, kind):
    path = tmp_path / "synthetic.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for _ in range(2 if kind == "duplicate" else 1):
            entry = tarfile.TarInfo("other.json" if kind == "missing" else "wanted.json")
            if kind == "symlink":
                entry.type = tarfile.SYMTYPE
                entry.linkname = "../unread.json"
                archive.addfile(entry)
            else:
                entry.size = 2
                archive.addfile(entry, io.BytesIO(b"{}"))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="exactly one regular"):
        local._read_json_members(path, digest, ["wanted.json"])


def test_scoped_acceptance_is_separate_from_historical_candidate(analytic):
    payloads = {
        "representation.json": {
            "accepted_for_scoped_local_use": True,
            "no_extrapolation": True,
            "domain_bohr": [1.9, 2.5],
        },
        "inputs/frozen-c2-model.json": {
            "accepted_representation": False,
            "family": "CubicSpline_not_a_knot_Er_logFullGamma",
            "development_R_bohr": analytic.knots_bohr.tolist(),
            "models": {
                name: {
                    "knots_bohr": analytic.knots_bohr.tolist(),
                    "coefficients_descending_local_power": coefficients.tolist(),
                }
                for name, coefficients in [
                    ("Er", analytic.er_coefficients),
                    ("log_full_Gamma", analytic.log_gamma_coefficients),
                ]
            },
        },
        "inputs/accepted-neutral-result.json": {
            "status": "pass",
            "reference_R_bohr": 2.1323,
            "fits": {"5z": {"coefficients_hartree": analytic.neutral_coefficients.tolist()}},
        },
    }
    np.testing.assert_allclose(
        local._from_payloads(payloads).er(2.2), analytic.er(2.2), rtol=0, atol=1e-15
    )
    payloads["representation.json"]["accepted_for_scoped_local_use"] = False
    with pytest.raises(ValueError, match="scoped acceptance"):
        local._from_payloads(payloads)
    payloads["representation.json"]["accepted_for_scoped_local_use"] = True
    payloads["inputs/frozen-c2-model.json"]["accepted_representation"] = True
    with pytest.raises(ValueError, match="historical candidate"):
        local._from_payloads(payloads)


def test_plot_roles_and_no_overwrite(analytic, tmp_path):
    pytest.importorskip("matplotlib")
    from projects.ukrmol_co.local_representation_plot import plot_saved_curves

    pairs = [
        {
            "R_bohr": r,
            "Er_ev": hartree_to_ev(analytic.er(r)),
            "full_Gamma_ev": hartree_to_ev(analytic.full_gamma(r)),
            "role": role,
        }
        for role, radii in [
            ("C2_development", np.linspace(1.9, 2.5, 8)),
            ("C2_prospective_validation", [2.02, 2.27]),
        ]
        for r in radii
    ]
    evidence = {
        "ten-qualified-pairs.json": {"pairs": pairs},
        "inputs/accepted-neutral-result.json": {
            "fits": {
                "5z": {
                    "data": [
                        {
                            "R_bohr": 2.2,
                            "actual_hartree": analytic.v0(2.2),
                            "error_mev": 0,
                            "historical_role": "historical_knot",
                        }
                    ]
                }
            }
        },
    }
    for label, r in [("lower", 2.02), ("upper", 2.27)]:
        evidence[f"inputs/{label}-frozen-score.json"] = {
            "geometry_bohr": r,
            "sectors": {
                s: {
                    "position_ev": hartree_to_ev(analytic.er(r)),
                    "full_width_ev": hartree_to_ev(analytic.full_gamma(r)),
                }
                for s in ["B1", "B2"]
            },
        }
    paths = plot_saved_curves(analytic, evidence, tmp_path / "plots")
    assert len(paths) == 4 and all(p.stat().st_size > 0 for p in paths)
    svg = (tmp_path / "plots/local-curves.svg").read_text()
    for text in [
        "8 development nodes",
        "2 prospective validation points",
        "Mixed-method",
        "FULL width",
        "NOT a correlated anion total",
    ]:
        assert text in svg
    with pytest.raises(FileExistsError):
        plot_saved_curves(analytic, evidence, tmp_path / "plots")
    with pytest.raises(ValueError, match="outside"):
        plot_saved_curves(analytic, evidence, Path(__file__).parent / "unwritten-plot-test")


def test_saved_delivery_regression():
    path = os.environ.get("CO_LOCAL_REPRESENTATION_CAPSULE")
    if not path:
        pytest.skip("Set CO_LOCAL_REPRESENTATION_CAPSULE to the pinned capsule; no auto-fetch")
    curve = local.load_representation(path)
    evidence = local.load_plot_evidence(path)
    model = evidence["inputs/frozen-c2-model.json"]
    neutral = evidence["inputs/accepted-neutral-result.json"]
    assert model["accepted_representation"] is False
    rows = neutral["fits"]["5z"]["data"]
    np.testing.assert_allclose(
        curve.v0([r["R_bohr"] for r in rows]),
        [r["fitted_hartree"] for r in rows],
        rtol=0,
        atol=1e-12,
    )
    np.testing.assert_allclose(curve.v0(2.1323), 0, rtol=0, atol=1e-12)
    pairs = evidence["ten-qualified-pairs.json"]["pairs"]
    development = [p for p in pairs if p["role"] == "C2_development"]
    validation = [p for p in pairs if p["role"] == "C2_prospective_validation"]
    assert len(development) == 8 and len(validation) == 2
    assert [r["R_bohr"] for r in development] == model["development_R_bohr"]
    np.testing.assert_allclose(
        curve.incident_pole([p["R_bohr"] for p in development]),
        [complex(*p["incident_z_hartree"]) for p in development],
        rtol=0,
        atol=1e-13,
    )
    for prediction in model["predictions"]:
        radius = prediction["validation_R_bohr"]
        er = hartree_to_ev(curve.er(radius))
        gamma = hartree_to_ev(curve.full_gamma(radius))
        np.testing.assert_allclose(
            [er, gamma], [prediction["Er_ev"], prediction["full_Gamma_ev"]], rtol=0, atol=1e-12
        )
        label = "lower" if radius == min(p["R_bohr"] for p in validation) else "upper"
        score = evidence[f"inputs/{label}-frozen-score.json"]
        assert score["status"] == "pass" and score["refits"] == score["coefficient_changes"] == 0
        for row in score["sectors"].values():
            error = [row["position_ev"] - er, row["full_width_ev"] - gamma]
            np.testing.assert_allclose(
                error,
                [row["signed_error_ev"]["position_ev"], row["signed_error_ev"]["full_width_ev"]],
                rtol=0,
                atol=1e-12,
            )
            assert abs(error[0]) <= 0.005 and abs(error[1]) <= max(0.01 * gamma, 0.001)
