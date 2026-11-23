"""Selected-pair provenance and the midpoint derivative, not molecular fitting."""

import cmath
from types import SimpleNamespace

import numpy as np
import pytest
from qscat.model import N2

from projects.potential_factory import tracker


def diagonal_pair():
    """The selected pair is deliberately not nearest either matrix's midpoint."""
    ra = np.array([1.0, 2.0], complex)
    rb = np.array([1.0 + 0.3j, 2.0 + 0.3j])
    reference = tracker.well_potential(1, 2.4, 0.48, None)
    ta = np.diag([-0.2, 0.03]) - np.diag(reference(ra))
    tb = np.diag([0.2, 0.04]) - np.diag(reference(rb))
    return SimpleNamespace(
        grid_a=SimpleNamespace(points=ra),
        grid_b=SimpleNamespace(points=rb),
        hamiltonians=lambda v: (ta + np.diag(v(ra)), tb + np.diag(v(rb))),
    )


def test_midpoint_gradient_uses_selected_eigenvalues_and_each_complex_grid():
    pair = diagonal_pair()
    pole = tracker.Pole(0j, 0.4, angle_energies=(-0.2 + 0j, 0.2 + 0j))
    actual = tracker.well_pole_gradient(pair, 1, tracker.WellParams(2.4, 0.48), pole)
    ra, rb = 1.0, 1.0 + 0.3j
    ga, gb = cmath.exp(-0.48 * ra**2), cmath.exp(-0.48 * rb**2)
    expected = (-(ga + gb) / 2, 2.4 * (ra**2 * ga + rb**2 * gb) / 2)
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)


@pytest.mark.parametrize("case", ["missing", "mean", "residual", "stale"])
def test_midpoint_gradient_rejects_missing_inconsistent_or_stale_provenance(case):
    pair = diagonal_pair()
    if case == "missing":
        pole = tracker.Pole(0j, 0.4)
    elif case == "mean":
        pole = tracker.Pole(0.1 + 0j, 0.4, angle_energies=(-0.2 + 0j, 0.2 + 0j))
    elif case == "residual":
        pole = tracker.Pole(0j, 0.3, angle_energies=(-0.2 + 0j, 0.2 + 0j))
    else:
        pole = tracker.Pole(0.5 + 0j, 0.4, angle_energies=(0.3 + 0j, 0.7 + 0j))
    with pytest.raises(ValueError, match=r"matched|provenance|eigenvalue"):
        tracker.well_pole_gradient(pair, 1, tracker.WellParams(2.4, 0.48), pole)


def test_pole_retains_the_literal_selected_angle_pair(monkeypatch):
    spectra = iter(
        [np.array([0.04 - 0.01j, 0.08 - 0.01j]), np.array([0.0402 - 0.01j, 0.0801 - 0.01j])]
    )
    pair = object.__new__(tracker.ElectronicPair)
    monkeypatch.setattr(pair, "hamiltonians", lambda _: (np.eye(2), np.eye(2)))
    monkeypatch.setattr(tracker, "eigen", lambda _: (next(spectra), np.eye(2)))
    pole = pair.pole(lambda r: r, tracker.DEFAULT_RES_WINDOW)
    assert pole is not None
    assert pole.angle_energies == (0.08 - 0.01j, 0.0801 - 0.01j)
    assert sum(pole.angle_energies) / 2 == pytest.approx(pole.energy, rel=0, abs=1e-14)


def test_joint_polish_uses_the_shared_midpoint_gradient(monkeypatch):
    from projects.potential_factory import fit as fitting
    from projects.potential_factory.ansatz import from_diatomic

    model = from_diatomic(N2)
    radii = np.array([2.0, 2.2])
    returned = tracker.Pole(0.08 - 0.01j, 0.0)
    pair = SimpleNamespace(pole=lambda *_: returned)
    calls = []
    captured = {}
    derivatives = (2.0 + 0.3j, 4.0 + 0.2j)

    def gradient(*args, **kwargs):
        calls.append(args)
        return derivatives

    def least_squares(residual, x0, *, jac, **kwargs):
        captured["jacobian"] = jac(x0)
        assert residual(x0).shape == (4,)
        return SimpleNamespace(x=x0, status=1, nfev=1)

    monkeypatch.setattr(fitting, "well_pole_gradient", gradient, raising=False)
    monkeypatch.setattr(fitting, "least_squares", least_squares)
    fitted, _ = fitting._joint_polish(
        model,
        pair,
        radii,
        lambda _: returned.energy,
        np.full(2, 0.02),
        np.ones(2, bool),
        [None, None],
        None,
    )
    assert len(calls) == 2
    complex_rows = np.array(
        [
            np.r_[
                derivatives[0] * fitting._smooth_grad(fitted.lam, r),
                derivatives[1] * fitting._smooth_grad(fitted.alpha, r),
            ]
            for r in radii
        ]
    )
    np.testing.assert_allclose(
        captured["jacobian"],
        np.vstack([complex_rows.real, complex_rows.imag]),
        rtol=1e-13,
        atol=1e-13,
    )
