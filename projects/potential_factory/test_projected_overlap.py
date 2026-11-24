"""Analytic piecewise-polynomial controls for different-grid real-region identity."""

import numpy as np
import pytest
from numpy.polynomial import Polynomial
from qscat.dvr import ElementSpec, FemDvrEcsGrid, GridSpec

from projects.potential_factory import projected_overlap
from projects.potential_factory.projected_overlap import projected_real_region_gram


def _grids():
    ra = FemDvrEcsGrid(GridSpec(6, [ElementSpec(1)] * 3 + [ElementSpec(1, 30)]))
    rb = FemDvrEcsGrid(GridSpec(8, [ElementSpec(0.5)] * 6 + [ElementSpec(2, 35)]))
    na = FemDvrEcsGrid(GridSpec(6, [ElementSpec(0.1)] * 6, x_min=1.9))
    nb = FemDvrEcsGrid(GridSpec(8, [ElementSpec(0.14)] * 4, x_min=1.92))
    return ra, rb, na, nb


def _fields(odd=False):
    ra, rb, na, nb = _grids()
    enda = 3 + np.exp(1j * np.pi / 6)
    endb = 3 + 2 * np.exp(1j * np.deg2rad(35))
    radial = Polynomial.fromroots([0, enda, endb])
    nuclear = Polynomial.fromroots([-0.3, -0.28, 0.28, 0.3])
    fields = []
    for i, (r, n) in enumerate([(ra, na), (rb, nb)]):
        f = nuclear(n.points.real - 2.2)
        if odd and i == 1:
            f *= n.points.real - 2.2
        b = f * np.sqrt(n.weights.real)
        norm = np.linalg.norm(b)
        B = (b / norm)[:, None]
        phase = 1 if i == 0 else 1.3 * np.exp(0.7j)
        V = (phase * norm * radial(r.points) * np.sqrt(r.weights))[:, None]
        fields.append((r, n, B, V))
    return fields, radial, nuclear


def test_same_polynomial_has_unit_overlap_across_boxes_and_ecs_tails():
    fields, radial, nuclear = _fields()
    G = projected_real_region_gram(*fields, radial_stop=2.5)
    r_integral = (radial * Polynomial(radial.coef.conjugate())).integ()
    R_integral = (nuclear * nuclear).integ()
    norm = (r_integral(2.5) - r_integral(0)).real * (R_integral(0.28) - R_integral(-0.28))
    np.testing.assert_allclose(G.diagonal().real, [norm, 1.3**2 * norm], rtol=1e-10, atol=1e-14)
    np.testing.assert_allclose(G[0, 1], 1.3 * np.exp(0.7j) * norm, rtol=1e-10, atol=1e-14)
    score = abs(G[0, 1]) ** 2 / (G[0, 0].real * G[1, 1].real)
    np.testing.assert_allclose(score, 1, rtol=0, atol=1e-12)
    np.testing.assert_allclose(G, G.conj().T, rtol=0, atol=1e-14)


def test_odd_nuclear_field_is_orthogonal_on_common_symmetric_interval():
    fields, _, _ = _fields(odd=True)
    G = projected_real_region_gram(*fields, radial_stop=2.5)
    assert abs(G[0, 1]) ** 2 / (G[0, 0].real * G[1, 1].real) < 1e-18


def test_streamed_and_unchunked_quadrature_agree(monkeypatch):
    fields, _, _ = _fields()
    expected = projected_real_region_gram(*fields, radial_stop=2.5)
    monkeypatch.setattr(projected_overlap, "_FIELD_POINT_BUDGET", 300)
    actual = projected_real_region_gram(*fields, radial_stop=2.5)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("which", ["shape", "quadrature", "complex_nuclear"])
def test_invalid_comparison_rejects(which):
    fields, _, _ = _fields()
    kwargs = dict(radial_stop=2.5)
    if which == "shape":
        r, n, B, V = fields[0]
        fields[0] = (r, n, B, V[:-1])
    elif which == "quadrature":
        kwargs["quadrature"] = 2
    else:
        r, _, B, V = fields[0]
        n = FemDvrEcsGrid(GridSpec(6, [ElementSpec(0.1, 10)] * 6, x_min=1.9))
        fields[0] = (r, n, B, V)
    with pytest.raises(ValueError):
        projected_real_region_gram(*fields, **kwargs)
