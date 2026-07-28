from __future__ import annotations

from qscat.core.grids import electronic_grid, nuclear_grid, segmented_grid
from qscat.model import N2
from qscat.tuning import (
    probe_channel_representation,
    probe_electronic,
    probe_nuclear,
    refine,
)

# The N2 nuclear real region always sums to 12.0 bohr (see
# `qscat.core.grids.nuclear_grid`'s `_REAL_SEGMENTS`), independent of
# `r_max`/`n_complex` -- a reasonable stand-in R for the dissociation limit
# used by `probe_electronic` below.
_R_DISSOC = 12.0


def test_refine_doubles_real_elements_keeps_tail_and_span():
    g = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)
    r = refine(g)
    assert len(r.spec.elements) == len(g.spec.elements) + 23  # 23 real elements doubled
    # real span unchanged
    real_len_g = sum(e.length for e in g.spec.elements if e.angle_deg == 0.0)
    real_len_r = sum(e.length for e in r.spec.elements if e.angle_deg == 0.0)
    assert abs(real_len_g - real_len_r) < 1e-12
    # ECS tail untouched (same lengths/angles, same count)
    tail_g = [(e.length, e.angle_deg) for e in g.spec.elements if e.angle_deg != 0.0]
    tail_r = [(e.length, e.angle_deg) for e in r.spec.elements if e.angle_deg != 0.0]
    assert tail_g == tail_r
    assert r.spec.quadrature == g.spec.quadrature
    assert r.n > g.n


def test_probe_nuclear_converged_on_known_good_n2_grid():
    g = nuclear_grid(r_max=22.0, n_complex=8, quadrature=12)
    result = probe_nuclear(N2, g, 3, rtol=1e-3)
    assert result.converged, result.detail
    assert result.cost == g.n
    assert len(result.value) == 3


def test_probe_electronic_converged_on_known_good_grid():
    g = electronic_grid(r_max=16.0, order=7, n_complex=6)
    result = probe_electronic(N2, g, _R_DISSOC, window=None, rtol=1e-3)
    assert result.converged, result.detail
    assert result.cost == g.n


def test_probe_channel_representation_converged_modest_k():
    # k=1.0 (wavelength ~6.3 bohr) is well within what electronic_grid's
    # default 0.2-4 bohr real-region elements resolve at order 7.
    g = electronic_grid(r_max=16.0, order=7, n_complex=6)
    result = probe_channel_representation(g, k=1.0, l=0, rtol=1e-3)
    assert result.converged, result.detail


def test_probe_channel_representation_fails_on_coarse_grid_fast_wave():
    # A deliberately coarse grid: uniform 1.0-bohr real elements out to 10
    # bohr, order-10 quadrature -- and a K~58 outgoing wave (the F2 DA wave
    # from `test_tuning_ecs.test_tail_absorbs_fast_wave`). Wavelength
    # 2*pi/58 ~ 0.108 bohr is ~9x smaller than the element -- unresolvable.
    coarse = segmented_grid(
        real_segments=[(10, 10.0)],
        complex_segments=[(4, 20.0)],
        angle_deg=35.0,
        quadrature=10,
    )
    result = probe_channel_representation(coarse, k=58.0, l=0, rtol=1e-3)
    assert not result.converged
    assert result.detail["rel_error"] > 1e-2  # >20x rtol: badly wrong, not a marginal miss


def test_probe_channel_representation_honors_nonunit_mass():
    # charge=0 with a non-electron mass (e.g. a nuclear dissociation wave) --
    # riccati_bessel_en_mass's r-dependence is mass-independent (only the
    # overall normalization scales), so this is a smoke test that the mass
    # path runs and still converges on a grid well-resolved for k=1.0.
    g = electronic_grid(r_max=16.0, order=7, n_complex=6)
    result = probe_channel_representation(g, k=1.0, l=0, mass=918.25, rtol=1e-3)
    assert result.converged, result.detail


def test_probe_channel_representation_fine_grid_resolves_same_fast_wave():
    # Same K~58 wave, but on a grid whose elements (0.05 bohr, roughly half
    # the K~58 wavelength of 2*pi/58~0.108 bohr) are sized for it -- the
    # probe must report converged=True, showing this is a genuine
    # resolution effect and not an artifact of the K=58 value itself.
    fine = segmented_grid(
        real_segments=[(200, 10.0)],
        complex_segments=[(4, 20.0)],
        angle_deg=35.0,
        quadrature=10,
    )
    result = probe_channel_representation(fine, k=58.0, l=0, rtol=1e-3)
    assert result.converged, result.detail
