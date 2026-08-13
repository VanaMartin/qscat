"""qscat.viz: domain colouring + the equidistant projector."""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.grids import electronic_grid, nuclear_grid
from qscat.dvr import TensorGrid
from qscat.viz import EquidistantProjector, complex_to_hsv, complex_to_rgb

# --- coloring ---------------------------------------------------------------


def test_complex_to_hsv_phase_is_hue() -> None:
    # Positive real -> hue 0; +i -> 0.25; -1 -> 0.5; -i -> 0.75.
    z = np.array([1.0, 1.0j, -1.0, -1.0j], dtype=complex)
    hsv = complex_to_hsv(z, mag=1.0)
    assert np.allclose(hsv[:, 0], [0.0, 0.25, 0.5, 0.75])


def test_complex_to_hsv_magnitude_maps_to_value() -> None:
    # Below mag: value = |z|/mag, full saturation. Above mag: value clips to 1,
    # saturation drops (white-wash).
    z = np.array([0.5, 2.0], dtype=complex)
    hsv = complex_to_hsv(z, mag=1.0)
    assert np.isclose(hsv[0, 2], 0.5) and np.isclose(hsv[0, 1], 1.0)  # |z|=0.5
    assert np.isclose(hsv[1, 2], 1.0) and np.isclose(hsv[1, 1], 0.5)  # |z|=2


def test_complex_to_rgb_zero_is_black_and_large_is_white() -> None:
    rgb = complex_to_rgb(np.array([0.0 + 0j, 1e6 + 0j]), mag=1.0)
    assert np.allclose(rgb[0], [0.0, 0.0, 0.0])  # zero -> black
    assert np.allclose(rgb[1], [1.0, 1.0, 1.0], atol=1e-3)  # huge -> white


def test_complex_to_hsv_rejects_real_input() -> None:
    with pytest.raises(ValueError, match="complex"):
        complex_to_hsv(np.array([1.0, 2.0]))


# --- projector --------------------------------------------------------------


def _tgrid() -> TensorGrid:
    return TensorGrid([
        electronic_grid(r_max=16.0, order=7, n_complex=4),
        nuclear_grid(r_max=22.0, quadrature=8, n_complex=4),
    ])


def test_projector_reproduces_separable_polynomial_exactly() -> None:
    # A separable low-degree polynomial (degree <= nq-1 per axis) is reproduced
    # to round-off by the DVR interpolation in the REAL region (z(x)=x there):
    # field(r,R) = gr(r)*hR(R) exactly. Build state_ij = sqrt(w_i w_j) f(r_i,R_j).
    tg = _tgrid()
    g0, g1 = tg.grids
    m0 = g0.real_points <= g0.R0
    m1 = g1.real_points <= g1.R0

    def gr(r: np.ndarray) -> np.ndarray:
        return 1.0 + 0.1 * r - 0.02 * r**2  # degree 2 <= 6 (electronic order 7)

    def hR(R: np.ndarray) -> np.ndarray:
        return 2.0 - 0.3 * R + 0.05 * R**2  # degree 2 <= 7 (nuclear order 8)

    vr = np.where(m0, np.sqrt(g0.weights) * gr(g0.real_points), 0.0)
    vR = np.where(m1, np.sqrt(g1.weights) * hR(g1.real_points), 0.0)
    M = np.outer(vr, vR)

    proj = EquidistantProjector(tg, samples=(40, 40), extent=((1.0, 6.0), (1.0, 4.0)))
    field = proj.project(M.reshape(-1))
    exact = np.outer(gr(proj.axis0), hR(proj.axis1))
    assert field.shape == (40, 40)
    assert np.max(np.abs(field.real - exact)) < 1e-10


def test_projector_exact_at_nodes() -> None:
    # Sampling at real node coordinates recovers state_ij / sqrt(w_i w_j) exactly
    # (Lagrange cardinal property), a differential check on the projector.
    tg = _tgrid()
    g0, g1 = tg.grids
    r_nodes = g0.real_points[g0.real_points <= g0.R0]
    R_nodes = g1.real_points[g1.real_points <= g1.R0]
    proj = EquidistantProjector(
        tg,
        samples=(r_nodes.size, R_nodes.size),
        extent=((r_nodes[0], r_nodes[-1]), (R_nodes[0], R_nodes[-1])),
    )
    # linspace over [first,last] node won't hit interior nodes; instead build the
    # projector's 1-D operators directly at the exact node coords via the public
    # dvr_interpolation_matrix and check the cardinal property.
    from qscat.dvr import dvr_interpolation_matrix

    P0 = dvr_interpolation_matrix(g0, r_nodes)
    rng = np.random.default_rng(1)
    state = rng.standard_normal(g0.n) + 1j * rng.standard_normal(g0.n)
    got = P0 @ state
    idx = np.flatnonzero(g0.real_points <= g0.R0)
    assert np.allclose(got, (state / np.sqrt(g0.weights))[idx], atol=1e-12)
    assert proj.axis0.size == r_nodes.size  # sanity on construction


def test_projector_rejects_non_2d_grid() -> None:
    from qscat.exceptions import GridError

    tg1 = TensorGrid([electronic_grid(r_max=16.0, order=6, n_complex=3)])
    with pytest.raises(GridError, match="2-D"):
        EquidistantProjector(tg1, samples=10)


def test_plot_wavefunction_2d_writes_png(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    from qscat.viz import plot_wavefunction_2d

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(30, 30))
    rng = np.random.default_rng(0)
    state = rng.standard_normal(tg.grids[0].n * tg.grids[1].n) + 0j
    out = tmp_path / "psi.png"
    plot_wavefunction_2d(proj, state, mag=0.05, path=out, ylabel="r", xlabel="R")
    assert out.exists() and out.stat().st_size > 0


def test_projector_project_values_interpolates_nodal_field() -> None:
    # project_values interpolates a nodal-VALUE field (no 1/sqrt(w)): a separable
    # low-degree polynomial is reproduced to round-off in the real region.
    tg = _tgrid()
    g0, g1 = tg.grids
    m0 = g0.real_points <= g0.R0
    m1 = g1.real_points <= g1.R0

    def vr(r: np.ndarray) -> np.ndarray:
        return 0.5 - 0.1 * r + 0.02 * r**2

    def vR(R: np.ndarray) -> np.ndarray:
        return -1.0 + 0.3 * R

    V = np.outer(np.where(m0, vr(g0.real_points), 0.0), np.where(m1, vR(g1.real_points), 0.0))
    proj = EquidistantProjector(tg, samples=(30, 30), extent=((1.0, 6.0), (1.0, 4.0)))
    got = np.real(proj.project_values(V))
    exact = np.outer(vr(proj.axis0), vR(proj.axis1))
    assert np.max(np.abs(got - exact)) < 1e-10


def test_contours_magnitude_thin_white_06(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.contour import QuadContourSet
    from qscat.viz import plot_wavefunction_2d

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(40, 40))
    rng = np.random.default_rng(0)
    state = rng.standard_normal(tg.grids[0].n * tg.grids[1].n) + 0j
    _, ax = plt.subplots()
    plot_wavefunction_2d(proj, state, mag=0.05, ax=ax, contours=True)
    csets = [c for c in ax.get_children() if isinstance(c, QuadContourSet)]
    assert len(csets) == 1
    cset = csets[0]
    # Thin white lines at 0.6 opacity, and magnitude levels derived from mag.
    assert np.allclose(cset.get_alpha(), 0.6)
    assert np.allclose(cset.get_linewidths(), 0.6)
    assert np.allclose(cset.get_edgecolor()[0][:3], [1.0, 1.0, 1.0])  # white
    # magnitude levels are k*mag/5 -> evenly spaced
    assert np.allclose(np.diff(cset.levels), cset.levels[1] - cset.levels[0])
    plt.close("all")


def test_contours_potential_array_and_callable_both_draw() -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.contour import QuadContourSet
    from qscat.viz import plot_wavefunction_2d

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(30, 30), extent=((1.0, 6.0), (1.0, 4.0)))
    state = np.zeros(tg.grids[0].n * tg.grids[1].n, dtype=complex)

    # Nodal-array potential on the same grid (the numerically-evaluated path).
    V = np.outer(tg.grids[0].real_points, tg.grids[1].real_points).astype(complex)
    _, ax1 = plt.subplots()
    plot_wavefunction_2d(
        proj, state, mag=0.05, ax=ax1, contours=6,
        contour_field="potential", potential=V, contour_color="grey",
    )
    assert len([c for c in ax1.get_children() if isinstance(c, QuadContourSet)]) == 1

    # Callable potential (analytic path): V(r, R) on the sampling meshgrid.
    _, ax2 = plt.subplots()
    plot_wavefunction_2d(
        proj, state, mag=0.05, ax=ax2, contours=6,
        contour_field="potential", potential=lambda r, R: r + R,
    )
    assert len([c for c in ax2.get_children() if isinstance(c, QuadContourSet)]) == 1
    plt.close("all")


def test_contours_potential_requires_potential() -> None:
    pytest.importorskip("matplotlib")
    from qscat.viz import plot_wavefunction_2d

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(20, 20))
    state = np.zeros(tg.grids[0].n * tg.grids[1].n, dtype=complex)
    with pytest.raises(ValueError, match="requires potential"):
        plot_wavefunction_2d(proj, state, mag=0.05, contours=True, contour_field="potential")


def test_contours_false_draws_nothing() -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.contour import QuadContourSet
    from qscat.viz import plot_wavefunction_2d

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(20, 20))
    state = np.zeros(tg.grids[0].n * tg.grids[1].n, dtype=complex)
    _, ax = plt.subplots()
    plot_wavefunction_2d(proj, state, mag=0.05, ax=ax)  # contours default False
    assert not [c for c in ax.get_children() if isinstance(c, QuadContourSet)]
    plt.close("all")


# --- combined potential + wavefunction overlay ------------------------------


def test_energy_contour_levels_selection() -> None:
    from qscat.viz import energy_contour_levels

    eps = np.array([0.01, 0.03, 0.06, 0.10])
    # thresholds only
    lv = energy_contour_levels(eps=eps)
    assert lv == [0.01, 0.03, 0.06, 0.10]
    # thresholds + total energies eps[0] + E
    lv = energy_contour_levels(eps=eps, v_init=0, energies=[0.05, 0.09])
    assert 0.06 in lv and pytest.approx(0.10, abs=1e-9) in lv  # 0.01+0.05, 0.01+0.09
    # e_range clip
    lv = energy_contour_levels(eps=eps, e_range=(0.02, 0.07))
    assert lv == [0.03, 0.06]
    # min_spacing thinning + sorted
    lv = energy_contour_levels(eps=np.array([0.0, 0.005, 0.02, 0.021, 0.05]), min_spacing=0.01)
    assert lv == [0.0, 0.02, 0.05]
    # max_levels cap
    lv = energy_contour_levels(eps=np.linspace(0.0, 1.0, 50), max_levels=6)
    assert len(lv) <= 6 and lv == sorted(lv)
    # energies without eps raises
    with pytest.raises(ValueError, match="requires eps"):
        energy_contour_levels(energies=[0.1])


def test_combined_overlays_dotted_potential_plus_psi() -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.contour import QuadContourSet
    from qscat.viz import plot_wavefunction_2d

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(40, 40), extent=((1.0, 6.0), (1.0, 4.0)))
    rng = np.random.default_rng(0)
    state = rng.standard_normal(tg.grids[0].n * tg.grids[1].n) + 0j
    # A smooth "potential" on the same grid so auto levels land inside its range.
    V = np.add.outer(tg.grids[0].real_points, tg.grids[1].real_points).astype(complex)

    _, ax = plt.subplots()
    plot_wavefunction_2d(
        proj, state, mag=0.05, ax=ax,
        contours=True,  # solid white |psi|
        potential=V, potential_levels="auto",  # dotted potential overlay
        eps=np.array([2.5, 4.0, 6.0]), v_init=0, energies=[1.0, 2.0],
    )
    csets = [c for c in ax.get_children() if isinstance(c, QuadContourSet)]
    assert len(csets) == 2  # |psi| + potential
    # Inline energy labels were drawn by clabel (Text children on the axes).
    assert any(t.get_text() for t in ax.texts)
    plt.close("all")


def test_combined_potential_only_overlay_independent_of_contours() -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.contour import QuadContourSet
    from qscat.viz import plot_wavefunction_2d

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(30, 30), extent=((1.0, 6.0), (1.0, 4.0)))
    state = np.zeros(tg.grids[0].n * tg.grids[1].n, dtype=complex)
    V = np.add.outer(tg.grids[0].real_points, tg.grids[1].real_points).astype(complex)
    _, ax = plt.subplots()
    # contours=False, but potential overlay still draws (independent).
    plot_wavefunction_2d(
        proj, state, mag=0.05, ax=ax,
        potential=V, potential_levels=[4.0, 6.0], potential_labels=False,
    )
    assert len([c for c in ax.get_children() if isinstance(c, QuadContourSet)]) == 1
    plt.close("all")


# --- artist + animation -----------------------------------------------------


def test_wavefunction_artist_update_changes_image_keeps_static_potential() -> None:
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.contour import QuadContourSet
    from qscat.viz import WavefunctionArtist

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(30, 30), extent=((1.0, 6.0), (1.0, 4.0)))
    V = np.add.outer(tg.grids[0].real_points, tg.grids[1].real_points).astype(complex)
    n = tg.grids[0].n * tg.grids[1].n
    rng = np.random.default_rng(0)

    _, ax = plt.subplots()
    artist = WavefunctionArtist(
        ax, proj, mag=0.5, contours=True,
        potential=V, potential_levels=[4.0, 6.0], potential_labels=False,
    )
    # static potential overlay drawn at construction (before any state)
    assert len([c for c in ax.get_children() if isinstance(c, QuadContourSet)]) == 1

    img1 = artist.update(rng.standard_normal(n) + 0j)[0].get_array().copy()
    img2 = artist.update(2.0 * (rng.standard_normal(n) + 0j))[0].get_array()
    assert not np.allclose(img1, img2)  # image updated for the new state
    # one static potential set + exactly one (refreshed) |psi| set -> 2 total
    assert len([c for c in ax.get_children() if isinstance(c, QuadContourSet)]) == 2
    plt.close("all")


def test_animate_wavefunction_writes_gif(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    pytest.importorskip("PIL")  # PillowWriter
    from qscat.viz import animate_wavefunction

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(24, 24))
    n = tg.grids[0].n * tg.grids[1].n
    rng = np.random.default_rng(0)
    frames = [np.exp(1j * k) * (rng.standard_normal(n) + 0j) for k in range(3)]
    out = tmp_path / "psi.gif"
    anim = animate_wavefunction(
        proj, frames, mag=0.5, times=[0.0, 1.0, 2.0], outfile=out, fps=5, contours=True
    )
    assert out.exists() and out.stat().st_size > 0
    assert anim is not None


def test_animate_wavefunction_empty_frames_raises() -> None:
    pytest.importorskip("matplotlib")
    from qscat.viz import animate_wavefunction

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(10, 10))
    with pytest.raises(ValueError, match="empty"):
        animate_wavefunction(proj, [], mag=0.5)


def test_animate_artists_multi_panel(tmp_path) -> None:
    pytest.importorskip("matplotlib")
    pytest.importorskip("PIL")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from qscat.viz import WavefunctionArtist, animate_artists

    tg = _tgrid()
    proj = EquidistantProjector(tg, samples=(20, 20))
    n = tg.grids[0].n * tg.grids[1].n
    rng = np.random.default_rng(1)
    frames = [rng.standard_normal(n) + 0j for _ in range(2)]

    fig, (axl, axr) = plt.subplots(1, 2)
    a_left = WavefunctionArtist(axl, proj, mag=0.5)
    a_right = WavefunctionArtist(axr, proj, mag=0.5, contours=True)
    out = tmp_path / "two.gif"
    animate_artists(fig, [(a_left, frames), (a_right, frames)], outfile=out, fps=5)
    assert out.exists() and out.stat().st_size > 0
    plt.close("all")
