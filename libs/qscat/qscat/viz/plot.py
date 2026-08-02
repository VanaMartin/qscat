"""Render a projected 2-D wavefunction as a domain-coloured image.

matplotlib is imported lazily (the optional ``qscat[plot]`` extra), mirroring
`qscat.core.plot`. `plot_wavefunction_2d` projects a state through an
`EquidistantProjector`, domain-colours it (`complex_to_rgb`), and draws/saves the
image -- the single-frame building block for animations.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from .coloring import complex_to_rgb
from .projector import EquidistantProjector

__all__ = ["plot_wavefunction_2d"]

_PathLike = str | os.PathLike[str]


def _contour_levels(
    contours: bool | int | Sequence[float],
    contour_field: str,
    mag: float,
    z: npt.NDArray[np.float64],
) -> list[float]:
    """Resolve contour levels; magnitude levels are derived from `mag`."""
    zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
    if contours is True or isinstance(contours, int):  # bool True or a count
        n = 20 if contours is True else int(contours)
        if contour_field == "magnitude":
            levels = [k * mag / 5.0 for k in range(1, n + 1)]  # tied to the colour scale
        else:  # potential: mag is a wavefunction scale, so span the field's range
            if zmax <= zmin:
                return []
            levels = list(np.linspace(zmin, zmax, n + 2)[1:-1])
    else:  # explicit sequence
        levels = sorted(float(v) for v in contours)
    # Keep only levels that fall within the data range (avoids empty-contour warns).
    return [v for v in levels if zmin < v < zmax]


def plot_wavefunction_2d(
    projector: EquidistantProjector,
    state: npt.NDArray[np.complex128],
    *,
    mag: float,
    path: _PathLike | None = None,
    inverse: bool = False,
    title: str | None = None,
    xlabel: str = "axis 1",
    ylabel: str = "axis 0",
    ax: Any = None,
    contours: bool | int | Sequence[float] = False,
    contour_field: Literal["magnitude", "potential"] = "magnitude",
    potential: npt.NDArray[np.complex128]
    | Callable[[npt.NDArray[np.float64], npt.NDArray[np.float64]], npt.NDArray[np.float64]]
    | None = None,
    contour_color: str = "white",
    contour_alpha: float = 0.6,
    contour_linewidth: float = 0.6,
) -> Any:
    """Project, domain-colour, and draw a 2-D state; save to ``path`` if given.

    Parameters
    ----------
    projector : EquidistantProjector
        A cached projector for the state's tensor grid.
    state : ndarray
        The flat 2-D complex state to render.
    mag : float
        Magnitude mapping to full brightness (see `complex_to_hsv`).
    path : path-like, optional
        If given, the figure is saved here (PNG).
    inverse : bool, optional
        Use the light-background (inverse) colour mapping.
    title, xlabel, ylabel : str, optional
        Axis labels; ``axis 0`` is the projector's first grid, ``axis 1`` the
        second (e.g. ``ylabel="electronic r"``, ``xlabel="nuclear R"``).
    ax : matplotlib Axes, optional
        Draw into an existing Axes (for composition / animation); a new figure is
        created when None.
    contours : bool or int or sequence of float, optional
        Contour overlay: ``False`` (default) draws none; ``True`` uses default
        levels; an ``int`` requests that many levels; a sequence gives explicit
        levels. Drawn as thin lines (see the ``contour_*`` params).
    contour_field : {"magnitude", "potential"}, optional
        What the contours trace: ``"magnitude"`` = ``|psi|`` from the projected
        state (levels derived from ``mag``: ``k*mag/5``); ``"potential"`` =
        the field passed as ``potential`` (levels span its range).
    potential : ndarray or callable, optional
        Required when ``contour_field="potential"``. Either a real potential as
        nodal values on the SAME tensor grid as ``state`` (shape ``(n0_nodes,
        n1_nodes)`` or flat), projected via `EquidistantProjector.project_values`
        -- the robust path for a numerically-evaluated potential known only on
        the grid -- OR a callable ``V(r, R)`` (analytic only) taking 2-D
        meshgrids ``(axis0, axis1)`` and returning the surface directly.
    contour_color, contour_alpha, contour_linewidth : optional
        Contour style; defaults are thin white lines at 0.6 opacity.

    Returns
    -------
    matplotlib.image.AxesImage
        The drawn image handle (useful for `set_data` in an animation loop).

    Raises
    ------
    ModuleNotFoundError
        If matplotlib is not installed (`pip install 'qscat[plot]'`).
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:  # pragma: no cover - trivial guard
        raise ModuleNotFoundError(
            "qscat.viz.plot_wavefunction_2d requires matplotlib. "
            "Install the plotting extra: pip install 'qscat[plot]'."
        ) from exc

    field = projector.project(state)
    rgb = complex_to_rgb(field, mag, inverse=inverse)

    a0, a1 = projector.axis0, projector.axis1
    # rows = axis 0 (drawn vertically), cols = axis 1 (horizontal); origin upper.
    extent = [float(a1[0]), float(a1[-1]), float(a0[-1]), float(a0[0])]

    created = ax is None
    if created:
        _, ax = plt.subplots(figsize=(8, 6))
    image = ax.imshow(rgb, origin="upper", aspect="auto", extent=extent)

    if contours is not False:
        if contour_field == "magnitude":
            z = np.abs(field)  # |psi|, contoured on the same sampling grid
        else:  # potential (can be negative -> use .real, not abs)
            if potential is None:
                raise ValueError(
                    "contour_field='potential' requires potential= (a nodal field "
                    "on the same tensor grid, or a callable V(r, R))"
                )
            if callable(potential):
                # Analytic potential: evaluate directly on the sampling meshgrid.
                r_mesh, big_r_mesh = np.meshgrid(a0, a1, indexing="ij")
                z = np.real(np.asarray(potential(r_mesh, big_r_mesh), dtype=float))
            else:
                # Numerically-evaluated potential known ONLY on the grid nodes:
                # project the nodal values through the same operator.
                z = np.real(projector.project_values(potential))
        levels = _contour_levels(contours, contour_field, mag, z)
        if levels:
            # Explicit (a1, a0) coords share the imshow-inverted axes, so the
            # lines land on the coloured features (imshow+contour alignment).
            ax.contour(
                a1,
                a0,
                z,
                levels=levels,
                colors=contour_color,
                alpha=contour_alpha,
                linewidths=contour_linewidth,
            )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if path is not None:
        ax.figure.tight_layout()
        ax.figure.savefig(path)
        if created:
            plt.close(ax.figure)
    return image
