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


def _potential_surface(
    projector: EquidistantProjector,
    potential: Any,
    a0: npt.NDArray[np.float64],
    a1: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """The potential on the sampling grid (can be negative -> real, not abs)."""
    if callable(potential):
        # Analytic potential: evaluate directly on the sampling meshgrid.
        r_mesh, big_r_mesh = np.meshgrid(a0, a1, indexing="ij")
        return np.real(np.asarray(potential(r_mesh, big_r_mesh), dtype=float))
    # Numerically-evaluated potential known only on the grid: project as values.
    return np.real(projector.project_values(potential))


def _resolve_potential_levels(
    potential_levels: Sequence[float] | Literal["auto"],
    z: npt.NDArray[np.float64],
    eps: npt.ArrayLike | None,
    v_init: int,
    energies: npt.ArrayLike | None,
) -> list[float]:
    """Resolve the dotted potential-overlay levels; ``"auto"`` uses the energies."""
    from .levels import energy_contour_levels

    zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
    if isinstance(potential_levels, str):  # "auto"
        if eps is None and energies is None:
            raise ValueError(
                "potential_levels='auto' needs eps= (vibrational energies) and/or "
                "energies= (collision energies) to place the turning-surface levels"
            )
        levels = energy_contour_levels(
            eps=eps, v_init=v_init, energies=energies, e_range=(zmin, zmax)
        )
    else:
        levels = sorted(float(v) for v in potential_levels)
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
    potential_levels: Sequence[float] | Literal["auto"] | None = None,
    potential_style: str = ":",
    potential_color: str = "0.75",
    potential_alpha: float = 0.7,
    potential_linewidth: float = 0.6,
    potential_labels: bool = True,
    potential_label_fmt: str = "%.3f",
    eps: npt.ArrayLike | None = None,
    v_init: int = 0,
    energies: npt.ArrayLike | None = None,
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
        |psi|-contour style; defaults are thin white lines at 0.6 opacity.
    potential_levels : sequence of float or "auto", optional
        Enables the DEDICATED dotted potential overlay (drawn in addition to the
        ``|psi|`` contours) when both this and ``potential`` are given. An
        explicit list of energies (Hartree), or ``"auto"`` to derive turning-
        surface levels from ``eps``/``energies`` via `energy_contour_levels`.
        The physical (turning-surface) reading holds only if ``potential`` is the
        FULL 2-D PES in Hartree -- pass ``model.surface`` (``v0(R) +
        ell(ell+1)/2r^2 + v_int(r,R)``), NOT just the nuclear ``v0`` or the
        interaction-only ``interaction_diag``.
    potential_style, potential_color, potential_alpha, potential_linewidth : optional
        Style of the potential overlay; defaults are dotted grey lines at 0.7 α.
    potential_labels : bool, optional
        Inline-label each dotted line with its energy (``potential_label_fmt``).
    eps, v_init, energies : optional
        For ``potential_levels="auto"``: the vibrational energies, initial level,
        and collision energies (levels ``eps_v`` and ``eps[v_init] + E``).

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

    # Primary overlay: |psi| (default) or potential-as-primary.
    if contours is not False:
        if contour_field == "magnitude":
            z = np.abs(field)  # |psi|, contoured on the same sampling grid
        else:  # potential-as-primary (can be negative -> _potential_surface uses .real)
            if potential is None:
                raise ValueError(
                    "contour_field='potential' requires potential= (a nodal field "
                    "on the same tensor grid, or a callable V(r, R))"
                )
            z = _potential_surface(projector, potential, a0, a1)
        levels = _contour_levels(contours, contour_field, mag, z)
        if levels:
            # Explicit (a1, a0) coords share the imshow-inverted axes, so the
            # lines land on the coloured features (imshow+contour alignment).
            ax.contour(
                a1, a0, z, levels=levels, colors=contour_color,
                alpha=contour_alpha, linewidths=contour_linewidth,
            )

    # Dedicated dotted potential overlay at energy-relevant (turning-surface)
    # levels -- drawn IN ADDITION to the |psi| contours above.
    if potential is not None and potential_levels is not None:
        zp = _potential_surface(projector, potential, a0, a1)
        plevels = _resolve_potential_levels(potential_levels, zp, eps, v_init, energies)
        if plevels:
            pset = ax.contour(
                a1, a0, zp, levels=plevels, colors=potential_color,
                linestyles=potential_style, alpha=potential_alpha,
                linewidths=potential_linewidth,
            )
            if potential_labels:
                ax.clabel(pset, fmt=potential_label_fmt, fontsize=7)

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
