"""Render a projected 2-D wavefunction as a domain-coloured image.

matplotlib is imported lazily (the optional ``qscat[plot]`` extra), mirroring
`qscat.core.plot`. `plot_wavefunction_2d` projects a state through an
`EquidistantProjector`, domain-colours it (`complex_to_rgb`), and draws/saves the
image -- the single-frame building block for animations.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import numpy.typing as npt

from .coloring import complex_to_rgb
from .projector import EquidistantProjector

__all__ = ["plot_wavefunction_2d"]

_PathLike = str | os.PathLike[str]


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
