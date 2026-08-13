"""One-shot render of a 2-D wavefunction as a domain-coloured image.

`plot_wavefunction_2d` is a thin wrapper over `WavefunctionArtist` (which does
the actual drawing into an axes): it creates a figure/axes when needed, draws one
state, and optionally saves. For composed multi-panel figures or animations, use
`WavefunctionArtist` directly (see `qscat.viz.animate_wavefunction`).

matplotlib is imported lazily (the optional ``qscat[plot]`` extra).
"""

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from .artist import WavefunctionArtist
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
    phase: float = 0.0,
) -> Any:
    """Project, domain-colour, and draw a 2-D state; save to ``path`` if given.

    A convenience wrapper: it makes a figure/axes if ``ax`` is None, builds a
    `WavefunctionArtist`, draws ``state`` once, and returns the image handle.

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
        Light-background colour mapping (swaps saturation/value). Distinct from
        the pending print-mode *value* inversion (see `complex_to_hsv`).
    title, xlabel, ylabel : str, optional
        Axis title/labels; ``axis 0`` is the projector's first grid, ``axis 1``
        the second (e.g. ``ylabel="electronic r"``, ``xlabel="nuclear R"``).
    ax : matplotlib Axes, optional
        Draw into an existing Axes (composition); a new figure is made when None.
    contours : bool or int or sequence of float, optional
        ``|psi|`` contour overlay: ``False`` (default) none; ``True`` default
        levels; ``int`` that many; a sequence explicit levels.
    contour_field : {"magnitude", "potential"}, optional
        Primary contour source: ``|psi|`` (levels ``k*mag/5``) or the
        ``potential`` field.
    potential : ndarray or callable, optional
        The potential field: nodal values on the SAME tensor grid (projected via
        `EquidistantProjector.project_values`) or a callable ``V(r, R)``. For the
        FULL 2-D PES pass ``model.surface`` (``v0(R) + ell(ell+1)/2r^2 +
        v_int(r,R)``), NOT ``v0`` or ``interaction_diag`` alone.
    contour_color, contour_alpha, contour_linewidth : optional
        ``|psi|``-contour style; defaults thin white lines at 0.6 opacity.
    potential_levels : sequence of float or "auto", optional
        Enables the dotted potential overlay (in addition to ``|psi|``) when this
        and ``potential`` are given. Explicit energies (Hartree) or ``"auto"``
        (turning surfaces via `energy_contour_levels` from ``eps``/``energies``).
    potential_style, potential_color, potential_alpha, potential_linewidth : optional
        Potential-overlay style; defaults dotted grey at 0.7 alpha.
    potential_labels : bool, optional
        Inline-label each dotted line with its energy (``potential_label_fmt``).
    eps, v_init, energies : optional
        For ``potential_levels="auto"``: vibrational energies, initial level, and
        collision energies (levels ``eps_v`` and ``eps[v_init] + E``).
    phase : float, optional
        Global hue-rotation (radians) applied to the colouring only
        (``psi -> e^{i*phase}*psi``); leaves ``|psi|`` unchanged. For animation,
        `animate_wavefunction`'s ``phase_reference`` sets ``phase = E_ref*t`` to
        view the phase relative to a channel base energy.

    Returns
    -------
    matplotlib.image.AxesImage
        The drawn image handle.

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

    created = ax is None
    if created:
        _, ax = plt.subplots(figsize=(8, 6))

    artist = WavefunctionArtist(
        ax, projector, mag=mag, inverse=inverse, title=title, xlabel=xlabel,
        ylabel=ylabel, contours=contours, contour_field=contour_field,
        contour_color=contour_color, contour_alpha=contour_alpha,
        contour_linewidth=contour_linewidth, potential=potential,
        potential_levels=potential_levels, potential_style=potential_style,
        potential_color=potential_color, potential_alpha=potential_alpha,
        potential_linewidth=potential_linewidth, potential_labels=potential_labels,
        potential_label_fmt=potential_label_fmt, eps=eps, v_init=v_init,
        energies=energies,
    )
    (image, *_) = artist.update(np.asarray(state), phase=phase)

    if path is not None:
        ax.figure.tight_layout()
        ax.figure.savefig(path)
        if created:
            plt.close(ax.figure)
    return image
