"""`WavefunctionArtist` -- draw/update one HSV wavefunction panel in an axes.

The rendering is deliberately ISOLATED from figure creation: a
`WavefunctionArtist` is bound to a caller-supplied matplotlib ``Axes`` and knows
only how to draw a state into it and, via ``update``, restyle it for a new state.
So the same primitive serves a one-shot figure (`plot_wavefunction_2d`), a frame
of an animation (`animate_wavefunction`), and one panel of a composed multi-axes
figure -- the caller owns the figure and layout.

matplotlib is NOT imported at module scope (the ``Axes`` is passed in), so
``import qscat.viz`` stays matplotlib-free; only the higher-level entry points
that CREATE a figure import it (lazily).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from .coloring import complex_to_rgb
from .projector import EquidistantProjector

__all__ = ["WavefunctionArtist"]

# z-order: image at the bottom, potential turning surfaces above it, |psi| on top.
_Z_IMAGE = 0
_Z_POTENTIAL = 1
_Z_PSI = 2


def contour_levels(
    contours: bool | int | Sequence[float],
    contour_field: str,
    mag: float,
    z: npt.NDArray[np.float64],
) -> list[float]:
    """Resolve |psi|/primary contour levels; magnitude levels are derived from `mag`."""
    zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
    if contours is True or isinstance(contours, int):
        n = 20 if contours is True else int(contours)
        if contour_field == "magnitude":
            levels = [k * mag / 5.0 for k in range(1, n + 1)]
        else:
            if zmax <= zmin:
                return []
            levels = list(np.linspace(zmin, zmax, n + 2)[1:-1])
    else:
        levels = sorted(float(v) for v in contours)
    return [v for v in levels if zmin < v < zmax]


def potential_surface(
    projector: EquidistantProjector,
    potential: Any,
    a0: npt.NDArray[np.float64],
    a1: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """The potential on the sampling grid (nodal array -> project_values, or callable)."""
    if callable(potential):
        r_mesh, big_r_mesh = np.meshgrid(a0, a1, indexing="ij")
        return np.real(np.asarray(potential(r_mesh, big_r_mesh), dtype=float))
    return np.real(projector.project_values(potential))


def resolve_potential_levels(
    potential_levels: Sequence[float] | Literal["auto"],
    z: npt.NDArray[np.float64],
    eps: npt.ArrayLike | None,
    v_init: int,
    energies: npt.ArrayLike | None,
) -> list[float]:
    """Resolve the dotted potential-overlay levels; ``"auto"`` uses the energies."""
    from .levels import energy_contour_levels

    zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
    if isinstance(potential_levels, str):
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


class WavefunctionArtist:
    """Draw and update one domain-coloured wavefunction panel in a given axes.

    The potential turning-surface overlay (if any) is STATIC -- drawn once at
    construction. `update(state)` refreshes only the state-dependent layers (the
    domain-coloured image and the ``|psi|`` contours), so animation redraws the
    minimum per frame. See `plot_wavefunction_2d` for the parameter meanings; the
    style kwargs are identical, including the caveat that an array ``mag``'s
    magnitude contour levels key off its maximum, not its per-point values.
    """

    def __init__(
        self,
        ax: Any,
        projector: EquidistantProjector,
        *,
        mag: float | npt.NDArray[np.float64],
        inverse: bool = False,
        title: str | None = None,
        xlabel: str = "axis 1",
        ylabel: str = "axis 0",
        contours: bool | int | Sequence[float] = False,
        contour_field: Literal["magnitude", "potential"] = "magnitude",
        contour_color: str = "white",
        contour_alpha: float = 0.6,
        contour_linewidth: float = 0.6,
        potential: Any = None,
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
    ) -> None:
        self.ax = ax
        self.projector = projector
        self.mag = mag
        # Contour levels must be scalars, so they key off the brightest
        # region's scale; the per-point array still drives the image.
        self._contour_mag = float(np.max(np.asarray(mag, dtype=np.float64)))
        self.inverse = inverse
        self.contours = contours
        self.contour_field = contour_field
        self._contour_style = dict(
            colors=contour_color, alpha=contour_alpha, linewidths=contour_linewidth
        )
        self._potential = potential

        a0, a1 = projector.axis0, projector.axis1
        self._a0, self._a1 = a0, a1
        self._extent = [float(a1[0]), float(a1[-1]), float(a0[-1]), float(a0[0])]
        self._image: Any = None
        self._psi_contours: Any = None

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title)

        # Static potential turning-surface overlay (state-independent) -- once.
        if potential is not None and potential_levels is not None:
            zp = potential_surface(projector, potential, a0, a1)
            plevels = resolve_potential_levels(potential_levels, zp, eps, v_init, energies)
            if plevels:
                pset = ax.contour(
                    a1,
                    a0,
                    zp,
                    levels=plevels,
                    colors=potential_color,
                    linestyles=potential_style,
                    alpha=potential_alpha,
                    linewidths=potential_linewidth,
                    zorder=_Z_POTENTIAL,
                )
                if potential_labels:
                    ax.clabel(pset, fmt=potential_label_fmt, fontsize=7)

    def update(self, state: npt.NDArray[np.complex128], *, phase: float = 0.0) -> list[Any]:
        """Refresh the image + ``|psi|`` contours for a new state; return changed artists.

        `phase` (radians) applies a global hue rotation ``psi -> e^{i*phase}*psi``
        to the COLOURING only -- used to view the phase relative to a reference
        energy (e.g. ``phase = E_ref * t`` removes the channel base-energy spin).
        Since it is a global phase it leaves ``|psi|`` (brightness, contours,
        the potential overlay) unchanged, and by linearity of the projector it is
        equivalent to phasing the state -- applied here on the small projected
        field, so it costs nothing extra over the per-frame recolour.
        """
        field = self.projector.project(state)
        colour_field = field if phase == 0.0 else field * np.exp(1j * phase)
        rgb = complex_to_rgb(colour_field, self.mag, inverse=self.inverse)
        if self._image is None:
            self._image = self.ax.imshow(
                rgb, origin="upper", aspect="auto", extent=self._extent, zorder=_Z_IMAGE
            )
        else:
            self._image.set_data(rgb)
        changed: list[Any] = [self._image]

        if self.contours is not False:
            if self._psi_contours is not None:
                self._psi_contours.remove()  # QuadContourSet.remove()
                self._psi_contours = None
            if self.contour_field == "magnitude":
                z = np.abs(field)
            elif self._potential is None:
                raise ValueError(
                    "contour_field='potential' requires potential= (a nodal field "
                    "on the same tensor grid, or a callable V(r, R))"
                )
            else:
                z = potential_surface(self.projector, self._potential, self._a0, self._a1)
            levels = contour_levels(self.contours, self.contour_field, self._contour_mag, z)
            if levels:
                self._psi_contours = self.ax.contour(
                    self._a1, self._a0, z, levels=levels, zorder=_Z_PSI, **self._contour_style
                )
                changed.append(self._psi_contours)
        return changed
