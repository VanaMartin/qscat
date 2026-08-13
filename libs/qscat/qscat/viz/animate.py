"""Animate a time-dependent 2-D wavefunction (domain-coloured frames).

Decoupled from the physics: `animate_wavefunction` takes a sequence of state
vectors (the propagated ``psi(t)`` snapshots you produced) and drives a
`WavefunctionArtist` frame by frame. The magnitude scale and the potential
turning-surface overlay are FIXED across frames (one brightness scale, static
contours); only the domain-coloured image and the ``|psi|`` contours update.

`animate_artists` animates several panels (each its own `WavefunctionArtist` and
frame sequence) in one figure, for composed / side-by-side views.

Output: ``.mp4`` (ffmpeg) or ``.gif`` (pillow), picked by extension; or return
the `FuncAnimation` unsaved. matplotlib is imported lazily.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import numpy.typing as npt

from .artist import WavefunctionArtist
from .projector import EquidistantProjector

__all__ = ["animate_wavefunction", "animate_artists"]

_PathLike = str | os.PathLike[str]


def _lazy_mpl() -> Any:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ModuleNotFoundError as exc:  # pragma: no cover - trivial guard
        raise ModuleNotFoundError(
            "qscat.viz animation requires matplotlib. Install the plotting extra: "
            "pip install 'qscat[plot]'."
        ) from exc


def _pick_writer(outfile: _PathLike, fps: int, writer: Any) -> Any:
    """Choose a matplotlib animation writer by file extension (or pass one through)."""
    if writer is not None:
        return writer
    from matplotlib import animation

    ext = os.path.splitext(str(outfile))[1].lower()
    if ext == ".gif":
        return animation.PillowWriter(fps=fps)
    if ext in (".mp4", ".m4v", ".mov"):
        if not animation.FFMpegWriter.isAvailable():
            raise RuntimeError(
                f"writing {ext} needs ffmpeg on PATH (not found). Install ffmpeg, "
                "or use a .gif output (pillow), or pass writer=."
            )
        return animation.FFMpegWriter(fps=fps)
    raise ValueError(f"unsupported animation extension {ext!r}; use .mp4 or .gif")


def _save(anim: Any, outfile: _PathLike | None, fps: int, writer: Any) -> None:
    if outfile is not None:
        anim.save(str(outfile), writer=_pick_writer(outfile, fps, writer))


def animate_wavefunction(
    projector: EquidistantProjector,
    frames: Iterable[npt.NDArray[np.complex128]],
    *,
    mag: float,
    times: Sequence[float] | None = None,
    time_fmt: str = "t = {:.1f}",
    outfile: _PathLike | None = None,
    fps: int = 15,
    writer: Any = None,
    ax: Any = None,
    figsize: tuple[float, float] = (8, 6),
    title: str | None = None,
    **style: Any,
) -> Any:
    """Animate a sequence of states through one `WavefunctionArtist`.

    Parameters
    ----------
    projector : EquidistantProjector
        Cached projector for the states' tensor grid.
    frames : iterable of ndarray
        The states ``psi(t)`` to animate, one per frame.
    mag : float
        Fixed brightness scale across all frames.
    times : sequence of float, optional
        Per-frame times; when given the title shows ``time_fmt.format(t)``.
    outfile : path-like, optional
        Save target: ``.mp4`` (ffmpeg) or ``.gif`` (pillow). If None, the
        `FuncAnimation` is returned unsaved.
    fps, writer, ax, figsize, title : optional
        Frame rate, an explicit matplotlib writer, an existing Axes (else a new
        figure), figure size, and a base title.
    **style
        Forwarded to `WavefunctionArtist` (contours/potential/colours/…).

    Returns
    -------
    matplotlib.animation.FuncAnimation
    """
    plt = _lazy_mpl()
    from matplotlib.animation import FuncAnimation

    frame_list = [np.asarray(f) for f in frames]
    if not frame_list:
        raise ValueError("frames is empty; nothing to animate")

    created = ax is None
    if created:
        _, ax = plt.subplots(figsize=figsize)
    artist = WavefunctionArtist(ax, projector, mag=mag, title=title, **style)

    def _update(i: int) -> list[Any]:
        arts = artist.update(frame_list[i])
        if times is not None:
            label = time_fmt.format(times[i])
            ax.set_title(f"{title}  {label}" if title else label)
        return arts

    anim = FuncAnimation(ax.figure, _update, frames=len(frame_list), blit=False)
    _save(anim, outfile, fps, writer)
    if outfile is not None and created:
        plt.close(ax.figure)
    return anim


def animate_artists(
    fig: Any,
    panels: Sequence[tuple[WavefunctionArtist, Sequence[npt.NDArray[np.complex128]]]],
    *,
    outfile: _PathLike | None = None,
    fps: int = 15,
    writer: Any = None,
) -> Any:
    """Animate several panels (artist + its frames) together in one figure.

    Each panel is ``(artist, frames)``; all are advanced in lockstep (the frame
    count is the shortest panel's). Use this for composed / side-by-side views:
    build a multi-axes figure, a `WavefunctionArtist` per axes, and pass them all.
    """
    _lazy_mpl()
    from matplotlib.animation import FuncAnimation

    prepared = [(a, [np.asarray(f) for f in seq]) for a, seq in panels]
    if not prepared:
        raise ValueError("panels is empty")
    n = min(len(seq) for _, seq in prepared)

    def _update(i: int) -> list[Any]:
        arts: list[Any] = []
        for artist, seq in prepared:
            arts += artist.update(seq[i])
        return arts

    anim = FuncAnimation(fig, _update, frames=n, blit=False)
    _save(anim, outfile, fps, writer)
    return anim
