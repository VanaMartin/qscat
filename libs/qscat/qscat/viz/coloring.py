"""Domain coloring: map a complex field to RGB (phase -> hue, magnitude -> value).

Ported from eMoScat's ``display_wf.py`` (``complex_to_hsv``/``complex_to_rgb``).
The mapping is the standard complex-plane "domain coloring":

  * hue     = arg(z) / 2pi           (the phase; 0 at zero angle)
  * value   = |z| / Mag              (brightness grows to full at |z| = Mag)
  * for |z| > Mag the colour desaturates toward white instead of clipping

``Mag`` sets the magnitude that maps to full brightness -- a scalar, or a
per-point array the same shape as the field (see `region_magnitudes`, which
builds one by normalising disjoint regions of the field independently, so a
region orders of magnitude weaker than the field's brightest feature still
renders visibly instead of black). ``inverse=True`` swaps the roles so large
magnitudes go dark on a white background (for light figures).

Pure numpy (no matplotlib) so the RGB array can be produced and tested without
the plotting extra; `qscat.viz.plot` renders it.

TODO (print mode): add an explicit inverse VALUE (brightness) scaling -- map the
value channel ``black->white`` to ``white->black`` so features render dark on a
white ground, useful for printed figures. This is distinct from the existing
``inverse`` (which swaps the saturation/value roles); it is a straight
brightness inversion (``v -> 1 - v``, with saturation handled to keep hue
readable). Deferred to a follow-on; see the roadmap.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

__all__ = ["complex_to_hsv", "complex_to_rgb", "hsv_to_rgb", "region_magnitudes"]


def region_magnitudes(
    magnitude: npt.ArrayLike,
    *,
    axis: int,
    boundaries: Sequence[int],
    percentile: float = 99.5,
    floor: float = 1e-300,
) -> npt.NDArray[np.float64]:
    """Per-point brightness scale that normalises each region separately.

    A single global ``mag`` ties the whole field to its largest feature, so a
    region whose amplitude is orders of magnitude smaller renders black. This
    splits ``axis`` at ``boundaries`` and gives every point the scale of its own
    region, so each region is visible on its own terms.

    The scale is a percentile rather than the maximum, so one outlying point
    cannot flatten the region it sits in.

    Parameters
    ----------
    magnitude : array_like
        Non-negative field, e.g. ``np.abs(psi)``.
    axis : int
        Axis the split runs along.
    boundaries : sequence of int
        Strictly increasing split indices along ``axis``, each in
        ``1 .. n-1``. An empty sequence means one region (the whole field).
    percentile : float
        Percentile of each region's magnitudes mapped to full brightness.
    floor : float
        Lower clamp, so an all-zero region yields a positive scale rather than
        a divide-by-zero downstream.

    Returns
    -------
    ndarray
        Same shape as ``magnitude``; every point carries its region's scale.
        Pass it straight to ``complex_to_rgb`` as ``mag``.
    """
    m = np.asarray(magnitude, dtype=np.float64)
    n = m.shape[axis]
    bounds = [int(b) for b in boundaries]
    if any(b < 1 or b > n - 1 for b in bounds):
        raise ValueError(
            f"boundaries must lie in 1..{n - 1} along axis {axis} (length {n}); got {bounds}"
        )
    # bounds/bounds[1:] and edges/edges[1:] are sliding-window pairs, always
    # one element shorter on the right by construction -- strict=False, not
    # an omission.
    if any(b >= c for b, c in zip(bounds, bounds[1:], strict=False)):
        raise ValueError(f"boundaries must be strictly increasing; got {bounds}")

    edges = [0, *bounds, n]
    out = np.empty_like(m)
    for lo, hi in zip(edges, edges[1:], strict=False):
        sl: list[slice] = [slice(None)] * m.ndim
        sl[axis] = slice(lo, hi)
        block = m[tuple(sl)]
        scale = float(np.percentile(block, percentile)) if block.size else 0.0
        out[tuple(sl)] = max(scale, floor)
    return out


def complex_to_hsv(
    z: npt.ArrayLike, mag: float | npt.NDArray[np.float64] = 1.0, *, inverse: bool = False
) -> npt.NDArray[np.float64]:
    """Complex array -> HSV array (shape ``z.shape + (3,)``), all channels in [0, 1]."""
    c = np.asarray(z)
    if not np.iscomplexobj(c):
        raise ValueError(f"input must be complex, got dtype {c.dtype}")

    hsv = np.zeros(c.shape + (3,), dtype=np.float64)
    # Hue from the phase, wrapped from [-pi, pi) to [0, 1).
    h = np.angle(c) / (2.0 * np.pi)
    hsv[..., 0] = np.where(h < 0.0, 1.0 + h, h)

    scale = np.asarray(mag, dtype=np.float64)
    if scale.ndim and scale.shape != c.shape:
        try:
            scale = np.broadcast_to(scale, c.shape)
        except ValueError as exc:
            raise ValueError(
                f"mag of shape {scale.shape} is not broadcastable to the field shape {c.shape}"
            ) from exc
    r = np.abs(c) / scale
    if inverse:
        hsv[..., 1] = np.minimum(1.0, r)  # saturation
        hsv[..., 2] = 1.0 / np.maximum(1.0, r)  # value (dark where large)
    else:
        hsv[..., 1] = 1.0 / np.maximum(1.0, r)  # saturation (white-washes when large)
        hsv[..., 2] = np.minimum(1.0, r)  # value (dark where small)
    return hsv


def hsv_to_rgb(hsv: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Vectorised HSV -> RGB (channels last), matching matplotlib's convention."""
    h = hsv[..., 0] * 6.0
    s = hsv[..., 1]
    v = hsv[..., 2]
    i = np.floor(h).astype(int)
    f = h - i
    i = i % 6
    p = v * (1.0 - s)
    q = v * (1.0 - s * f)
    t = v * (1.0 - s * (1.0 - f))
    rgb = np.zeros(hsv.shape, dtype=np.float64)
    conditions = [
        (v, t, p),
        (q, v, p),
        (p, v, t),
        (p, q, v),
        (t, p, v),
        (v, p, q),
    ]
    for idx, (rr, gg, bb) in enumerate(conditions):
        m = i == idx
        rgb[m, 0], rgb[m, 1], rgb[m, 2] = rr[m], gg[m], bb[m]
    return rgb


def complex_to_rgb(
    z: npt.ArrayLike, mag: float | npt.NDArray[np.float64] = 1.0, *, inverse: bool = False
) -> npt.NDArray[np.float64]:
    """Complex array -> RGB array (shape ``z.shape + (3,)``), channels in [0, 1]."""
    return hsv_to_rgb(complex_to_hsv(z, mag, inverse=inverse))
