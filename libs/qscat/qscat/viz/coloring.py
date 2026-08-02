"""Domain coloring: map a complex field to RGB (phase -> hue, magnitude -> value).

Ported from eMoScat's ``display_wf.py`` (``complex_to_hsv``/``complex_to_rgb``).
The mapping is the standard complex-plane "domain coloring":

  * hue     = arg(z) / 2pi           (the phase; 0 at zero angle)
  * value   = |z| / Mag              (brightness grows to full at |z| = Mag)
  * for |z| > Mag the colour desaturates toward white instead of clipping

``Mag`` sets the magnitude that maps to full brightness. ``inverse=True`` swaps
the roles so large magnitudes go dark on a white background (for light figures).

Pure numpy (no matplotlib) so the RGB array can be produced and tested without
the plotting extra; `qscat.viz.plot` renders it.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["complex_to_hsv", "complex_to_rgb", "hsv_to_rgb"]


def complex_to_hsv(
    z: npt.ArrayLike, mag: float = 1.0, *, inverse: bool = False
) -> npt.NDArray[np.float64]:
    """Complex array -> HSV array (shape ``z.shape + (3,)``), all channels in [0, 1]."""
    c = np.asarray(z)
    if not np.iscomplexobj(c):
        raise ValueError(f"input must be complex, got dtype {c.dtype}")

    hsv = np.zeros(c.shape + (3,), dtype=np.float64)
    # Hue from the phase, wrapped from [-pi, pi) to [0, 1).
    h = np.angle(c) / (2.0 * np.pi)
    hsv[..., 0] = np.where(h < 0.0, 1.0 + h, h)

    r = np.abs(c) / mag
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
    z: npt.ArrayLike, mag: float = 1.0, *, inverse: bool = False
) -> npt.NDArray[np.float64]:
    """Complex array -> RGB array (shape ``z.shape + (3,)``), channels in [0, 1]."""
    return hsv_to_rgb(complex_to_hsv(z, mag, inverse=inverse))
