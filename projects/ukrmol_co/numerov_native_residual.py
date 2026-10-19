"""Experimental residual primitive with retained NumPy differential reference."""

from __future__ import annotations

import sys

import numpy as np


def numpy_residual(a, b, x):
    """Original same-dtype product then subtraction, without symmetrization."""
    return b - a @ x


def native_residual(a, b, x):
    """Native64-significand-bit residual; preserve array values, shapes and dtype."""
    arrays = [np.asarray(value) for value in (a, b, x)]
    a, b, x = arrays
    if (
        any(value.ndim != 2 or min(value.shape) == 0 for value in arrays)
        or a.shape[1] != x.shape[0]
        or b.shape != (a.shape[0], x.shape[1])
        or any(value.dtype != a.dtype or not np.isfinite(value).all() for value in arrays)
    ):
        raise ValueError("Finite compatible same-dtype matrices required")
    if (
        a.dtype != np.dtype(np.longdouble)
        or not a.dtype.isnative
        or np.finfo(a.dtype).nmant + 1 != 64
        or a.dtype.itemsize != 16
        or sys.byteorder != "little"
    ):
        raise ValueError("Native residual requires the qualified extended scalar format")
    import qscat_extended_kernels as native

    if tuple(native.capability()) != (16, 64, 16384, True):
        raise ValueError("Native and NumPy scalar formats differ")
    result = native.residual_bytes(
        *(np.ascontiguousarray(value).tobytes() for value in arrays),
        a.shape[0],
        a.shape[1],
        x.shape[1],
    )
    return np.frombuffer(result, dtype=a.dtype).reshape(b.shape).copy()
