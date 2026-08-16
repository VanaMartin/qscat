"""The c-product: the bilinear inner product exterior complex scaling requires.

Under ECS the Hamiltonian is complex SYMMETRIC (`H = H^T`), not Hermitian, so
the natural pairing is the bilinear `sum_i a_i b_i` with NO complex conjugate
-- not `numpy.vdot`'s sesquilinear `sum_i conj(a_i) b_i`.

Getting this wrong is a recurring, quiet failure mode: it produces
plausible-looking complex "cross sections" with the wrong phase rather than an
obvious error. It has already bitten this repo once (the N2 S-matrix,
where the Hermitian convention gave negative sigma), and the reference C++ code
uses `cblas_zdotc` here -- formally wrong, and correct in practice only because
it zeroes every channel function on the complex-scaled tail. Naming the
operation makes the choice explicit at every call site.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["c_product"]


def c_product(a: npt.ArrayLike, b: npt.ArrayLike) -> complex:
    """`sum_i a_i b_i` -- the bilinear (NOT conjugated) inner product.

    Shapes are compared BEFORE flattening: `c_product(psi_(n0, n1),
    chi_(n1, n0))` with `n0 != n1` raises rather than silently ravelling both
    down to the same total element count and returning a plausible-looking
    but physically wrong number (a transposed-axis bug hiding behind a
    reshape that `ravel()` alone would never catch, since `ravel` only cares
    about total size, not per-axis shape).
    """
    a_arr = np.asarray(a, dtype=np.complex128)
    b_arr = np.asarray(b, dtype=np.complex128)
    if a_arr.shape != b_arr.shape:
        raise ValueError(f"shape mismatch: {a_arr.shape} vs {b_arr.shape}")
    return complex(np.dot(a_arr.ravel(), b_arr.ravel()))
