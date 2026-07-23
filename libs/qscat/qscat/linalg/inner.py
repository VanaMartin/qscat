"""The c-product: the bilinear inner product exterior complex scaling requires.

Under ECS the Hamiltonian is complex SYMMETRIC (`H = H^T`), not Hermitian, so
the natural pairing is the bilinear `sum_i a_i b_i` with NO complex conjugate
-- not `numpy.vdot`'s sesquilinear `sum_i conj(a_i) b_i`.

Getting this wrong is a recurring, quiet failure mode: it produces
plausible-looking complex "cross sections" with the wrong phase rather than an
obvious error. It has already bitten this repo once (sub-project #3's S-matrix,
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
    """`sum_i a_i b_i` -- the bilinear (NOT conjugated) inner product."""
    av = np.asarray(a, dtype=np.complex128).ravel()
    bv = np.asarray(b, dtype=np.complex128).ravel()
    if av.shape != bv.shape:
        raise ValueError(f"shape mismatch: {av.shape} vs {bv.shape}")
    return complex(np.dot(av, bv))
