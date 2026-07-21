"""General two-spectrum resonance-pole matcher.

Promoted from `projects/n2_resonance/pole.find_pole` (sub-project #2, Task 2)
once validated there: the physics-specific part (assembling `H_el(R)` on the
N2 grids) stays in `projects/n2_resonance/pole.py`, while the pure
eigenvalue-matching core -- generic to any two-angle ECS resonance search --
lives here so other projects can reuse it.

See `docs/physics/femdvr-ecs.md` and `docs/physics/n2-resonance.md` for the
method: an ECS-rotated Hamiltonian's discretized continuum eigenvalues rotate
with the rotation angle `theta`, while a true resonance pole is (nearly)
`theta`-independent. Diagonalizing the same physical Hamiltonian at two
different ECS angles and finding the pair of eigenvalues -- one from each
spectrum -- that agree most closely therefore isolates the resonance.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["find_resonance_pole"]


def _filter_window(
    E: npt.NDArray[np.complex128], window: tuple[float, float, float, float]
) -> npt.NDArray[np.complex128]:
    re_lo, re_hi, im_lo, im_hi = window
    mask = (E.real >= re_lo) & (E.real <= re_hi) & (E.imag >= im_lo) & (E.imag <= im_hi)
    return E[mask]


def find_resonance_pole(
    eigs_a: npt.ArrayLike,
    eigs_b: npt.ArrayLike,
    window: tuple[float, float, float, float],
) -> tuple[complex, float]:
    """Match the angle-stable resonance pole between two eigenvalue spectra.

    `eigs_a`/`eigs_b` are complex eigenvalue arrays of the same Hamiltonian
    computed at two different ECS rotation angles (or otherwise perturbed in
    a way that moves the discretized continuum but not a true pole). Both
    spectra are restricted to `window = (re_lo, re_hi, im_lo, im_hi)`, and the
    pair `(ea, eb)` -- one eigenvalue from each, restricted set -- with the
    smallest `|ea - eb|` is returned as the matched pole: `E_pole =
    0.5*(ea+eb)`, `residual = |ea-eb|`. A small residual (<< the resonance
    width) is the signature of a genuine angle-stable pole; discretized
    continuum eigenvalues rotate with the angle and do not match this
    closely.

    Raises `ValueError` if `window` contains no eigenvalues in either input
    spectrum (window too tight, or grid too coarse to resolve the pole).
    """
    fa = _filter_window(np.asarray(eigs_a, dtype=np.complex128), window)
    fb = _filter_window(np.asarray(eigs_b, dtype=np.complex128), window)

    if fa.size == 0 or fb.size == 0:
        raise ValueError(
            f"find_resonance_pole: window {window} contains no eigenvalues in "
            f"{'eigs_a' if fa.size == 0 else 'eigs_b'} "
            f"(found {fa.size} in A, {fb.size} in B) -- window too tight "
            "or spectrum too coarse."
        )

    # For each candidate in A, distance to nearest candidate in B.
    diffs = np.abs(fa[:, None] - fb[None, :])
    idx = np.unravel_index(np.argmin(diffs), diffs.shape)
    i, j = int(idx[0]), int(idx[1])
    ea, eb = fa[i], fb[j]
    residual = float(np.abs(ea - eb))
    E_pole = complex(0.5 * (ea + eb))
    return E_pole, residual
