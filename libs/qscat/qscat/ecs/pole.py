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

__all__ = ["find_resonance_pole", "match_angle_stable"]


def _window_indices(
    E: npt.NDArray[np.complex128], window: tuple[float, float, float, float]
) -> npt.NDArray[np.intp]:
    re_lo, re_hi, im_lo, im_hi = window
    mask = (E.real >= re_lo) & (E.real <= re_hi) & (E.imag >= im_lo) & (E.imag <= im_hi)
    return np.flatnonzero(mask)


def _paired(
    eigs_a: npt.ArrayLike,
    eigs_b: npt.ArrayLike,
    window: tuple[float, float, float, float],
    caller: str,
) -> tuple[
    npt.NDArray[np.complex128],
    npt.NDArray[np.complex128],
    npt.NDArray[np.intp],
    npt.NDArray[np.intp],
    npt.NDArray[np.float64],
]:
    """Window-filter both spectra and pair each surviving `a` with its nearest `b`.

    Returns `(fa, fb, ia, nearest, dist)`: the windowed values, `ia` their
    indices into the ORIGINAL `eigs_a`, `nearest[k]` the index into `fb` closest
    to `fa[k]`, and `dist[k]` that distance. Shared by `find_resonance_pole`
    (which takes the global argmin) and `match_angle_stable` (which applies a
    tolerance cut) so there is exactly one implementation of the criterion.
    """
    a = np.asarray(eigs_a, dtype=np.complex128)
    b = np.asarray(eigs_b, dtype=np.complex128)
    ia = _window_indices(a, window)
    ib = _window_indices(b, window)
    if ia.size == 0 or ib.size == 0:
        raise ValueError(
            f"{caller}: window {window} contains no eigenvalues in "
            f"{'eigs_a' if ia.size == 0 else 'eigs_b'} "
            f"(found {ia.size} in A, {ib.size} in B) -- window too tight "
            "or spectrum too coarse."
        )
    fa, fb = a[ia], b[ib]
    diffs = np.abs(fa[:, None] - fb[None, :])
    nearest = np.asarray(np.argmin(diffs, axis=1), dtype=np.intp)
    dist = np.asarray(diffs[np.arange(fa.size), nearest], dtype=np.float64)
    return fa, fb, ia, nearest, dist


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
    fa, fb, _ia, nearest, dist = _paired(eigs_a, eigs_b, window, "find_resonance_pole")
    i = int(np.argmin(dist))
    ea, eb = fa[i], fb[nearest[i]]
    return complex(0.5 * (ea + eb)), float(np.abs(ea - eb))


def match_angle_stable(
    eigs_a: npt.ArrayLike,
    eigs_b: npt.ArrayLike,
    window: tuple[float, float, float, float],
    *,
    rel_tol: float = 1e-4,
    atol: float = 1e-8,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64], npt.NDArray[np.intp]]:
    """Every angle-stable eigenvalue shared by two ECS spectra, not just one.

    The multi-state generalization of `find_resonance_pole`: `eigs_a`/`eigs_b`
    are complex eigenvalue arrays of the same Hamiltonian at two different ECS
    rotation angles. An eigenvalue of `eigs_a` inside `window` is ACCEPTED when
    its nearest `eigs_b` partner satisfies

        ``|E_a - E_b| < max(rel_tol * |E_a|, atol)``

    -- eMoScat's `DiscreteStates` criterion, vectorized. Discretized continuum
    eigenvalues rotate with the angle and fail it; bound and resonance states do
    not. Returns `(energies, residuals, indices)`, ascending in `Re E`:
    `energies` the midpoints `(E_a + E_b)/2` (matching `find_resonance_pole`'s
    convention), `residuals` the `|E_a - E_b|` per accepted state, and `indices`
    the positions in the ORIGINAL `eigs_a` -- so a caller holding the grid-`a`
    eigenvectors can pull the matching columns straight out.

    An empty result is a normal outcome (no stable state in `window`), NOT an
    error. `ValueError` is raised only when `window` catches nothing at all in
    one of the two spectra, mirroring `find_resonance_pole`.
    """
    fa, fb, ia, nearest, dist = _paired(eigs_a, eigs_b, window, "match_angle_stable")
    keep = dist < np.maximum(rel_tol * np.abs(fa), atol)
    energies = 0.5 * (fa[keep] + fb[nearest[keep]])
    residuals = dist[keep]
    indices = ia[keep]
    order = np.argsort(energies.real)
    return (
        np.asarray(energies[order], dtype=np.complex128),
        np.asarray(residuals[order], dtype=np.float64),
        np.asarray(indices[order], dtype=np.intp),
    )
