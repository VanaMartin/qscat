"""Choose physically relevant potential-contour levels for a 2-D wavefunction.

A state's total energy is ``E_tot = eps_v + E_collision``; the equipotential
``V(r,R) = E_tot`` is the classical turning surface. So the useful potential
levels to overlay on a wavefunction are the total PES sampled at the energies in
play: the vibrational thresholds ``eps_v`` and/or ``eps_{v_init} + E`` over the
collision-energy range shown. `energy_contour_levels` assembles those, clipped
and de-cluttered. Pure array logic (model-independent) -- feed it
``ScatteringProblem.eps``.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["energy_contour_levels"]


def energy_contour_levels(
    *,
    eps: npt.ArrayLike | None = None,
    v_init: int = 0,
    energies: npt.ArrayLike | None = None,
    include_thresholds: bool = True,
    e_range: tuple[float, float] | None = None,
    max_levels: int = 12,
    min_spacing: float | None = None,
) -> list[float]:
    """Assemble energy-relevant potential-contour levels (Hartree), sorted.

    Parameters
    ----------
    eps : array_like, optional
        Vibrational energies. Their values become threshold levels (when
        `include_thresholds`), and ``eps[v_init]`` is the base for total energies.
    v_init : int, optional
        Initial vibrational level; ``eps[v_init] + E`` is the total energy.
    energies : array_like, optional
        Collision energies; each adds a total-energy level ``eps[v_init] + E``
        (requires `eps`).
    include_thresholds : bool, optional
        Include the vibrational-threshold levels ``eps_v``.
    e_range : tuple[float, float], optional
        Keep only levels within ``[lo, hi]`` (e.g. the shown potential's range).
    max_levels : int, optional
        Cap the number of levels, keeping an evenly-spread subset.
    min_spacing : float, optional
        Drop levels closer than this to the previous kept level (de-clutter).

    Returns
    -------
    list[float]
        Sorted, de-cluttered levels.
    """
    levels: set[float] = set()
    eps_arr = None if eps is None else np.atleast_1d(np.asarray(eps, dtype=np.float64))
    if include_thresholds and eps_arr is not None:
        levels |= {float(e) for e in eps_arr}
    if energies is not None:
        if eps_arr is None:
            raise ValueError("energies requires eps (to add E_tot = eps[v_init] + E)")
        base = float(eps_arr[v_init])
        levels |= {base + float(e) for e in np.atleast_1d(np.asarray(energies, dtype=np.float64))}

    out = sorted(levels)
    if e_range is not None:
        lo, hi = e_range
        out = [v for v in out if lo <= v <= hi]
    if min_spacing is not None and out:
        thinned = [out[0]]
        for v in out[1:]:
            if v - thinned[-1] >= min_spacing:
                thinned.append(v)
        out = thinned
    if max_levels is not None and len(out) > max_levels:
        idx = sorted(set(np.linspace(0, len(out) - 1, max_levels).round().astype(int)))
        out = [out[i] for i in idx]
    return out
