"""Model-independent local complex potential `(V_d(R), Gamma(R))`.

Sub-project B: the local-complex-potential (LCP) approximation of
dissociative attachment (and vibrational excitation) reduces the full 2-D
(electronic r x nuclear R) resonance problem to a 1-D nuclear problem by
replacing the fixed-R electronic resonance with a single complex number at
each R -- eMoScat's `ModelLCP` (`v0(R) + E_res(R)` real part, width
`-2*Im(E_res(R))`). This is the RESEARCH-PROGRAM "approximation under test":
the exact 2-D solver (`qscat.core.dissociation`/`driven`) is the oracle, and
`local_complex_potential` is the reduction whose accuracy against that oracle
is what sub-project B is built to measure -- not a description of the "real"
physics.

`V_d(R) = Re(E_pole(R))`, `Gamma(R) = max(0, -2*Im(E_pole(R)))`, where
`E_pole(R)` is the two-ECS-angle-matched resonance pole (`qscat.ecs.
find_resonance_pole`) of the fixed-R electronic Hamiltonian
`H_el(R) = -1/2 d^2/dr^2 + model.surface(r, R)`. Because `model.surface`
ALREADY INCLUDES `v0(R)` (the neutral-molecule channel potential), `V_d =
Re(E_pole)` directly -- adding `v0(R)` again would double-count it. (Contrast
`projects/n2_ti_cross_section/vres.py`, whose `v_eff_el` EXCLUDES `v0`, so
that code adds `v0(R) + E_res(R)` separately; the two bookkeeping schemes
must and do agree on the observable `V_d(R)`, checked in
`test_lcp.py::test_matches_n2_vres_oracle`.)

The continuation is seeded from the bound anion state at `R_inf =
nuclear_grid.R0` (`qscat.core.dissociation.anion_electronic_states`, the
same electronic bound state the exact-2D DA solver uses for its exit
channel) and walked INWARD (decreasing R), matching
`projects/n2_resonance/pole.py`'s `resonance_curve` /
`projects/n2_ti_cross_section/vres.py`'s `vres_on_grid` continuation, made
model-independent. The first rejected step (pole finder raises, or residual
>= `resid_tol`) freezes the electronic shift `s = V_d - v0(R)` at its last
accepted value for all remaining smaller-R points (`vres.py`'s documented
small-R breakdown -- physically irrelevant there since `v0(R)` is already
Ha above threshold and the vibrational wavefunction has negligible density).
The complex ECS tail (R.imag != 0) clamps `Gamma = 0` and freezes the shift
at its outermost-real-R (asymptotic) value.

See `docs/superpowers/specs/2026-07-27-da-cross-sections-design.md`.

`qscat.core` never imports `qscat.model`/`projects.*` at runtime: `model` is
typed against the `ResonanceModel` protocol under `TYPE_CHECKING` only,
exactly like `driven.py`/`dissociation.py`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt
import scipy.sparse as sp

from qscat.dvr import FemDvrEcsGrid, eigen, kinetic, kinetic_sparse
from qscat.ecs import find_resonance_pole
from qscat.evolution import make_pade_stepper

from .dissociation import anion_electronic_states

if TYPE_CHECKING:
    from qscat.model import ResonanceModel

__all__ = ["local_complex_potential", "lcp_da_cross_section"]


def _h_el(model: ResonanceModel, R: complex, g: FemDvrEcsGrid) -> npt.NDArray[np.complex128]:
    return kinetic(g, 1.0) + np.diag(model.surface(g.points, R))


def local_complex_potential(
    model: ResonanceModel,
    nuclear_grid: FemDvrEcsGrid,
    elec_grid_a: FemDvrEcsGrid,
    elec_grid_b: FemDvrEcsGrid,
    *,
    re_half_width: float = 0.05,
    im_half_width: float = 0.05,
    resid_tol: float = 1e-3,
) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.float64]]:
    """Local complex potential `(V_d(R), Gamma(R))` for the LCP DA/VE model.

    `V_d(R) = Re(E_pole(R))`, `Gamma(R) = max(0, -2 Im(E_pole(R)))`, `E_pole`
    the two-angle-matched resonance pole of `-1/2 d^2/dr^2 + model.surface(r,R)`
    (surface includes v0(R), so V_d = Re(E_pole) directly). Seeded from the
    bound anion at `R_inf = nuclear_grid.R0` (`anion_electronic_states`) and
    continued INWARD; small-R breakdown freezes the electronic shift; the
    complex tail clamps Gamma=0. See module docstring.
    """
    R_inf = nuclear_grid.R0
    eps_e, _ = anion_electronic_states(elec_grid_a, model, R_inf, 1)
    window = (eps_e[0] - re_half_width, eps_e[0] + re_half_width, -im_half_width, im_half_width)

    pts = nuclear_grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    order = np.argsort(pts[real_idx].real)[::-1]          # descending R: outer -> inner
    walk = real_idx[order]
    R_real = pts[walk].real

    shift = np.empty(walk.size, dtype=np.float64)          # s = V_d - v0(R)
    gamma_w = np.empty(walk.size, dtype=np.float64)
    last_s: float | None = None
    last_g = 0.0
    broken = False
    for j in range(walk.size):
        R = float(R_real[j])
        if not broken:
            try:
                E_pole, resid = find_resonance_pole(
                    eigen(_h_el(model, R, elec_grid_a))[0],
                    eigen(_h_el(model, R, elec_grid_b))[0],
                    window,
                )
            except (ValueError, np.linalg.LinAlgError):
                resid = np.inf
            else:
                if resid < resid_tol:
                    v0R = float(np.real(model.v0(np.asarray(R))))
                    last_s = E_pole.real - v0R
                    last_g = max(0.0, -2.0 * E_pole.imag)
                    window = (E_pole.real - re_half_width, E_pole.real + re_half_width,
                              E_pole.imag - im_half_width, E_pole.imag + im_half_width)
                    shift[j], gamma_w[j] = last_s, last_g
                    continue
            broken = True
        if last_s is None:
            raise RuntimeError("local_complex_potential: pole finder failed at the seed edge")
        shift[j], gamma_w[j] = last_s, last_g

    Vd = np.empty(nuclear_grid.n, dtype=np.complex128)
    Gamma = np.zeros(nuclear_grid.n, dtype=np.float64)
    Vd[walk] = model.v0(R_real) + shift
    Gamma[walk] = gamma_w

    tail = np.flatnonzero(pts.imag != 0.0)
    if tail.size:
        assert last_s is not None
        s_asym = shift[0]                                  # shift at the largest real R
        Vd[tail] = model.v0(pts[tail]) + s_asym
    return Vd, Gamma


def _quadrature_weights(n_t: int) -> npt.NDArray[np.float64]:
    """Composite Simpson weights (unscaled by dt); trapezoidal fallback for even n_t.
    Copied from projects/n2_td_cross_section/td_cross_section.py."""
    if n_t < 2:
        raise ValueError("need at least 2 time samples")
    if n_t % 2 == 1:
        w = np.ones(n_t)
        w[1:-1:2] = 4.0
        w[2:-1:2] = 2.0
        w /= 3.0
    else:
        w = np.full(n_t, 2.0)
        w[0] = w[-1] = 1.0
        w /= 2.0
    return np.asarray(w, dtype=np.float64)


def lcp_da_cross_section(
    nuclear_grid: FemDvrEcsGrid,
    mu: float,
    Vd: npt.NDArray[np.complex128],
    Gamma: npt.NDArray[np.float64],
    eps: npt.NDArray[np.float64],
    chi: npt.NDArray[np.complex128],
    v_init: int,
    eps_e: float,
    E: float | npt.ArrayLike,
    *,
    dt: float = 1.0,
    n_steps: int = 1500,
    order: int = 3,
) -> npt.NDArray[np.float64]:
    """LCP dissociative-attachment sigma_DA(E) (bohr^2): the 1-D nuclear
    doorway propagated on V_d - i Gamma/2, DA = flux at the dissociation
    boundary. eMoScat ModelLCP/SMatrix.cpp. The approximation under test vs the
    exact-2D `da_cross_section` oracle."""
    pts = nuclear_grid.points
    real_idx = np.flatnonzero(pts.imag == 0.0)
    b = int(real_idx[np.argmax(pts[real_idx].real)])   # outermost real point index
    X = float(pts[b].real)

    doorway = np.sqrt(Gamma / (2.0 * np.pi)).astype(np.complex128) * chi[v_init]
    H_res = kinetic_sparse(nuclear_grid, mu) + sp.diags(Vd - 0.5j * Gamma)
    step = make_pade_stepper(H_res.tocsc(), dt, order)

    n_t = n_steps + 1
    t = np.arange(n_t, dtype=np.float64) * dt
    psi_X = np.empty(n_t, dtype=np.complex128)
    psi = doorway.copy()
    for n in range(n_t):
        psi_X[n] = psi[b]
        if n < n_steps:
            psi = step(psi)

    w = _quadrature_weights(n_t)
    e_arr = np.atleast_1d(np.asarray(E, dtype=np.float64))
    out = np.zeros(e_arr.size, dtype=np.float64)
    for ie, e in enumerate(e_arr):
        if e <= 0.0:
            continue
        e_tot = float(e) + eps[v_init]
        e_dr = e_tot - eps_e
        if e_dr <= 0.0:
            continue
        K = float(np.sqrt(2.0 * mu * e_dr))
        phase = np.exp(1j * e_tot * t)
        S = np.sqrt(K / (2.0 * np.pi * mu)) * np.exp(-1j * K * X) * np.sum(w * phase * psi_X) * dt
        out[ie] = 4.0 * np.pi**3 * abs(S) ** 2 / (2.0 * float(e))

    scalar = np.isscalar(E) or (isinstance(E, np.ndarray) and np.ndim(E) == 0)
    return np.asarray(out[0] if scalar else out, dtype=np.float64)
