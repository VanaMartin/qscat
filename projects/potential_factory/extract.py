"""Build a Target FROM an existing model with the repo's own forward models.

This is the round-trip data path: the extracted curves are exact properties of
the model, so fitting them back must recover the model."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.dissociation import anion_electronic_states
from qscat.core.nrm.coupling import gamma_from_coupling, v_dk_plus
from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
from qscat.model import ResonanceModel

from projects.potential_factory.target import (
    CouplingTarget,
    Curve,
    NeutralTarget,
    Provenance,
    ResonanceTarget,
    Target,
)
from projects.potential_factory.tracker import (
    DEFAULT_BOUND_WINDOW,
    DEFAULT_RES_WINDOW,
    MAX_STEP,
    ElectronicPair,
    Pole,
    Window,
)

__all__ = ["extract_target", "walk_t1"]


def walk_t1(
    model: ResonanceModel,
    pair: ElectronicPair,
    R: npt.NDArray[np.float64],
    *,
    seed_energy: complex | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """A per-node gated pole walk, rather than `resonance_pole_walk`
    (which freezes permanently on the first breakdown, measured).
    Each node is searched independently: recentre on the last
    FOUND pole (not necessarily the immediately preceding node -- a dropped
    node leaves the last found pole in place for the next node to recentre
    on), falling back to the two default windows if that finds nothing. A
    node with no gated pole anywhere is DROPPED, not frozen.

    `seed_energy` (electronic shift scale) recentres the FIRST node's search
    on a known pole -- the target's own -- so that a well holding more than
    one angle-stable state is walked on the intended one rather than on the
    default window's global residual-argmin.

    A fallback-window candidate is subject to the same continuity
    guard as the recentred window (only accepted if within `MAX_STEP` of the
    last found pole), since the fallback windows are wide and
    `find_resonance_pole` has no built-in preference for the tracked state.
    """
    last_pole: Pole | None = None
    R_ok: list[float] = []
    shift_ok: list[float] = []
    gamma_ok: list[float] = []
    for Rj in R:
        Rj = float(Rj)
        v0Rj = float(model.v0(np.asarray(Rj)).real)

        def v_fn(
            r: npt.NDArray[np.complex128], Rj: float = Rj, v0Rj: float = v0Rj
        ) -> npt.NDArray[np.complex128]:
            return np.asarray(model.surface(r, Rj) - v0Rj, dtype=np.complex128)

        pole: Pole | None = None
        if last_pole is None and seed_energy is not None and not np.isnan(seed_energy.real):
            E0 = complex(seed_energy)
            pole = pair.pole(
                v_fn,
                (E0.real - MAX_STEP, E0.real + MAX_STEP, min(E0.imag - MAX_STEP, -MAX_STEP), 1e-6),
            )
        if pole is None and last_pole is not None:
            E = last_pole.energy
            window: Window = (
                E.real - MAX_STEP,
                E.real + MAX_STEP,
                min(E.imag - MAX_STEP, -MAX_STEP),
                1e-6,
            )
            pole = pair.pole(v_fn, window)
        # The DEFAULT_* fallback windows are wide, and
        # `find_resonance_pole` returns the global residual-argmin over the
        # whole window with no preference for the state actually being
        # tracked (same reasoning as tracker.py's continuity guard
        # on Newton candidates) -- so once a previous pole exists, a
        # fallback candidate is accepted only if it is within one
        # recentred-window step of it; otherwise it is discarded as if no
        # pole had been found (the node may still be dropped, or a later
        # fallback may find a continuous one).
        if pole is None:
            cand = pair.pole(v_fn, DEFAULT_BOUND_WINDOW)
            if cand is not None and (
                last_pole is None or abs(cand.energy - last_pole.energy) <= MAX_STEP
            ):
                pole = cand
        if pole is None:
            cand = pair.pole(v_fn, DEFAULT_RES_WINDOW)
            if cand is not None and (
                last_pole is None or abs(cand.energy - last_pole.energy) <= MAX_STEP
            ):
                pole = cand
        if pole is not None:
            last_pole = pole
            R_ok.append(Rj)
            shift_ok.append(pole.shift)
            gamma_ok.append(pole.gamma)
    if len(R_ok) < 4:
        raise ValueError(
            f"resonance walk gated out too many nodes: only {len(R_ok)}/{R.size} survived"
        )
    return (
        np.asarray(R_ok, dtype=np.float64),
        np.asarray(shift_ok, dtype=np.float64),
        np.asarray(gamma_ok, dtype=np.float64),
    )


def extract_target(
    model: ResonanceModel,
    *,
    pair: ElectronicPair,
    R_desc: npt.ArrayLike,
    R_inf: float = 10.0,
    eps_window: tuple[float, float] = (0.002, 0.25),
    n_eps: int = 12,
    name: str = "model",
) -> Target:
    R = np.asarray(R_desc, dtype=np.float64)
    if R.size > 1 and np.any(np.diff(R) >= 0.0):
        raise ValueError("R_desc must be strictly descending")
    v0 = model.v0(R).real

    # T1: a per-node gated pole walk, not `resonance_pole_walk` (which
    # freezes permanently on the first breakdown, measured). Nodes with
    # no gated pole are dropped.
    R_ok, shift_ok, gamma_ok = walk_t1(model, pair, R)
    v0_ok = model.v0(R_ok).real
    v_ion_ok = v0_ok + shift_ok
    eps_inf, _ = anion_electronic_states(pair.grid_a, model, R_inf, 1)
    ea = -(eps_inf[0] - model.v0(R_inf).real)

    # T3: the model's own Gamma~(eps, R) with the R-independent discrete state.
    # Uses ALL of R, not just the T1 survivors.
    phi_d = AsymptoticDiscreteState(pair.grid_a, model, R_inf)
    eps = np.geomspace(eps_window[0], eps_window[1], n_eps)
    table = np.empty((eps.size, R.size))
    R_asc = R[::-1]
    for i, e in enumerate(eps):
        table[i] = gamma_from_coupling(v_dk_plus(pair.grid_a, model, phi_d, R_asc, float(e)))

    rng = (float(R.min()), float(R.max()))
    # ResonanceTarget's R_range is the range of the SURVIVING T1
    # nodes the v_ion/gamma curves are actually built from, not the full
    # requested R -- NeutralTarget keeps the full range since v0 is
    # evaluated on every node, T1 dropouts notwithstanding.
    rng_ok = (float(R_ok.min()), float(R_ok.max()))
    return Target(
        name=name,
        mu=model.mu,
        ell=model.ell,
        charge=model.charge,
        coordinates=("R",),
        neutral=NeutralTarget(Curve.from_table(R, v0), {}, rng),
        resonance=ResonanceTarget(
            Curve.from_table(R_ok, v_ion_ok),
            Curve.from_table(R_ok, gamma_ok),
            float(ea),
            rng_ok,
            R_inf=float(R_inf),
        ),
        coupling=CouplingTarget.from_table(eps, R_asc, table, alpha=model.ell + 0.5),
        provenance={"all": Provenance(f"extract_target({name})", "computed, not published")},
    )
