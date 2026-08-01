"""H2+ TIME-DEPENDENT dissociative-recombination (TD-DR) experiment.

The time-dependent counterpart of `validation.h2plus.dr` (which runs the exact
TI `qscat.core.dr_cross_section` oracle). TD-DR propagates an incident electron
wavepacket under `H_2D` and reads sigma_DR(E) off the outgoing NUCLEAR flux with
the SP2 nuclear-axis extractors (`Flux`/`Dirac`/`TannorWeeks`,
`qscat.core.time_dependent.td_da_cross_sections_all`).

There is nothing DA-specific about that machinery for a CHARGED target: it is
the SAME nuclear-dissociation extraction, specialized to the ion by
`H2P.charge == -1` (so the outgoing wave is the Coulomb `coulomb_h1_en`, not a
free Hankel) and by `n_channels` (the exit states are the Rydberg series via
`anion_electronic_states`, exactly as the TI `dr_cross_section` uses them). So
TD-DR = `td_da_cross_sections_all(..., H2P, n_channels=N_CHANNELS)`.

## eMoScat parametrisation (input/experimental/H2p.json, transcribed here)

- Grid: `config.full_grid()` -- electronic real->1300 bohr + 5-deg ECS tail,
  nuclear real->14 bohr + 22-deg ECS tail, order 8 (the ~1.15M-unknown deck ->
  Docker/MUMPS; `config.proxy_grid()`/a small grid for a laptop).
- Incident: electronic Gaussian at `r0=800`, `p0=-0.25`, `sigma=8.0`, channel 0
  (the H2+ ground vibrational level). On a reduced grid `r0` is scaled inside
  the electronic box (an off-box incident lands in the ECS tail -> divergence;
  see docs/physics/td-da.md's KEY LESSON).
- DR extraction: nuclear, outgoing, method "all" (flux+delta+tw), 3 channels,
  nuclear test packet `R=12`, `K=15`, `sigma=0.4`.
- Evolution: `dt=10.0`, order-3 (diagonal Pade == eMoScat "Crank-Nicolson
  order 3"). eMoScat propagates until `||psi|| < 1e-20`; here `n_steps` is an
  explicit argument (the norm-cutoff loop is a follow-on).

## Running in Docker (MUMPS)

The full deck needs the MUMPS backend (`qscat[mumps]`, provisioned in the
docker `test` stage). Build + run:

    docker/build.sh test
    docker run --rm <image> uv run python -m validation.h2plus.td_dr

`SparseLU`/the sparse Pade stepper auto-select MUMPS when present; on a laptop
without MUMPS use a small grid (SuperLU) via `compute_td_dr(small_grid, ...)`.

NOTE (SP2 tech-debt, relevant here): the FLUX method's charged-Coulomb outgoing
DERIVATIVE (`outgoing_surface_wave`'s `dphi_out` at `charge != 0`) is a finite
difference that was never exercised for a neutral target -- H2+ is the first
real charged use. `delta`/`tw` (which need only the Coulomb VALUE, not its
derivative) are the safer reads; treat a flux/delta disagreement as a signal to
scrutinize that derivative, not as physics.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from qscat.core.time_dependent import td_da_cross_sections_all
from qscat.core.vibrational import vibrational_states
from qscat.dvr import TensorGrid
from qscat.model import H2P

from validation.h2plus.config import N_CHANNELS, energy_grid, full_grid

__all__ = ["compute_td_dr", "WP_OUT_DR", "DT", "ORDER"]

# eMoScat H2p.json evolution + DR test-function parametrisation.
DT = 10.0
ORDER = 3
P0_INCIDENT = -0.25
SIGMA_INCIDENT = 8.0
R0_INCIDENT_FULL = 800.0  # electronic launch position on the full (1300-bohr) deck
# nuclear DR outgoing test packet (Tannor-Weeks): position R=12, impulse 15, thickness 0.4
WP_OUT_DR = {"r0_out": 12.0, "p0_out": 15.0, "sigma_out": 0.4}
SURFACE_R = 12.0  # flux surface / delta point (nuclear), at the eMoScat test-function position


def _nuclear_index_near(grid: TensorGrid, r_value: float) -> int:
    nuc = grid.grids[1]
    rp = nuc.real_points
    masked = np.where(rp <= nuc.R0, rp, np.inf)
    return int(np.argmin(np.abs(masked - r_value)))


def compute_td_dr(
    grid: TensorGrid,
    *,
    energies: npt.ArrayLike,
    n_steps: int,
    r0_incident: float = R0_INCIDENT_FULL,
    n_channels: int = N_CHANNELS,
    dt: float = DT,
    order: int = ORDER,
) -> tuple[npt.NDArray[np.float64], dict[str, npt.NDArray[np.float64]]]:
    """sigma_DR(E) for H2+ on `grid`, via the three TD nuclear extractors from
    ONE shared propagation.

    Diagonalizes the H2+ ION vibrational states on `H2P.v0` (`v_init == 0`),
    builds the incident electron wavepacket (`r0_incident`/`P0_INCIDENT`/
    `SIGMA_INCIDENT`), propagates `n_steps` of order-`order` Pade at `dt`, and
    reads sigma_DR off the outgoing nuclear flux/delta/tw against the Rydberg
    exit channels (`n_channels`). Returns `(energies, {"flow","delta","tw"})`,
    each sigma shape `(len(energies), n_channels)`.
    """
    eps, chi = vibrational_states(grid.grids[1], H2P.mu, max(3, n_channels + 1), H2P.v0)
    surface = _nuclear_index_near(grid, SURFACE_R)
    e_arr = np.asarray(energies, dtype=np.float64)
    sigmas = td_da_cross_sections_all(
        grid,
        H2P,
        eps,
        chi,
        0,
        e_arr,
        dt=dt,
        n_steps=n_steps,
        wp_in={"r0": r0_incident, "p0": P0_INCIDENT, "sigma": SIGMA_INCIDENT},
        surface=surface,
        position=surface,
        wp_out=WP_OUT_DR,
        n_channels=n_channels,
        order=order,
    )
    return e_arr, sigmas


def main() -> None:
    """Docker/MUMPS smoke run: the full 1300-bohr H2+ deck, the eMoScat incident
    (r0=800), a couple of energies, a modest `n_steps`. NOT run in the test
    suite (~1.15M unknowns; MUMPS required). Full convergence needs many more
    steps (eMoScat propagates until ||psi|| < 1e-20) -- this is a "does it run +
    produce finite sigma_DR" smoke, comparable to `validation.h2plus.dr.main`.
    """
    grid = full_grid()
    energies = energy_grid()[:2]
    n_steps = 2000
    print(
        f"H2+ TD-DR smoke run: grid.size={grid.size}, energies={energies}, "
        f"n_steps={n_steps}, dt={DT}, r0={R0_INCIDENT_FULL}",
        flush=True,
    )
    e_out, sigmas = compute_td_dr(grid, energies=energies, n_steps=n_steps)
    out_path = "h2plus_td_dr_smoke.npz"
    np.savez(
        out_path,
        energies=e_out,
        sigma_flow=sigmas["flow"],
        sigma_delta=sigmas["delta"],
        sigma_tw=sigmas["tw"],
    )
    for method, sigma in sigmas.items():
        print(f"sigma_DR[{method}] shape={sigma.shape}:\n{sigma}", flush=True)
    print(f"saved to {out_path}", flush=True)


if __name__ == "__main__":
    main()
