"""TD-matches-TI + convergence tests for the time-dependent VE cross section
(sub-project #4, Task 2 -- THE CRUX).

`docs/superpowers/specs/2026-07-22-n2-td-cross-section-design.md` ("Validation")
and the task brief: because
`S_TD(E) = (1/i)*integral_0^inf exp(i*E_tot*t)*<d_v'|exp(-i*H_res*t)|d_v> dt`
equals `S_TI(E) = <d_v'|(E_tot-H_res)^-1|d_v>` in the long-time limit, TD sigma
(`projects.n2_td_cross_section.td_cross_section.td_ve_cross_section`) is
checked against the already-validated TI oracle
(`projects.n2_ti_cross_section.cross_section.ve_cross_section`) -- an EXACT
differential oracle, not a loose cross-model comparison.

- **V1 (TD ~= TI):** at (E=0.1 Ha, v'=1) and (E=0.2 Ha, v'=2), sigma_TD at a
  tuned, converged (dt, n_steps) matches sigma_TI to `rtol <= 0.10`; sigma_TD
  is also real and >= 0.
- **V2 (convergence + depletion):** sigma_TD at (E=0.1, v'=1) agrees between
  `dt` and `dt/2` (same total propagation time `T`, so `n_steps` doubles) to
  `rtol <= 0.05`; and `||psi(T)|| < 0.1*||psi(0)||` (the resonance has
  depleted, so the finite-time energy transform is not truncated).

`dt`/`n_steps` tuning (see the development notes): the
resonance's own eigenmodes of `H_res` sit at `Re(E) ~ -0.7..-0.4 Ha`
(v_init=0's `eps[0] ~ -0.745` shifted by the ~2.3-2.5 eV Pi_g resonance), and
the Crank-Nicolson Cayley-transform phase error per step grows as
`~(E*dt)^3`, accumulating over `n_steps` steps -- so accuracy requires
`dt` small relative to `1/|E|`, not just `n_steps` large. `T = n_steps*dt =
1500` a.u. is long enough to deplete the resonance (`Gamma(R0) ~ 0.017 Ha`
gives a decay time ~`1/(Gamma/2) ~ 120` a.u.); `dt = 0.025` keeps the
per-step Cayley phase error for these modes under ~1e-4 rad/step.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from projects.n2_td_cross_section.td_cross_section import td_ve_cross_section
from projects.n2_ti_cross_section.cross_section import ve_cross_section
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states
from projects.n2_ti_cross_section.vres import vres_on_grid

_CONFIG = json.loads(
    (Path(__file__).resolve().parents[2] / "validation" / "n2" / "config.json").read_text()
)
MU = _CONFIG["reduced_mass"]  # N2 nuclear reduced mass (a.u.), 12766.36
N_VIB = 6  # v=0..5, matches the TI oracle's test setup

V_INIT = 0
# (collision energy E in Ha, final vibrational channel v') anchors from the
# task brief.
ANCHORS: list[tuple[float, int]] = [(0.1, 1), (0.2, 2)]

# Converged propagation: T = n_steps*dt = 1500 a.u. See module docstring and
# the development notes for the tuning study.
DT_CONVERGED = 0.025
N_STEPS_CONVERGED = 60000  # T = 1500 a.u.
# V2's "dt" reference point: same total T, half the time resolution.
DT_COARSE = 0.05
N_STEPS_COARSE = 30000  # T = 1500 a.u.


@pytest.fixture(scope="module")
def system():
    """Build the shared grid / vibrational states / V_d(R),Gamma(R) once.

    `vres_on_grid` walks a two-angle-matched electronic pole-finder
    continuation across ~300 nuclear grid points (~7s) -- module-scoped so
    the whole test file pays that cost exactly once.
    """
    grid = n2_nuclear_grid()
    eps, chi = vibrational_states(grid, MU, N_VIB)
    Vd, Gamma = vres_on_grid(grid)
    return grid, eps, chi, Vd, Gamma


@pytest.fixture(scope="module")
def td_converged(system):
    """sigma_TD at both anchors, computed once at the converged (dt, n_steps).

    Requests `E=[0.1, 0.2]` against `vprimes=[1, 2]` in a single propagation
    (the correlation functions `c_v'(t)` are E-independent, so this reuses
    one CN propagation for both anchors); returns `(sigma, norm_ratio)` with
    `sigma` shape `(2, 2)` -- anchor `i` is the diagonal element
    `sigma[i, i]` (E[i] paired with vprimes[i]).
    """
    grid, eps, chi, Vd, Gamma = system
    E = np.array([e for e, _ in ANCHORS])
    vprimes = [vp for _, vp in ANCHORS]
    sigma, norm_ratio = td_ve_cross_section(
        grid,
        MU,
        Vd,
        Gamma,
        eps,
        chi,
        V_INIT,
        vprimes,
        E,
        dt=DT_CONVERGED,
        n_steps=N_STEPS_CONVERGED,
        return_norm_ratio=True,
    )
    return sigma, norm_ratio


def test_v1_td_matches_ti_and_is_physical(system, td_converged):
    grid, eps, chi, Vd, Gamma = system
    sigma_td, _norm_ratio = td_converged

    print("\nTD vs TI at the V1 anchors:")
    for i, (e, vp) in enumerate(ANCHORS):
        sigma_td_i = sigma_td[i, i]
        sigma_ti_i = float(ve_cross_section(grid, MU, Vd, Gamma, eps, chi, V_INIT, [vp], e)[0])
        ratio = sigma_td_i / sigma_ti_i
        print(f"  E={e} Ha, v'={vp}: TD={sigma_td_i:.6e}  TI={sigma_ti_i:.6e}  ratio={ratio:.4f}")

        # Physical: non-negative.
        assert sigma_td_i.real >= -1e-12

        # TD ~= TI: the exact differential oracle, within the propagation's
        # finite-dt/finite-T discretization.
        assert sigma_td_i == pytest.approx(sigma_ti_i, rel=0.10)


def test_v2_convergence_in_dt_and_depletion(system, td_converged):
    grid, eps, chi, Vd, Gamma = system
    sigma_td, norm_ratio = td_converged
    e, vp = ANCHORS[0]
    sigma_fine = float(sigma_td[0, 0])  # dt=DT_CONVERGED result, already computed

    sigma_coarse, norm_ratio_coarse = td_ve_cross_section(
        grid,
        MU,
        Vd,
        Gamma,
        eps,
        chi,
        V_INIT,
        [vp],
        e,
        dt=DT_COARSE,
        n_steps=N_STEPS_COARSE,
        return_norm_ratio=True,
    )
    sigma_coarse_val = float(sigma_coarse[0])

    print(
        f"\nV2 convergence at (E={e}, v'={vp}): "
        f"dt={DT_COARSE} sigma={sigma_coarse_val:.6e} (norm_ratio={norm_ratio_coarse:.4e}); "
        f"dt/2={DT_CONVERGED} sigma={sigma_fine:.6e} (norm_ratio={norm_ratio:.4e})"
    )

    # Halving dt (same total T, so n_steps doubles) changes sigma negligibly.
    assert sigma_fine == pytest.approx(sigma_coarse_val, rel=0.05)

    # The resonance has depleted: the finite-time energy transform is not
    # truncated.
    assert norm_ratio < 0.1
    assert norm_ratio_coarse < 0.1
