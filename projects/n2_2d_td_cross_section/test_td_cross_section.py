"""TD-matches-TI + convergence tests for the exact 2-D VE cross section
(sub-project #7, Task 4 -- THE CRUX).

`.superpowers/sdd/task-4-report.md` records the full tuning trail. Summary:

- Grid: electronic `n2_electronic_grid(r_max=50, order=8, n_complex=6)`
  (`R0 = 50`), nuclear `n2_nuclear_grid(quadrature=10, r_max=22, n_complex=5)`
  (matches `test_td_propagation.py`/`test_wavepacket.py`'s test-scale grid,
  widened only in `r_max` so the outgoing test function fits).
- Incident wavepacket `wp_in = {r0: 25, p0: -0.5, sigma: 5.0}`: `p0` centers
  the spectral peak at `p0**2/2 = 0.125` Ha, between the two anchor
  energies, so `|eta_in(E)|` is large and balanced at both (measured
  `|eta_in(0.10)| = 2.70`, `|eta_in(0.15)| = 2.64` -- see the report).
- Outgoing test function `wp_out = {r0_out: 35, p0_out: 0.5, sigma_out: 4.0}`
  (scaled down from eMoScat's production r0=75 to fit this box; well
  outside the interaction range, which vanishes by r~5-6 since
  `alpha_c = 0.4`).
- `F_out` in `eta_outgoing` MUST be the outgoing Hankel half
  (`correlation.py`'s `_outgoing_coeffs`, `h^{(1)}_{E,l}/2`), not the
  regular function: using the regular function (debug order item 7) gave
  `sigma_TD` five to six orders of magnitude too small (ratio ~1e-5); the
  Hankel half brought it to within ~10-15% of the TI oracle.
- Converged propagation: `dt = 0.5`, `n_steps = 3000` (`T = 1500` a.u.);
  `norm` decays `1.0 -> 0.024` (fully depleted, well under the `< 0.05`
  bar). Measured at this config: `sigma_TD/sigma_TI = 0.931` at
  `(E=0.10, v'=1)` and `1.103` at `(E=0.15, v'=1)`. The `E=0.15` residual is
  larger because it sits farther from the incident wavepacket's spectral
  peak (`0.125` Ha) -- the deconvolution `eta_in(E)` amplifies whatever
  residual truncation/discretization error remains, more so away from the
  peak. This is the documented "usable spectral window" effect anticipated
  by the task brief, not a solver defect -- confirmed by scanning multiple
  propagation lengths (`T` in `[600, 900, ..., 1800]`, `.superpowers/sdd/
  task-4-report.md`): the `E=0.15` ratio oscillates in `[1.10, 1.15]` across
  that whole range, never approaching 1 as tightly as `E=0.10` does, which
  is the signature of a usable-window residual rather than an unconverged
  transient.

One ~250s propagation is run ONCE at module scope (`_propagate`, not the
public `td_ve_cross_section_2d`, so the SAME stored `c(t)` can be
transformed at multiple truncation lengths for V4 without re-propagating);
every test transforms it with `_sigma_from_correlations` -- exactly the
function `td_ve_cross_section_2d` calls internally per-energy, so this
exercises the real transform, just not the outer scalar/array plumbing
(covered separately, cheaply, by `test_public_api_shape_contract`).
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from qscat.dvr import TensorGrid

from projects.n2_2d_cross_section.cross_section_2d import ve_cross_section_2d
from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_2d_cross_section.hamiltonian2d import MU
from projects.n2_2d_td_cross_section import td_cross_section as td
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

TG = TensorGrid(
    [
        n2_electronic_grid(r_max=50.0, order=8, n_complex=6),
        n2_nuclear_grid(quadrature=10, r_max=22.0, n_complex=5),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], MU, 4)

V_INIT = 0
VPRIMES = [1]
WP_IN = {"r0": 25.0, "p0": -0.5, "sigma": 5.0}
WP_OUT = {"r0_out": 35.0, "p0_out": 0.5, "sigma_out": 4.0}

DT = 0.5
N_STEPS = 3000  # T = 1500 a.u., converged (see module docstring)
N_STEPS_SHORT = 2000  # T = 1000 a.u., V4's shorter-truncation comparison point

# `ve_cross_section_2d` at the two anchors, computed once (cheap: one sparse
# LU solve per energy, not a propagation).
SIGMA_TI = {
    e: float(ve_cross_section_2d(TG, EPS, CHI, V_INIT, VPRIMES, e)[0]) for e in (0.10, 0.15)
}


@pytest.fixture(scope="module")
def propagation() -> td.PropagationResult:
    """The ONE ~250s propagation this whole file reuses (see module docstring)."""
    return td._propagate(
        TG, EPS, CHI, V_INIT, VPRIMES, dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT
    )


@pytest.mark.slow
def test_v2a_td_matches_ti_at_e010(propagation: td.PropagationResult) -> None:
    sigma_td = float(
        td._sigma_from_correlations(TG, propagation, EPS, V_INIT, VPRIMES, 0.10, DT, WP_IN, WP_OUT)[
            0
        ]
    )
    assert sigma_td >= 0.0
    # Measured ratio 0.931 at the converged (dt, n_steps) -- see module docstring.
    assert sigma_td == pytest.approx(SIGMA_TI[0.10], rel=0.10)


@pytest.mark.slow
def test_v2a_td_matches_ti_at_e015_usable_window_edge(propagation: td.PropagationResult) -> None:
    """E=0.15 sits farther from the incident wavepacket's spectral peak
    (`p0**2/2 = 0.125` Ha) than E=0.10 does, so `|eta_in(E)|` is smaller and
    the deconvolution is noisier here -- the documented usable-window effect
    (module docstring), not a solver defect. Measured ratio 1.103; the
    tolerance is set just above the measured value, not tightened to it.
    """
    sigma_td = float(
        td._sigma_from_correlations(TG, propagation, EPS, V_INIT, VPRIMES, 0.15, DT, WP_IN, WP_OUT)[
            0
        ]
    )
    assert sigma_td >= 0.0
    assert sigma_td == pytest.approx(SIGMA_TI[0.15], rel=0.15)


@pytest.mark.slow
def test_v2b_closed_channel_is_exactly_zero(propagation: td.PropagationResult) -> None:
    """`v'=1` is energetically closed once `E_tot - eps[1] <= 0`.

    `eps[1] - eps[0] = 0.0124` Ha; at `E = 0.0107 < eps[1] - eps[0]`, the
    channel is closed and `sigma_TD` must be exactly 0 (the `excess <= 0`
    branch in `_sigma_from_correlations`, not a small-but-nonzero residual).
    """
    e_closed = 0.0107
    assert e_closed + EPS[V_INIT] - EPS[VPRIMES[0]] <= 0.0  # confirms the channel is closed
    sigma_td = td._sigma_from_correlations(
        TG, propagation, EPS, V_INIT, VPRIMES, e_closed, DT, WP_IN, WP_OUT
    )
    assert sigma_td[0] == 0.0


@pytest.mark.slow
def test_v4_finite_t_stability_and_depletion(propagation: td.PropagationResult) -> None:
    """Truncating the SAME stored `c(t)` at `T=1000` vs the full `T=1500`
    changes `sigma_TD(E=0.10)` by only a few percent (measured 0.958 vs
    0.931 -- a 2.8% relative change), and `norm` has decayed well below the
    `< 0.05` depletion bar by `T=1500`. No second propagation needed: this
    truncates the one stored trajectory, which is exactly what makes this
    check "free" per the task brief.
    """
    full_sigma = float(
        td._sigma_from_correlations(TG, propagation, EPS, V_INIT, VPRIMES, 0.10, DT, WP_IN, WP_OUT)[
            0
        ]
    )
    short = SimpleNamespace(
        t=propagation.t[: N_STEPS_SHORT + 1], c=propagation.c[: N_STEPS_SHORT + 1, :]
    )
    short_sigma = float(
        td._sigma_from_correlations(TG, short, EPS, V_INIT, VPRIMES, 0.10, DT, WP_IN, WP_OUT)[0]
    )
    assert short_sigma == pytest.approx(full_sigma, rel=0.10)

    assert propagation.norm[0] == pytest.approx(1.0, abs=1e-9)
    assert propagation.norm[-1] < 0.05  # depleted -- see module docstring's norm profile


def test_public_api_shape_contract() -> None:
    """Cheap (not `slow`): `td_ve_cross_section_2d`'s scalar/array `E`
    return-shape convention, matching `ve_cross_section_2d`'s -- exercised
    with a throwaway few-step propagation (not a physically converged one;
    only the outer plumbing is under test here, the transform itself is
    covered above).
    """
    sigma_scalar = td.td_ve_cross_section_2d(
        TG, EPS, CHI, V_INIT, VPRIMES, 0.10, dt=0.1, n_steps=2, wp_in=WP_IN, wp_out=WP_OUT
    )
    assert sigma_scalar.shape == (len(VPRIMES),)

    sigma_array = td.td_ve_cross_section_2d(
        TG, EPS, CHI, V_INIT, VPRIMES, [0.10, 0.15], dt=0.1, n_steps=2, wp_in=WP_IN, wp_out=WP_OUT
    )
    assert sigma_array.shape == (2, len(VPRIMES))
