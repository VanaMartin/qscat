"""TD-matches-TI + convergence tests for the exact 2-D VE cross section
(sub-project #7, Task 4 -- THE CRUX).

the development notes records the full tuning trail. Summary:

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
  Hankel half brought it to the right order of magnitude.
- Evolution operator: **order-3 diagonal Pade** (`qscat.evolution.make_pade_stepper`),
  `dt = 1.0`, `n_steps = 1500` (`T = 1500` a.u.) -- eMoScat's setting. Order-1
  Crank-Nicolson under-converges catastrophically over a multi-thousand-step
  run (~100% accumulated propagation error at dt=0.5-1.0, verified vs `expm`),
  which capped `sigma_TD/sigma_TI` at ~0.93/1.10 for `(E=0.10/0.15, v'=1)` and
  left the boomerang oscillations unresolved. The order-3 Pade operator
  (`O(dt^7)`) removes that: across 0.04-0.18 Ha `sigma_TD` matches the TI
  oracle to ~1-2% median for all channels (elastic + excitations), tracking
  the boomerang peaks point-by-point. `norm` still decays `~1.0 -> ~0.02` by
  `T=1500` (the resonance's formation and decay). See
  `docs/physics/n2-2d-td-cross-section.md`.

One ~5-min propagation is run ONCE at module scope (`td._propagate`, not the
public `td_ve_cross_section_2d`, so the SAME stored `c(t)` can be
transformed at multiple truncation lengths for V4 without re-propagating);
every test transforms it with the public `sigma_from_correlations` -- exactly
the function `td_ve_cross_section_2d` calls internally, so this exercises the
real transform, just not the outer scalar/array plumbing (covered separately,
cheaply, by `test_public_api_shape_contract`).
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

DT = 1.0
N_STEPS = 1500  # T = 1500 a.u., order-3 Pade (see module docstring)
N_STEPS_SHORT = 1000  # T = 1000 a.u., V4's shorter-truncation comparison point
PADE_ORDER = 3

# `ve_cross_section_2d` at the two anchors, computed once (cheap: one sparse
# LU solve per energy, not a propagation).
SIGMA_TI = {
    e: float(ve_cross_section_2d(TG, EPS, CHI, V_INIT, VPRIMES, e)[0]) for e in (0.10, 0.15)
}


@pytest.fixture(scope="module")
def propagation() -> td.PropagationResult:
    """The ONE ~5-min order-3 Pade propagation this file reuses (module docstring)."""
    return td._propagate(
        TG, EPS, CHI, V_INIT, VPRIMES,
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, order=PADE_ORDER,
    )


@pytest.mark.slow
def test_v2a_td_matches_ti_at_e010(propagation: td.PropagationResult) -> None:
    sigma_td = float(
        td.sigma_from_correlations(
            TG, propagation, EPS, V_INIT, VPRIMES, 0.10, dt=DT, wp_in=WP_IN, wp_out=WP_OUT
        )[0]
    )
    assert sigma_td >= 0.0
    # Order-3 Pade: measured ratio ~0.97 (was 0.931 with order-1 CN) -- the
    # rel=0.06 gate reflects the converged accuracy, not the old CN residual.
    assert sigma_td == pytest.approx(SIGMA_TI[0.10], rel=0.06)


@pytest.mark.slow
def test_v2a_td_matches_ti_at_e015(propagation: td.PropagationResult) -> None:
    """With the order-3 Pade operator the E=0.15 point is no longer a
    usable-window outlier (order-1 CN gave 1.103 here); measured ratio ~0.99.
    """
    sigma_td = float(
        td.sigma_from_correlations(
            TG, propagation, EPS, V_INIT, VPRIMES, 0.15, dt=DT, wp_in=WP_IN, wp_out=WP_OUT
        )[0]
    )
    assert sigma_td >= 0.0
    assert sigma_td == pytest.approx(SIGMA_TI[0.15], rel=0.06)


@pytest.mark.slow
def test_v2b_closed_channel_is_exactly_zero(propagation: td.PropagationResult) -> None:
    """`v'=1` is energetically closed once `E_tot - eps[1] <= 0`.

    `eps[1] - eps[0] = 0.0124` Ha; at `E = 0.0107 < eps[1] - eps[0]`, the
    channel is closed and `sigma_TD` must be exactly 0 (the `excess <= 0`
    branch in the underlying per-energy transform, not a small-but-nonzero
    residual).
    """
    e_closed = 0.0107
    assert e_closed + EPS[V_INIT] - EPS[VPRIMES[0]] <= 0.0  # confirms the channel is closed
    sigma_td = td.sigma_from_correlations(
        TG, propagation, EPS, V_INIT, VPRIMES, e_closed, dt=DT, wp_in=WP_IN, wp_out=WP_OUT
    )
    assert sigma_td[0] == 0.0


@pytest.mark.slow
def test_v4_finite_t_stability_and_depletion(propagation: td.PropagationResult) -> None:
    """Truncating the SAME stored `c(t)` at `T=1000` vs the full `T=1500`
    changes `sigma_TD(E=0.10)` by only a few percent, and `norm` has decayed
    well below the `< 0.05` depletion bar by `T=1500`. No second propagation
    needed: this truncates the one stored trajectory, which is exactly what
    makes this check "free" per the task brief. (Convergence with propagation
    length `T`; the order-3 Pade operator handles convergence with `dt`.)
    """
    full_sigma = float(
        td.sigma_from_correlations(
            TG, propagation, EPS, V_INIT, VPRIMES, 0.10, dt=DT, wp_in=WP_IN, wp_out=WP_OUT
        )[0]
    )
    short = SimpleNamespace(
        t=propagation.t[: N_STEPS_SHORT + 1], c=propagation.c[: N_STEPS_SHORT + 1, :]
    )
    short_sigma = float(
        td.sigma_from_correlations(
            TG, short, EPS, V_INIT, VPRIMES, 0.10, dt=DT, wp_in=WP_IN, wp_out=WP_OUT
        )[0]
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


# --- Elastic channel: the free-particle-reference subtraction (v' == v_init) ---
# The outgoing normalization factor C(E) multiplies every channel's S, so the
# inelastic sigma=|S|^2 already absorbs it -- but the diagonal sigma=|S-ref|^2
# only isolates genuine scattering if `ref` is the actual unscattered value
# S_free(E)=C(E), NOT a literal 1 (a free-particle V_int=0 propagation gives
# |S_elastic|~2*pi^2, not 1). Subtracting a V_int=0 reference restores it. See
# td_cross_section._sigma_one_energy and the `td-elastic-wavepacket-normalization`
# note. Gate energies sit in the RESOLVED part of the usable window (E>=0.13);
# below that the elastic channel degrades near threshold / on the steep resonance
# rise, the ordinary TD finite-T limit (module docstring), not this fix.
ELASTIC_GATE_ENERGIES = (0.14, 0.15)


@pytest.fixture(scope="module")
def elastic_propagations() -> tuple[td.PropagationResult, td.PropagationResult]:
    """The full and the V_int=0 free-reference elastic-channel propagations
    (~2x250s), reused across the elastic assertions -- one full run and one
    free (`_propagate(..., free=True)`) reference on the SAME wavepacket/grid.
    """
    full = td._propagate(
        TG, EPS, CHI, V_INIT, [V_INIT],
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, order=PADE_ORDER,
    )
    free = td._propagate(
        TG, EPS, CHI, V_INIT, [V_INIT],
        dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT, free=True, order=PADE_ORDER,
    )
    return full, free


@pytest.mark.slow
def test_v2c_td_elastic_matches_ti_with_free_reference(
    elastic_propagations: tuple[td.PropagationResult, td.PropagationResult],
) -> None:
    """Elastic sigma with the free-particle reference matches the exact TI oracle;
    the literal-1 subtraction is orders of magnitude too large. With the order-3
    Pade operator, measured elastic TD/TI at this config is ~1.01 (E=0.14),
    ~0.99 (E=0.15) -- rel=0.08 gates the converged accuracy, while the broken
    literal-1 result is >50x the TI value.
    """
    full, free = elastic_propagations
    for e in ELASTIC_GATE_ENERGIES:
        sigma_ti = float(ve_cross_section_2d(TG, EPS, CHI, V_INIT, [V_INIT], e)[0])
        sigma_fixed = float(
            td.sigma_from_correlations(
                TG, full, EPS, V_INIT, [V_INIT], e,
                dt=DT, wp_in=WP_IN, wp_out=WP_OUT, free_result=free,
            )[0]
        )
        sigma_literal1 = float(
            td.sigma_from_correlations(
                TG, full, EPS, V_INIT, [V_INIT], e, dt=DT, wp_in=WP_IN, wp_out=WP_OUT
            )[0]
        )
        assert sigma_fixed == pytest.approx(sigma_ti, rel=0.08)
        assert sigma_literal1 > 50.0 * sigma_ti  # the bug the free reference fixes
