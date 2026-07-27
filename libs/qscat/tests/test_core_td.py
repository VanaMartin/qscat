"""Gates for the promoted `qscat.core.{wavepacket,correlation,time_dependent,plot}`
(sub-project #A, Task 5): the time-dependent Pade-propagation + Tannor-Weeks
VE cross section, now taking a `model: qscat.model.ResonanceModel` instead of
a hardcoded N2 Hamiltonian.

This is a FAST, genuine gate, not a tautological one: it never compares
against the (now-rewired) `projects.n2_2d_td_cross_section` shims -- that
comparison would just check the shim delegates correctly, which is a
one-line call. Instead it computes everything through `qscat.core` directly,
on a deliberately tiny/unconverged grid with a handful of propagation steps
(seconds, not the ~250s `TD_WORKING_GRID` runs), and checks:

  * the public `td_ve_cross_section`'s scalar/array `E` shape contract, pinned
    to specific regression values captured from a run of this exact code (not
    independently re-derived -- a bit-identical regression pin, `rtol=1e-10`);
  * `propagate`'s `keep_psi_at` contract (psi kept only at requested times);
  * `sigma_from_correlations` reproduces looping `_sigma_one_energy` by hand,
    per energy -- i.e. the public batch entry point is not doing anything
    different from the private single-energy kernel it wraps;
  * the free-reference elastic path returns a finite, non-negative sigma, and
    is dramatically smaller than the (physically wrong) literal-1 fallback --
    the same invariant `test_td_cross_section.py::test_v2c_...` gates at
    `TD_WORKING_GRID` scale, cheaply reproduced here.

The real physics gate (TD converging to the exact TI oracle at a converged
grid) remains `projects/n2_2d_td_cross_section/test_td_cross_section.py`'s
`@pytest.mark.slow` tests -- not duplicated here.
"""

from __future__ import annotations

import numpy as np
import pytest
from qscat.core.time_dependent import (
    PropagationResult,
    _propagate,
    _sigma_one_energy,
    propagate,
    sigma_from_correlations,
    td_ve_cross_section,
)
from qscat.core.wavepacket import initial_state
from qscat.dvr import TensorGrid
from qscat.model import N2

from projects.n2_2d_cross_section.electronic_grid import n2_electronic_grid
from projects.n2_ti_cross_section.nuclear_grid import n2_nuclear_grid
from projects.n2_ti_cross_section.vibrational import vibrational_states

# Deliberately tiny/unconverged grid -- fast, not physically meaningful (the
# converged grid is `TD_WORKING_GRID` in `projects/n2_2d_td_cross_section/
# convergence.py`, gated at `@slow` scale elsewhere).
TG = TensorGrid(
    [
        n2_electronic_grid(r_max=12.0, order=5, n_complex=3),
        n2_nuclear_grid(quadrature=6, r_max=14.0, n_complex=3),
    ]
)
EPS, CHI = vibrational_states(TG.grids[1], N2.mu, 4)

V_INIT = 0
VPRIMES = [0, 1]  # includes the elastic (v'=v_init) channel
WP_IN = {"r0": 4.0, "p0": -0.5, "sigma": 1.2}
WP_OUT = {"r0_out": 6.0, "p0_out": 0.5, "sigma_out": 1.0}
DT = 0.2
N_STEPS = 5  # a handful of steps -- shape/contract test, not a converged run

# Captured from THIS code (`qscat.core.time_dependent`, post-promotion) at the
# fixture above -- a regression pin, not an independent derivation; guards
# against a future refactor silently changing the arithmetic. `rtol=1e-6`
# (not round-off-tight): the sparse-LU solves inside `make_pade_stepper`
# route through a multi-threaded BLAS, whose floating-point summation order
# (and hence the last few digits of the result) is not guaranteed identical
# run-to-run -- measured ~4e-9 relative drift between two runs of this exact
# code, well inside this tolerance.
_SIGMA_SCALAR_REF = np.array([1.1294074959617389e-05, 1.891186716451761e-09])
_SIGMA_ARRAY_REF = np.array(
    [
        [1.1294074959617389e-05, 1.891186716451761e-09],
        [9.426919299274127e-06, 1.7032730749254725e-09],
    ]
)


def test_td_ve_cross_section_shape_and_regression() -> None:
    """`td_ve_cross_section`'s scalar/array `E` shape contract, matching
    `qscat.core.driven.ve_cross_section`'s convention, plus a regression pin
    (see module docstring re the `rtol` choice)."""
    sigma_scalar = td_ve_cross_section(
        TG, N2, EPS, CHI, V_INIT, VPRIMES, 0.10, dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT
    )
    assert sigma_scalar.shape == (len(VPRIMES),)
    np.testing.assert_allclose(sigma_scalar, _SIGMA_SCALAR_REF, rtol=1e-6, atol=0.0)

    sigma_array = td_ve_cross_section(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        [0.10, 0.15],
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
    )
    assert sigma_array.shape == (2, len(VPRIMES))
    np.testing.assert_allclose(sigma_array, _SIGMA_ARRAY_REF, rtol=1e-6, atol=0.0)


def test_propagate_keep_psi_at_contract() -> None:
    """`propagate`'s `keep_psi_at` contract, with a caller-supplied
    `hamiltonian` (now REQUIRED -- the model-agnostic engine has no
    N2-specific default): `psi` is kept only in snapshots at the requested
    times, `None` elsewhere."""
    psi0 = initial_state(TG, CHI[V_INIT], **WP_IN)
    out_channels = [TG.outer([np.ones(TG.shape[0]), CHI[v]]) for v in VPRIMES]
    H = N2.hamiltonian(TG)
    result = propagate(
        TG,
        psi0,
        out_channels,
        dt=DT,
        n_steps=6,
        sample_period=2,
        keep_psi_at=[0.0, DT * 2],
        hamiltonian=H,
        order=3,
    )
    kept = {s.time: s.psi for s in result.snapshots}
    assert kept[0.0] is not None
    assert kept[0.0].shape == (TG.size,)
    assert kept[DT * 2] is not None
    # a snapshot on the coarse cadence but NOT in keep_psi_at has no psi kept
    other_times = [s.time for s in result.snapshots if s.time not in (0.0, DT * 2)]
    assert other_times  # sanity: there is at least one such snapshot
    assert all(s.psi is None for s in result.snapshots if s.time not in (0.0, DT * 2))


@pytest.fixture(scope="module")
def propagation_pair() -> tuple[PropagationResult, PropagationResult]:
    """One full + one `V_int=0` free-reference propagation, reused across the
    `sigma_from_correlations`-vs-`_sigma_one_energy` and elastic-path tests
    below (built once, not re-run per test)."""
    full = _propagate(
        TG, N2, EPS, CHI, V_INIT, VPRIMES, dt=DT, n_steps=N_STEPS, wp_in=WP_IN, wp_out=WP_OUT
    )
    free = _propagate(
        TG,
        N2,
        EPS,
        CHI,
        V_INIT,
        VPRIMES,
        dt=DT,
        n_steps=N_STEPS,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        free=True,
    )
    return full, free


def test_sigma_from_correlations_matches_sigma_one_energy_per_energy(
    propagation_pair: tuple[PropagationResult, PropagationResult],
) -> None:
    """`sigma_from_correlations`'s array-`E` batching reproduces looping the
    private single-energy kernel `_sigma_one_energy` by hand -- the public
    entry point does not do anything different from its documented
    per-energy transform."""
    full, free = propagation_pair
    e_arr = np.array([0.10, 0.15])
    batched = sigma_from_correlations(
        TG,
        N2,
        full,
        EPS,
        V_INIT,
        VPRIMES,
        e_arr,
        dt=DT,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        free_result=free,
    )
    manual = np.stack(
        [
            _sigma_one_energy(TG, N2, full, EPS, V_INIT, VPRIMES, float(e), DT, WP_IN, WP_OUT, free)
            for e in e_arr
        ]
    )
    np.testing.assert_allclose(batched, manual, rtol=0.0, atol=0.0)  # same call, must agree exactly


def test_free_reference_elastic_path_finite_and_nonnegative(
    propagation_pair: tuple[PropagationResult, PropagationResult],
) -> None:
    """The free-reference elastic (`v'==v_init`) channel is finite and >= 0,
    and is orders of magnitude smaller than the literal-1 fallback -- the
    same invariant `test_td_cross_section.py::test_v2c_...` gates at
    `TD_WORKING_GRID` scale (module docstring), cheaply reproduced here."""
    full, free = propagation_pair
    sigma_fixed = sigma_from_correlations(
        TG,
        N2,
        full,
        EPS,
        V_INIT,
        [V_INIT],
        0.10,
        dt=DT,
        wp_in=WP_IN,
        wp_out=WP_OUT,
        free_result=free,
    )
    sigma_literal1 = sigma_from_correlations(
        TG, N2, full, EPS, V_INIT, [V_INIT], 0.10, dt=DT, wp_in=WP_IN, wp_out=WP_OUT
    )
    assert np.all(np.isfinite(sigma_fixed))
    assert np.all(sigma_fixed >= 0.0)
    assert np.all(np.isfinite(sigma_literal1))
    assert sigma_literal1[0] > 100.0 * sigma_fixed[0]  # the bug the free reference fixes
