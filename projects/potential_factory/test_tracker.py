# projects/potential_factory/test_tracker.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import N2

from projects.potential_factory.tracker import (
    DEFAULT_RES_WINDOW,
    ElectronicPair,
    Pole,
    WellParams,
    pole_sensitivity,
    solve_pole_params,
    track_curve,
    well_potential,
)

EV = 27.211386


@pytest.fixture(scope="module")
def pair():
    return ElectronicPair()


def _v_n2(R):
    # Bare-well convention: ElectronicPair.pole windows are absolute, so the
    # potential passed in must be on the same scale as the window -- here the
    # bare well (v_int + centrifugal), WITHOUT v0(R).
    def v(r):
        return N2.surface(r, R) - N2.v0(R)

    return v


def test_n2_pole_at_R0_matches_documented_values(pair):
    # docs/physics/n2-resonance.md: E_res(R0) = 2.445 eV, Gamma(R0) = 0.455 eV.
    p = pair.pole(_v_n2(N2.R0), DEFAULT_RES_WINDOW)
    assert isinstance(p, Pole)
    shift_eV = p.shift * EV
    assert abs(shift_eV - 2.445) < 0.02
    assert abs(p.gamma * EV - 0.455) < 0.02
    assert p.residual < 0.05 * p.gamma


def test_gate_rejects_the_spike_fake_pole(pair):
    # Bare l=1 Gaussian well with lam=0.5, alpha=0.3 has NO resonance; the
    # 2026-08-24 spike saw find_resonance_pole return (0.11 eV, 0.26 eV) here.
    def v(r):
        rr = np.asarray(r, dtype=np.complex128)
        return -0.5 * np.exp(-0.3 * rr**2) + 1.0 / rr**2

    assert pair.pole(v, DEFAULT_RES_WINDOW) is None


def test_bound_state_is_accepted_through_bound_window(pair):
    from projects.potential_factory.tracker import DEFAULT_BOUND_WINDOW

    # N2 at R = 3.0 bohr: the anion is bound (Gamma == 0, shift < 0) relative
    # to the neutral curve v0(R) -- the bare-well convention means "bound"
    # shows up directly as a negative shift, not a comparison to v0(R).
    p = pair.pole(_v_n2(3.0), DEFAULT_BOUND_WINDOW)
    # A bound state's Im E is round-off (Linux BLAS gives ~1e-15), never exactly 0.
    assert p is not None and p.gamma < 1e-12 and p.shift < 0.0


def test_c_product_gradient_matches_finite_difference(pair):
    Ha, _ = pair.hamiltonians(_v_n2(N2.R0))
    p = pair.pole(_v_n2(N2.R0), DEFAULT_RES_WINDOW)
    assert p is not None
    dEdV = pole_sensitivity(Ha, p.energy)
    # find_resonance_pole returns the two-angle midpoint (Ea+Eb)/2, not an
    # eigenvalue of Ha itself; the derivative must be taken on ONE matrix's
    # own eigenvalue, so the finite difference is referenced to Ha's own
    # matched eigenvalue E0, not to p.energy.
    E0 = eigen_nearest(Ha, p.energy)
    # perturb one real-region diagonal entry and re-find the nearest eigenvalue
    i = int(np.argmin(np.abs(pair.grid_a.points - 1.5)))
    h = 1e-5
    Hp = Ha.copy()
    Hp[i, i] += h
    Ep = eigen_nearest(Hp, E0)
    fd = (Ep - E0) / h
    assert abs(fd - dEdV[i]) < 1e-3 * max(1.0, abs(dEdV[i]))


def eigen_nearest(H, E):
    from qscat.dvr import eigen

    vals = eigen(H)[0]
    return vals[int(np.argmin(np.abs(vals - E)))]


def test_newton_recovers_n2_well_parameters_from_a_perturbed_seed(pair):
    R = N2.R0
    lam_true = float(N2.lam(R).real)
    v_true = well_potential(N2.ell, lam_true, N2.alpha_c, None)
    target = pair.pole(v_true, DEFAULT_RES_WINDOW)
    assert target is not None
    sol, pole = solve_pole_params(pair, N2.ell, target.energy, WellParams(lam=4.5, alpha=0.5))
    assert abs(sol.lam - lam_true) < 1e-5
    assert abs(sol.alpha - N2.alpha_c) < 1e-5
    assert abs(pole.energy - target.energy) < 1e-8


def test_newton_bound_target_solves_lam_only(pair):
    R = 3.0
    lam_true = float(N2.lam(R).real)
    v_true = well_potential(N2.ell, lam_true, N2.alpha_c, None)
    from projects.potential_factory.tracker import DEFAULT_BOUND_WINDOW

    target = pair.pole(v_true, DEFAULT_BOUND_WINDOW)
    assert target is not None
    # 0.9*lam_true has NO bound state at R=3.0 (the d-wave well
    # unbinds somewhere between 0.93 and 0.94*lam_true, measured) -- 0.95 is
    # the smallest round perturbation that still seeds a real pole.
    seed = WellParams(lam=lam_true * 0.95, alpha=N2.alpha_c)
    sol, _ = solve_pole_params(pair, N2.ell, complex(target.energy.real, 0.0), seed)
    assert abs(sol.lam - lam_true) < 1e-5 and sol.alpha == N2.alpha_c


def test_newton_bound_target_converges_despite_imaginary_noise(pair):
    # A bound target's Im carries only the
    # ECS pole-finder's own irreducible numerical noise (measured ~1e-7 on a
    # genuinely real state, at F2's R=2.617) -- the bound branch's lam-only
    # step cannot and should not
    # try to correct it. Before this fix, `abs(f) < tol` on the FULL complex
    # residual made that noise floor unreachable and Newton stalled even
    # after `f.real` converged to round-off (the exact failure mode that
    # tripped `track_curve` there). Perturb a genuinely bound target by
    # `-1e-7j` (larger than the default `tol=1e-8`, so the OLD check could
    # never be satisfied) and confirm Newton converges instead of raising.
    R = 3.0
    lam_true = float(N2.lam(R).real)
    v_true = well_potential(N2.ell, lam_true, N2.alpha_c, None)
    from projects.potential_factory.tracker import DEFAULT_BOUND_WINDOW

    target = pair.pole(v_true, DEFAULT_BOUND_WINDOW)
    assert target is not None
    noisy_target = complex(target.energy.real, -1e-7)
    seed = WellParams(lam=lam_true * 0.95, alpha=N2.alpha_c)
    sol, pole = solve_pole_params(pair, N2.ell, noisy_target, seed)
    assert abs(sol.lam - lam_true) < 1e-5 and sol.alpha == N2.alpha_c
    assert abs(pole.energy.real - noisy_target.real) < 1e-8


@pytest.mark.slow
def test_track_curve_recovers_n2_lam_and_alpha_over_R(pair):
    from projects.potential_factory.tracker import DEFAULT_BOUND_WINDOW

    R_desc = np.linspace(3.0, 1.6, 15)

    def target(R):
        # N2's anion state is still bound at R=2.5 (it crosses into
        # the continuum near R~2.3-2.4), and right at the crossing the pole
        # sits below the gate's threshold floor (e_floor=0.006 Ha) so it is
        # legitimately gated out -- try bound first, then resonance, and
        # signal "no target" with NaN rather than asserting.
        #
        # Boundedness is decided by WHICH window found the pole, not by
        # `p.gamma == 0.0`: ElectronicPair.pole's bound branch only requires
        # abs(imag) < 1e-6, so a genuinely bound pole can carry floating-point
        # imaginary noise of EITHER sign -- when it's negative, gamma comes
        # out as a tiny positive float, not exactly 0.0, and the equality
        # check misclassifies a bound target as resonant (verified: R=2.9,
        # 2.8, 2.5 all hit this). Zeroing the imaginary part for anything
        # found via the bound window is the robust fix.
        v = well_potential(N2.ell, float(N2.lam(R).real), N2.alpha_c, None)
        p = pair.pole(v, DEFAULT_BOUND_WINDOW)
        if p is not None:
            return complex(p.energy.real, 0.0)
        p = pair.pole(v, DEFAULT_RES_WINDOW)
        if p is None:
            return complex(float("nan"), float("nan"))
        return p.energy

    res = track_curve(
        pair,
        N2.ell,
        R_desc,
        target,
        WellParams(lam=float(N2.lam(3.0).real) * 0.95, alpha=N2.alpha_c),
    )
    assert res.converged.sum() >= 12
    ok = res.converged
    np.testing.assert_allclose(res.lam[ok], N2.lam(R_desc[ok]).real, rtol=1e-4)
    np.testing.assert_allclose(res.alpha[ok], N2.alpha_c, rtol=1e-4)


def test_track_curve_alpha_of_r_overrides_bound_branch_alpha(pair):
    # A BOUND target's Newton step never moves
    # `alpha` (`solve_pole_params` solves `lam` only there), so without
    # `alpha_of_R` the whole bound branch would run at whatever `alpha` the
    # FIRST node's seed carried, not at a per-node `alpha(R)` curve --
    # invisible when that curve happens to be constant, wrong otherwise. Use
    # a non-constant `alpha_of_R` and check the SECOND (bound) node's tracked
    # alpha is exactly that curve's value there, not the seed's.
    from projects.potential_factory.tracker import DEFAULT_BOUND_WINDOW

    R_desc = np.array([3.0, 2.9])

    def alpha_of_R(R: float) -> float:
        return 0.40 + 0.01 * (R - 3.0)  # 0.40 at R=3.0, 0.399 at R=2.9 -- non-constant

    def target_of_R(R: float) -> complex:
        lam_true = float(N2.lam(R).real)
        v = well_potential(N2.ell, lam_true, alpha_of_R(R), None)
        p = pair.pole(v, DEFAULT_BOUND_WINDOW)
        assert p is not None  # both nodes are bound at this (lam, alpha)
        return complex(p.energy.real, 0.0)

    seed = WellParams(lam=float(N2.lam(3.0).real) * 0.95, alpha=alpha_of_R(3.0))
    res = track_curve(pair, N2.ell, R_desc, target_of_R, seed, alpha_of_R=alpha_of_R)
    assert res.converged.all()
    assert res.alpha[1] == alpha_of_R(float(R_desc[1]))


def test_continuity_guard_tracks_upper_bound_state_not_deep_one(pair):
    # A bare l=2 well with lam=12, alpha=0.4 has TWO bound states
    # (verified via eigen below: shallow/"upper" ~-0.0863 Ha, deep/"lower"
    # ~-2.881 Ha). Newton is seeded at the UPPER state's own (lam, alpha)
    # with a resonant target on that SAME state's own trajectory (measured
    # at lam=11.475, past the bound-resonance threshold). Without the
    # continuity guard, `_BRIDGE_WINDOW`'s wide union search has no
    # preference for the pole actually being tracked -- `find_resonance_pole`
    # returns the global residual-argmin over the whole window -- so a far,
    # non-continuous jump onto some other angle-stable state could in
    # principle be silently accepted instead of the upper state Newton is
    # walking. The guard restricts every accepted step to within `_MAX_STEP`
    # of the CURRENT pole, so the solution must stay close to the seed, not
    # land near whatever (lam, alpha) would put the DEEP state at the target.
    from qscat.dvr import eigen

    ell = 2
    lam0, alpha0 = 12.0, 0.4
    v0 = well_potential(ell, lam0, alpha0, None)
    Ha, _ = pair.hamiltonians(v0)
    vals, _ = eigen(Ha)
    bound = vals[(vals.real < 0.0) & (np.abs(vals.imag) < 1e-6)]
    assert bound.size == 2
    reals = sorted(bound.real, reverse=True)
    assert abs(reals[0] - (-0.0863417)) < 1e-3  # shallow/upper state
    assert abs(reals[1] - (-2.8809450)) < 1e-3  # deep/lower state

    target = complex(0.009645461633939076, -7.649505012750466e-05)
    seed = WellParams(lam=lam0, alpha=alpha0)
    sol, pole = solve_pole_params(pair, ell, target, seed)
    assert abs(pole.energy - target) < 1e-6
    # Within 5% of the seed's lam -- nowhere near a value that would put the
    # deep state (lam far from 12.0) at this target instead.
    assert abs(sol.lam - lam0) < 0.05 * lam0
