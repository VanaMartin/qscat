"""Propagation, the half-Fourier transform, and the packet diagnostics."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from qscat.core.grids import nuclear_grid
from qscat.core.nrm.extended import LaunchBasis
from qscat.core.nrm.propagation import propagate_nrm


def _launch(psi0: np.ndarray, e_total: np.ndarray) -> LaunchBasis:
    """A raw-column `LaunchBasis` with `coeffs = I` -- exercises propagation
    and the transform, not the SVD factorization `test_nrm_extended.py`
    already gates."""
    r = psi0.shape[1]
    return LaunchBasis(
        vectors=psi0.astype(np.complex128),
        coeffs=np.eye(r, dtype=np.complex128),
        energies=e_total,
        e_total=e_total,
        truncation_error=0.0,
    )


def test_transform_of_a_single_decaying_mode_is_the_resolvent():
    """A 1x1 'Hamiltonian' h has Psi(t) = e^{-iht}, whose transform is
    -i * i/(E-h) ... i.e. exactly (E-h)^-1. The smallest possible check that
    the -i prefactor and the quadrature weight are both right."""
    h = sp.csr_matrix(np.array([[0.20 - 0.05j]], dtype=np.complex128))
    psi0 = np.array([[1.0 + 0.0j]])
    e_total = np.array([0.10])
    res = propagate_nrm(h, _launch(psi0, e_total), nuclear_grid=None, dt=0.02, n_steps=40000)
    assert res.psi_d[0, 0] == pytest.approx(1.0 / (0.10 - (0.20 - 0.05j)), rel=1e-4)


def test_columns_do_not_talk_to_each_other():
    """Two energies propagated together == each propagated alone."""
    rng = np.random.default_rng(11)
    n = 12
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    # symmetric, absorbing: -6j margin (checked against this seed's spectrum)
    # comfortably clears a+a.T's largest imaginary eigenvalue (~4.96), so
    # every mode of h actually decays under d/dt psi = -i h psi -- a smaller
    # shift (e.g. -3j) is NOT guaranteed to for an arbitrary complex-normal
    # a+a.T and, for this rng, does not (see test_survival_... below).
    h = sp.csr_matrix(a + a.T - 6j * np.eye(n))
    psi0 = (rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))).astype(np.complex128)
    e = np.array([0.05, 0.09])
    both = propagate_nrm(h, _launch(psi0, e), nuclear_grid=None, dt=0.05, n_steps=400)
    for j in range(2):
        one = propagate_nrm(
            h, _launch(psi0[:, [j]], e[[j]]), nuclear_grid=None, dt=0.05, n_steps=400
        )
        assert np.allclose(both.psi_d[:, j], one.psi_d[:, 0], rtol=1e-10)


def test_survival_decays_and_unabsorbed_reports_it():
    rng = np.random.default_rng(3)
    n = 12
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    # -3j is NOT enough margin for this seed: a+a.T's largest imaginary
    # eigenvalue is ~4.58, so a -3j shift leaves a genuinely GROWING mode
    # (verified against scipy.linalg.expm) and survival diverges rather than
    # decays. -6j clears it (largest shifted Im eigenvalue ~-1.42).
    h = sp.csr_matrix(a + a.T - 6j * np.eye(n))
    psi0 = np.ones((n, 1), dtype=np.complex128)
    res = propagate_nrm(h, _launch(psi0, np.array([0.05])), nuclear_grid=None, dt=0.05, n_steps=400)
    assert res.survival[0, 0] > res.survival[-1, 0]
    assert res.unabsorbed[0] == pytest.approx(res.survival[-1, 0])


def test_diagnostics_match_analytic_gaussian_at_t0():
    """A real-envelope Gaussian times a plane wave `exp(i*p0*(R-r0))` has
    EXACT `<R>_0 = r0` and `<P>_0 = p0`: the envelope's own contribution to
    `<-i d/dR>` is purely imaginary-suppressed by the c-product's Hermitian
    structure on the real region (it integrates to zero by symmetry), so
    only the phase gradient `p0` survives. This is the diagnostics-ON path
    (`nuclear_grid` given) that the other three tests never exercise --
    before the `sqrt(w)` fix in `_record`'s Eq. (4.6) term, this measured
    5.93 instead of 1.70 (a ~3.5x, grid-dependent error), not a rounding
    discrepancy.
    """
    grid = nuclear_grid()
    mask = grid.real_points <= grid.R0
    r = grid.real_points[mask]
    r0, sigma, p0 = 6.0, 0.7, 1.7
    psi_val = np.exp(-((r - r0) ** 2) / (2 * sigma**2)) * np.exp(1j * p0 * (r - r0))
    coeff = psi_val * np.sqrt(grid.weights[mask].real)

    psi0 = np.zeros((grid.n, 1), dtype=np.complex128)
    psi0[mask, 0] = coeff
    h_ext = sp.identity(grid.n, format="csr", dtype=np.complex128) * (-3j)

    res = propagate_nrm(
        h_ext, _launch(psi0, np.array([0.05])), nuclear_grid=grid, dt=0.01, n_steps=2
    )

    assert res.survival[0, 0] > 0.0
    assert res.centroid[0, 0] == pytest.approx(r0, abs=1e-9)
    assert res.momentum[0, 0] == pytest.approx(p0, abs=1e-6)


def test_low_rank_reconstruction_matches_dense_per_energy_propagation():
    """A rank-2 launch with random COMPLEX `coeffs`, gated against a dense
    per-energy reconstruction propagated on its own.

    `test_columns_do_not_talk_to_each_other` uses `coeffs = I`, which cannot
    distinguish `coeffs` from `coeffs.conj()` in the reconstruction
    `d = psi[:n_r, :] @ coeffs` (a c-product violation that would silently
    conjugate the SVD coefficients) -- with `coeffs = I` both give the same
    answer. Random complex, non-identity `coeffs` breaks that degeneracy.
    """
    rng = np.random.default_rng(7)
    n = 10
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h = sp.csr_matrix(a + a.T - 6j * np.eye(n))
    vectors = (rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))).astype(np.complex128)
    coeffs = (rng.normal(size=(2, 3)) + 1j * rng.normal(size=(2, 3))).astype(np.complex128)
    e_total = np.array([0.05, 0.07, 0.09])
    launch = LaunchBasis(
        vectors=vectors, coeffs=coeffs, energies=e_total, e_total=e_total, truncation_error=0.0
    )
    batched = propagate_nrm(h, launch, nuclear_grid=None, dt=0.05, n_steps=300)
    for j in range(3):
        psi0_j = (vectors @ coeffs[:, [j]]).astype(np.complex128)
        one = propagate_nrm(
            h, _launch(psi0_j, e_total[[j]]), nuclear_grid=None, dt=0.05, n_steps=300
        )
        assert np.allclose(batched.psi_d[:, j], one.psi_d[:, 0], rtol=1e-8)


# --- the transform identity, on a COMPLETE fixture --------------------------
#
# `Psi_d^TI(R;E) = -i Int_0^inf dt e^{iEt} Psi_d(R,t)` is an identity, not a
# comparison: propagating and transforming must return the vector
# `_psi_d_for_energy` solves for from the SAME ingredients on the SAME grid.
# Two measured facts decide how it can be gated. Every number below was
# measured on 2026-08-19 by the campaign these two tests are the residue of;
# they are recorded HERE because no other shipped file carries them (the NRM
# physics note is still time-independent-only until the TD section lands).
#
# 1. THE PROJECTED-STATE SUM MUST BE COMPLETE (`n_states=None`). Truncating it
#    breaks the premise the transform rests on -- that every eigenmode of
#    `H_ext` decays. `V_dn` and `E_n` are complex (electronic ECS), so the
#    anti-Hermitian part of the arrow matrix is INDEFINITE and a truncated arm
#    set can leave genuinely GROWING eigenmodes: measured `max Im(E)` over the
#    spectrum of `H_ext` on the fixture below is +2.79e-3 at `n_states=3` (56
#    modes with `Im>1e-8`, eigenvalue condition numbers ~1, so they are real
#    eigenvalues and not `eig` noise), against +5.0e-13 for the complete sum
#    -- MARGINALLY stable, i.e. zero to solver tolerance rather than
#    comfortably negative, which is why F2 below can still fail to decay. With
#    a truncated set the propagation diverges instead of converging -- on F2
#    at `n_states=3`, `rel` walks 1.23 (T=2000) -> 1.4e3 (T=64000). The
#    complete sum also restores the autodetachment damping that makes
#    `Psi_d(t)` decay in the first place. The failure is NOT monotone in
#    `n_states` (a small F2 fixture measured here loses its growing modes
#    again by `n_states=12`), so the defensible rule is not "truncation always
#    diverges" but "only the complete sum needs no check": any other value
#    would have to be cleared by a dense `eig`. Dense, specifically --
#    ARPACK `eigs(which="LI")` under-reports `max Im` badly on these strongly
#    non-normal matrices (+4.64e-4 against dense `eig`'s +2.79e-3, 6x low, on
#    the identical 716-dimensional matrix). `propagate_nrm` warns at runtime
#    when the propagated packet ends heavier than it started.
# 2. THE ERROR IS `T`- AND `dt`-SEPARABLE, and the two add IN QUADRATURE:
#    `rel = sqrt(truncation(T)^2 + propagation(dt)^2)` with
#    `truncation(T) = 0.40*sqrt(S(T)/S(0))` and `propagation(dt=1) = 1.43e-4`.
#    Measured at `n_states=None`, `dt=1`: 2.9e-1 (T=500), 4.0e-2 (1e3),
#    4.0e-3 (2e3), 6.5e-4 (3e3), 1.7e-4 (4e3), 1.4e-4 (5e3) -- every row
#    reproduced by that budget to a few percent. Beyond T~5000 the `dt=1`
#    propagation error is all that is left; halving to `dt=0.5` divides it by
#    64 (`dt^6`, the order-3 Pade rate) and the T=5000 point then reads
#    1.65e-5, which is that run's TRUNCATION, no longer hidden. The gate below
#    sits at the T=4000, `dt=1` point.

_N2_REAL_SEGMENTS = ((3, 1.5), (8, 3.0), (2, 4.0), (4, 8.0))
_N2_COMPLEX_SEGMENTS = ((3, 20.0),)


@pytest.fixture(scope="module")
def n2_deck():
    """A small-but-COMPLETE N2 deck: 179 nuclear x 74 electronic points, all
    73 projected states kept.

    Deliberately coarser than the production `N2:emoscat` deck (251 x 107).
    The identity is algebraic -- both routes run on the same ingredients and
    the same grid -- so it does not require a physically converged
    discretisation, only one where the packet actually decays within an
    affordable propagation. The production deck was measured too (H_ext =
    26857): same behaviour, ~8x the cost per unit time, `rel` 7.1e-2 at
    T=1000 and 1.0e-2 at T=4000/`dt=2` (that last is the `dt=2` propagation
    error, not truncation).

    WHY N2 AND NOT F2, even at `n_states=None`. F2's packet does not leave.
    Measured on a coarse F2 nuclear deck with the complete arm set: the
    centroid creeps 2.66 -> 2.77 bohr and then oscillates, `<P>_t` oscillates
    about zero, and survival plateaus at 0.673 of `S(0)` at T=4000 (`rel`
    8.4e-1 / 7.0e-1 / 5.4e-1 at T = 1e3 / 2e3 / 4e3). Two causes, and the
    dominant one is not physics:
    (a) a genuine bound component -- F2's `V_d(R)` has a 0.0223 Ha well
        (minimum -0.149264 Ha at R=3.363 against the F+F- asymptote
        -0.126931), supporting >=24 near-real modes that carry 5.08e-3 of the
        launch norm and contribute ~2.0e-2 relative to the transform. No
        affordable `T` removes those (they would need T >~ 1e7);
    (b) dominant -- GRID-TRAPPED flux. `E_total` sits +0.0991 Ha above the
        anion asymptote, so the dissociating wave has `K_R = 58.5` and a
        wavelength of 0.107 bohr, against ~15 points/bohr on that deck where
        >~70 are needed. The outgoing wave cannot reach the ECS absorber and
        rattles in place. It LOOKS bound and is not.
    So "F2 leaves by dissociation rather than autodetachment" is the wrong
    description of what was measured: on a coarse deck it does not leave at
    all. A TD-DA F2 run needs the FINE production nuclear grid (which resolves
    K~58-78) and a much longer `T` than N2's VE case -- a real cost item, not
    a fixture-selection detail.
    """
    from qscat.core.grids import electronic_grid, segmented_grid
    from qscat.core.nrm.discrete_state import AsymptoticDiscreteState
    from qscat.core.nrm.ingredients import nrm_ingredients
    from qscat.core.vibrational import vibrational_states
    from qscat.model import N2

    nuc = segmented_grid(_N2_REAL_SEGMENTS, _N2_COMPLEX_SEGMENTS, angle_deg=35.0, quadrature=10)
    elec = electronic_grid(r_max=11.0, order=6, n_complex=3)
    phi_d = AsymptoticDiscreteState(elec, N2, R_inf=nuc.R0)
    # nrm_ingredients requires strictly DESCENDING R.
    r_values = nuc.points[nuc.points.imag == 0.0].real[::-1]
    ing = nrm_ingredients(elec, N2, phi_d, r_values)
    eps, chi = vibrational_states(nuc, N2.mu, 4, N2.v0)
    return N2, nuc, elec, phi_d, ing, eps, chi


@pytest.mark.slow
def test_propagated_psi_d_reproduces_the_time_independent_solution(n2_deck):
    """Psi_d^TD(R;E) == Psi_d^TI(R;E), VECTOR to vector.

    Not a cross-section comparison: the two routes must agree on the whole
    nuclear vector, which isolates propagation error from every downstream
    convention. Measured 2026-08-19 at these settings: rel = 1.7264e-04,
    unabsorbed/S(0) = 6.68e-08, 117 s wall on the 12-core dev machine. The
    gate sits at 5e-4 -- 2.9x the achieved error, headroom for BLAS/LU
    variation across platforms and nothing more.
    """
    from qscat.core.nrm.extended import extended_hamiltonian, initial_packet
    from qscat.core.nrm.vibrational_excitation import _psi_d_for_energy

    model, nuc, elec, phi_d, ing, eps, chi = n2_deck
    e_kin, v_init = 0.10, 0

    h = extended_hamiltonian(ing, nuc, model)
    launch = initial_packet(
        nuc, elec, model, phi_d, ing, eps, chi, v_init, np.array([e_kin]), rank_tol=1e-10
    )
    res = propagate_nrm(h, launch, nuc, dt=1.0, n_steps=4000, order=3)

    want = _psi_d_for_energy(nuc, elec, model, phi_d, eps, chi, v_init, e_kin, ing, None)
    got = res.psi_d[:, 0]

    assert res.unabsorbed[0] < 1e-6 * res.survival[0, 0], "packet not yet absorbed"
    rel = np.linalg.norm(got - want) / np.linalg.norm(want)
    assert rel < 5e-4, f"relative vector error {rel:.3e}"


def test_propagation_and_transform_match_the_exact_finite_time_transform(n2_deck):
    """The machinery alone, at a tolerance the physics gate cannot reach.

    `-i Int_0^T dt e^{iEt} Psi_d(t)` has a closed form in the eigenbasis of
    `H_ext`: `sum_j v_j a_j (1 - e^{i(E-E_j)T}) / (E - E_j)`. Comparing the
    propagated-and-transformed vector against THAT (rather than against the
    `T -> infinity` TI solve) removes time truncation from the comparison
    entirely, leaving only the Pade propagation and the Simpson quadrature.
    Any prefactor, weight, or reconstruction error survives here at 1e-8,
    five orders below what the identity gate above can resolve.

    Run at `n_states=3` deliberately: small enough for a dense `eig`, and the
    truncation's growing modes are irrelevant to a finite-T comparison. The
    measured errors (4.10e-5 / 6.48e-7 / 1.03e-8 at dt = 1 / 0.5 / 0.25, T
    fixed) fall as `dt^6` -- the order-3 diagonal Pade rate -- which the
    ratio assertion pins.
    """
    import numpy.linalg as nla
    from qscat.core.nrm.extended import extended_hamiltonian, initial_packet

    model, nuc, elec, phi_d, ing, eps, chi = n2_deck
    n_states = 3
    h = extended_hamiltonian(ing, nuc, model, n_states=n_states)
    launch = initial_packet(
        nuc,
        elec,
        model,
        phi_d,
        ing,
        eps,
        chi,
        0,
        np.array([0.10]),
        n_states=n_states,
        rank_tol=1e-10,
    )
    vals, vecs = nla.eig(np.asarray(h.todense()))
    # Complex-symmetric: expand the launch vector in the (c-orthogonal, but
    # not orthonormal) eigenbasis by solving rather than projecting.
    a = nla.solve(vecs, launch.vectors[:, 0] * launch.coeffs[0, 0])
    de = float(launch.e_total[0]) - vals

    t_max = 200.0
    exact = (vecs * (a * (1.0 - np.exp(1j * de * t_max)) / de)[None, :]).sum(axis=1)[: nuc.n]
    errs = []
    for dt in (0.5, 0.25):
        # `n_states=3` HAS growing modes here (that is the point of running it
        # truncated), so the runtime guard fires -- correctly. Asserting the
        # warning rather than silencing it makes this test also cover the
        # guard on a real truncated fixture: measured final/minimum survival
        # 1.3598 at both dt.
        with pytest.warns(UserWarning, match="growing eigenmode"):
            res = propagate_nrm(h, launch, nuc, dt=dt, n_steps=int(t_max / dt), order=3)
        errs.append(nla.norm(res.psi_d[:, 0] - exact) / nla.norm(exact))

    assert errs[0] < 2e-6, f"dt=0.5 transform error {errs[0]:.3e}"
    assert errs[1] < 5e-8, f"dt=0.25 transform error {errs[1]:.3e}"
    # dt^6 (order-3 diagonal Pade); measured ratio 63.
    assert errs[0] / errs[1] > 30.0, f"convergence ratio {errs[0] / errs[1]:.1f}"


def test_a_growing_mode_is_warned_about_not_returned_silently():
    """The runtime guard `_warn_if_diverging`.

    A truncated arm set can leave `H_ext` with an eigenvalue in the UPPER
    half-plane (see the comment block above for the measured spectra), and the
    resulting `psi_d` is exponentially wrong rather than under-converged --
    with nothing in the returned object saying so. A 1x1 `H` with `Im > 0` is
    the smallest instance: the packet ends heavier than it started, which is
    the condition the guard reports.
    """
    h = sp.csr_matrix(np.array([[0.20 + 0.01j]], dtype=np.complex128))
    launch = _launch(np.array([[1.0 + 0.0j]]), np.array([0.10]))
    with pytest.warns(UserWarning, match="growing eigenmode"):
        res = propagate_nrm(h, launch, nuclear_grid=None, dt=0.5, n_steps=200)
    assert res.unabsorbed[0] > res.survival[0, 0]


def test_growth_after_absorption_is_warned_about_too():
    """The second half of the guard: growth that ENDS BELOW where it started.

    A packet that absorbs first and grows afterwards leaves
    `survival[-1] < survival[0]` -- the F2 `n_states=3` run at T=8000 reads
    0.59 there while being flat wrong (`rel = 2.3`). Only the comparison
    against the run's MINIMUM catches it. Modelled here by two modes: a large
    one that decays fast and a small one that grows slowly.
    """
    h = sp.csr_matrix(np.diag(np.array([0.20 - 0.025j, 0.30 + 0.001j], dtype=np.complex128)))
    launch = _launch(np.array([[1.0], [0.01]], dtype=np.complex128), np.array([0.10]))
    with pytest.warns(UserWarning, match="growing eigenmode"):
        res = propagate_nrm(h, launch, nuclear_grid=None, dt=0.5, n_steps=3000)
    assert res.unabsorbed[0] < res.survival[0, 0], "this case must END below its start"
    assert res.unabsorbed[0] > 1.2 * res.survival[:, 0].min()
