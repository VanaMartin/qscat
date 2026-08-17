"""Tests for the NRM per-R ingredients (PRA 77 Eq. 56-59)."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.core.grids import electronic_grid
from qscat.core.nrm import ingredients as ingredients_module
from qscat.core.nrm.discrete_state import (
    AsymptoticDiscreteState,
    PhysicalDiscreteState,
    electronic_hamiltonian,
)
from qscat.core.nrm.ingredients import (
    _c_normalize_columns,
    _drop_null_mode,
    nrm_ingredients,
)
from qscat.dvr import eigen
from qscat.exceptions import ConvergenceError
from qscat.linalg import c_product
from qscat.model import F2


@pytest.fixture(scope="module")
def setup():
    g = electronic_grid(r_max=16.0, order=8, n_complex=6)
    ds = AsymptoticDiscreteState(g, F2, R_inf=g.R0)
    R = np.linspace(6.0, 1.8, 25)
    return g, ds, R


def test_shapes_and_null_mode_dropped(setup):
    g, ds, R = setup
    ing = nrm_ingredients(g, F2, ds, R)
    assert ing.E_n.shape == (R.size, g.n - 1)
    assert ing.V_dn.shape == (R.size, g.n - 1)
    assert ing.v_d_discrete.shape == (R.size,)
    # The dropped mode was the null one: nothing left is near zero AND
    # strongly overlapping phi_d.
    assert np.min(np.abs(ing.E_n)) > 1e-8


def test_v_d_discrete_is_eq_20(setup):
    """V_d(R) = V_0(R) + <phi_d|H_el|phi_d> (Eq. 20), against an INDEPENDENT
    oracle rather than the implementation's own expression.

    Recomputing `expected = v0(R) + d @ (h @ d)` with the exact same code
    path `nrm_ingredients` itself uses is a tautology -- it would still pass
    if that shared expression were simply wrong. Instead, at `R = R_inf`
    (choice B's bound-state R), `<phi_d|H_el|phi_d>` must equal the bound
    eigenvalue `eps_b` that `anion_electronic_states` finds completely
    independently (its own `qscat.dvr.eigen` call on the FULL, unprojected
    `H_el`, not `nrm_ingredients`'s `P H_el P` machinery at all).
    """
    g, ds, _R = setup
    R_inf = g.R0
    ing = nrm_ingredients(g, F2, ds, np.array([R_inf]))
    eps_b, _phi = anion_electronic_states(g, F2, R_inf, n_states=1)
    expected = complex(F2.v0(R_inf)) + complex(eps_b[0])
    # Measured agreement is ~2e-10 (independent eigensolves on different
    # Hamiltonians -- the full H_el here vs. P H_el P inside nrm_ingredients
    # -- so this is not a round-off-identical comparison); 1e-9 leaves margin
    # without loosening past what the two routes actually agree to.
    assert abs(ing.v_d_discrete[0] - expected) < 1e-9

    # R_inf alone does not exercise the `v0(R)` ADDITION: F2's v0 is defined
    # to vanish at large R (measured v0(R_inf) ~ -2e-10, already below the
    # tolerance above), so dropping the `+ v0(R)` term from Eq. (20) entirely
    # would NOT have failed the check above -- verified directly. Cover that
    # separately at an R where v0 is not negligible (measured v0(2.5) ~
    # -0.053 Ha): the `<phi_d|H_el|phi_d>` part is recomputed the same way
    # nrm_ingredients does internally (a partial tautology for that term,
    # unavoidable without a second independent H_el solve at this R), but the
    # v0(R) contribution itself is isolated by subtraction and checked
    # against the model's OWN v0 function directly.
    j = 5
    d = ds.phi_d(float(_R[j]))
    h = electronic_hamiltonian(g, F2, float(_R[j]))
    ing_mid = nrm_ingredients(g, F2, ds, _R)
    electronic_part = d @ (h @ d)
    v0_contribution = ing_mid.v_d_discrete[j] - electronic_part
    assert abs(v0_contribution - complex(F2.v0(float(_R[j])))) < 1e-10


def test_p_space_eigenvectors_are_c_normalized(setup):
    """The paper's explicit normalization convention (p. 012710-6, "we have
    to use for the wave functions the scalar product defined without complex
    conjugation") is checked directly, not inferred from finiteness.

    Swapping `_c_normalize_columns` for numpy's Hermitian `v^dagger v = 1`
    normalization (what `eigen` returns raw) still gives finite, plausible
    E_n/V_dn -- every OTHER test in this file still passes against that
    defect (verified: see the task report's fail/restore demonstration).
    Only checking `c_product(v, v) == 1` directly catches it.
    """
    g, ds, R = setup
    d = ds.phi_d(float(R[0]))
    h_el = electronic_hamiltonian(g, F2, float(R[0]))
    p = np.eye(g.n, dtype=np.complex128) - np.outer(d, d)
    php = p @ h_el @ p
    _evals, vecs = eigen(php)
    vecs = _c_normalize_columns(vecs)
    norms2 = np.array([c_product(vecs[:, i], vecs[:, i]) for i in range(vecs.shape[1])])
    assert np.allclose(norms2, 1.0, atol=1e-8)


def test_v_dn_reflects_c_normalized_eigenvectors(setup):
    """`nrm_ingredients`'s OUTPUT must actually reflect c-normalization, not
    just the standalone helper working when called in isolation.

    `test_p_space_eigenvectors_are_c_normalized` calls `_c_normalize_columns`
    directly, so it would still pass even if `nrm_ingredients`'s own loop
    silently stopped calling it (verified: see the task report's fail/
    restore demonstration). This test instead recomputes `V_dn` at one R
    node from the SAME raw `eigen()` output, explicitly c-normalizing, and
    requires the MAGNITUDES to match `nrm_ingredients`'s own value --
    comparing magnitudes (not the full complex value) so this is decoupled
    from adiabatic tracking/sign alignment (covered by other tests) and
    isolates the normalization SCALE specifically. A Hermitian-normalized
    eigenvector generally has a different c-product scale, so this fails if
    normalization is dropped.
    """
    g, ds, R = setup
    j = 5
    ing = nrm_ingredients(g, F2, ds, R)
    d = ds.phi_d(float(R[j]))
    h_el = electronic_hamiltonian(g, F2, float(R[j]))
    p = np.eye(g.n, dtype=np.complex128) - np.outer(d, d)
    php = p @ h_el @ p
    evals, vecs = eigen(php)
    vecs = _c_normalize_columns(vecs)
    evals, vecs = _drop_null_mode(evals, vecs, d)
    expected_mag = np.sort(np.abs(d @ (h_el @ vecs)))
    actual_mag = np.sort(np.abs(ing.V_dn[j]))
    assert np.allclose(actual_mag, expected_mag, rtol=1e-6, atol=1e-9)


def test_eigenvalues_are_continuous_in_R(setup):
    """Adiabatic tracking: E_n(R) must not jump between adjacent nodes.

    Without tracking, the eigensolver's ordering permutes states across R and
    E_n acquires discontinuities that corrupt F(E,R,R').
    """
    g, ds, R = setup
    ing = nrm_ingredients(g, F2, ds, R)
    low = ing.E_n[:, :20]
    jumps = np.abs(np.diff(low, axis=0))
    spread = np.abs(low).max()
    assert np.max(jumps) < 0.2 * spread


def test_v_dn_continuous_in_R(setup):
    """Sign-tracking: V_dn must not flip sign between adjacent R nodes.

    c-normalization fixes an eigenvector only up to c^2=1, i.e. +-1 -- the
    bilinear c-product admits no other unit-modulus solution. `_track`
    matches state IDENTITY by nearest eigenvalue but says nothing about that
    residual sign, so an untracked sign lets np.linalg.eig flip a state's
    orientation from one R to the next. Measured on this exact fixture
    pre-fix: c_product(phi_n(R_j), phi_n(R_{j+1})) = -0.995 to -1.000 at
    several (n, R) with |V_dn| ~ 0.1 (not tail noise) -- a full flip makes
    the jump ratio below ~2 (|v-(-v)|/|v| = 2). Post-fix the measured worst
    ratio on this fixture is ~0.51 (a genuine, smooth physical variation,
    not a flip), so 1.0 discriminates cleanly between the two.
    """
    g, ds, R = setup
    ing = nrm_ingredients(g, F2, ds, R)
    v = ing.V_dn[:, :100]
    jumps = np.abs(np.diff(v, axis=0))
    local_scale = np.maximum(np.abs(v[:-1]), np.abs(v[1:]))
    # Only compare where the coupling is not already near the numerical
    # floor (a near-zero state's "jump" is dominated by noise, not sign).
    mask = local_scale > 1e-3
    assert mask.sum() > 100, "fixture no longer exercises enough states"
    ratio = jumps[mask] / local_scale[mask]
    assert np.max(ratio) < 1.0, f"V_dn jump ratio {np.max(ratio):.3g} looks like a sign flip"


def test_v_dn_is_finite_and_decays_for_high_states(setup):
    """V_dn stays finite and its magnitude falls off for the highest,
    most-oscillatory states.

    (c-normalization itself is checked directly, against an independent
    Hermitian-normalized recomputation, by
    `test_p_space_eigenvectors_are_c_normalized` and
    `test_v_dn_reflects_c_normalized_eigenvectors` above -- this test does
    not exercise that property.)
    """
    g, ds, R = setup
    ing = nrm_ingredients(g, F2, ds, R)
    assert np.all(np.isfinite(ing.V_dn))
    # The coupling must decay for the highest (most oscillatory) states.
    assert np.abs(ing.V_dn[:, -1]).max() < np.abs(ing.V_dn[:, :10]).max()


def test_coupling_decays_at_large_R(setup):
    """Eq. (67) again, now for the discretized couplings."""
    g, ds, R = setup
    ing = nrm_ingredients(g, F2, ds, R)
    inner = np.abs(ing.V_dn[-1, :20]).max()  # R = 1.8
    outer = np.abs(ing.V_dn[0, :20]).max()  # R = 6.0
    assert outer < 0.5 * inner


def test_min_overlap_warns_on_a_synthetic_tracking_failure(setup, monkeypatch):
    """`nrm_ingredients` must warn when `_sign_align`'s overlap looks like a
    tracking failure (Task 2's fix: NO's physical discrete state hits
    `min|c_product| = 3.3e-15` at the crossing -- see
    docs/physics/nonlocal-resonance-model.md Sec. 10 -- and that failure
    previously had no signal at all, only a docstring caveat).

    Monkeypatching `_sign_align` to report a fabricated near-zero overlap,
    while leaving its actual (correct) vectors untouched, isolates the
    warning's WIRING in `nrm_ingredients` from needing a genuinely
    pathological electronic structure to reproduce the failure here.
    """
    g, ds, R = setup
    real_sign_align = ingredients_module._sign_align
    calls = [0]

    def fabricated_sign_align(prev_vecs, vecs):
        out, _overlap = real_sign_align(prev_vecs, vecs)
        calls[0] += 1
        overlap = 1e-10 if calls[0] == 1 else _overlap
        return out, overlap

    monkeypatch.setattr(ingredients_module, "_sign_align", fabricated_sign_align)
    with pytest.warns(UserWarning, match="tracking failure"):
        ing = nrm_ingredients(g, F2, ds, R)
    assert ing.min_overlap == pytest.approx(1e-10)


def test_min_overlap_does_not_warn_on_a_clean_run(setup):
    """The unmodified fixture must not warn -- `warnings.simplefilter("error")`
    turns any warning into a failure, so this fails loudly if the threshold in
    `_MIN_OVERLAP_WARN` is ever miscalibrated against a genuine, clean run."""
    g, ds, R = setup
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ing = nrm_ingredients(g, F2, ds, R)
    assert ing.min_overlap > 0.5


def test_rejects_ascending_R(setup):
    g, ds, _R = setup
    with pytest.raises(ValueError, match="descending"):
        nrm_ingredients(g, F2, ds, np.array([1.8, 3.0, 6.0]))


def test_p_projector_is_idempotent_and_annihilates_d(setup):
    """Verify the Q = outer(d, d) claim (Task 3's "cannot verify" item).

    A DVR coefficient is `c_j = phi(r_j) sqrt(w_j)`, so Eq. (58)'s projector
    is claimed to need no explicit weights: `Q = outer(d, d)` under the
    c-product alone. If that were wrong, `Q` would not be idempotent (a
    projector squared must equal itself) and `P = I - Q` would not annihilate
    `d`. This is the first eigenproblem in the package to actually rely on
    that claim, so it is checked directly rather than inherited.
    """
    g, ds, R = setup
    d = ds.phi_d(float(R[0]))
    q = np.outer(d, d)
    p = np.eye(g.n, dtype=np.complex128) - q

    # Q is idempotent under the c-product: Q @ Q == Q. This uses c_product's
    # bilinear pairing implicitly (outer(d, d) @ v == d * c_product(d, v)),
    # not np.vdot / conjugation.
    assert np.max(np.abs(q @ q - q)) < 1e-10 * max(1.0, np.max(np.abs(q)))

    # P annihilates d exactly (up to the c-normalization already enforced by
    # phi_d): (I - outer(d,d)) @ d == d - d * c_product(d, d) == 0 since
    # c_product(d, d) == 1.
    assert np.max(np.abs(p @ d)) < 1e-10
    assert abs(c_product(d, d) - 1.0) < 1e-10


def test_null_mode_ambiguity_raises_with_a_corrupted_projector(setup):
    """A degraded null-mode identification must raise, not silently guess.

    Feeding `_drop_null_mode` a `d` that is NOT the actual null eigenvector of
    `P H_el P` (here, a fabricated, unrelated unit vector) must not find a
    clean match: neither the smallest-|E| nor the largest-overlap criterion
    should agree with a real null mode, so the ambiguity/threshold guards
    must fire. This exercises the branch that a "does it run" test would
    never reach.
    """
    g, ds, R = setup
    d = ds.phi_d(float(R[0]))
    h_el = electronic_hamiltonian(g, F2, float(R[0]))
    p = np.eye(g.n, dtype=np.complex128) - np.outer(d, d)
    php = p @ h_el @ p
    evals, vecs = eigen(php)
    norms2 = np.einsum("ij,ij->j", vecs, vecs)
    vecs = vecs / np.sqrt(norms2)

    # A vector orthogonal (in the c-product sense, approximately) to d and
    # unrelated to any eigenvector -- e.g. a different DVR basis vector deep
    # in the ECS tail -- should not overlap ANY eigenvector strongly enough
    # to pass the _NULL_OVERLAP_MIN gate, let alone agree with the
    # smallest-|E| criterion.
    fake_d = np.zeros(g.n, dtype=np.complex128)
    fake_d[-1] = 1.0
    with pytest.raises(ConvergenceError):
        _drop_null_mode(evals, vecs, fake_d)


def test_tracking_and_sign_alignment_are_invariant_to_raw_eigenpair_order(setup, monkeypatch):
    """Feeding a randomly PERMUTED raw eigendecomposition at every R node
    AFTER the first must not change the final E_n/V_dn at all -- tracking
    plus sign alignment must fully undo it.

    The FIRST node's raw order is left untouched deliberately: with no
    predecessor to track against, `nrm_ingredients` has no way to anchor
    column `n` to a physical state at that node, so which physical state
    lands in which output column there is genuinely arbitrary (confirmed
    directly: permuting node 0 too changes the whole trajectory, because
    every later node is then correctly tracked onto that DIFFERENT, but
    self-consistent, labeling -- not a tracking bug). From node 1 onward the
    label IS pinned by the previous node, so a correct tracker must recover
    the exact same per-n assignment and sign as the unpermuted run.

    This is the "teeth" this suite was missing: on this fixture the raw
    ordering `eigen()` returns is already near-adiabatic (only 1 of 24 steps
    is a non-identity permutation), so silently no-op'ing `_track` still
    passed every other test here. Injecting an explicit random permutation
    of `eigen`'s output at every node but the first is a much stronger
    probe -- any tracking or sign bug now shows up as a real numerical
    difference against the unpermuted baseline, not a coincidence of this
    particular fixture.
    """
    g, ds, R = setup
    baseline = nrm_ingredients(g, F2, ds, R)

    rng = np.random.default_rng(0)
    real_eigen = ingredients_module.eigen
    call_count = [0]

    def shuffled_eigen(H):
        evals, vecs = real_eigen(H)
        call_count[0] += 1
        if call_count[0] == 1:
            return evals, vecs  # anchor node: leave the labeling untouched
        perm = rng.permutation(evals.size)
        return evals[perm], vecs[:, perm]

    monkeypatch.setattr(ingredients_module, "eigen", shuffled_eigen)
    shuffled = nrm_ingredients(g, F2, ds, R)

    assert np.allclose(shuffled.E_n, baseline.E_n, rtol=1e-8, atol=1e-10)
    assert np.allclose(shuffled.V_dn, baseline.V_dn, rtol=1e-8, atol=1e-10)
    assert np.allclose(shuffled.v_d_discrete, baseline.v_d_discrete, rtol=1e-8, atol=1e-10)


@pytest.fixture(scope="module")
def crossing_setup():
    """`PhysicalDiscreteState` across F2's bound/resonance crossing (R~2.6),
    where Task 4 found 1.08% of `phi_d`'s c-norm sitting on the ECS contour --
    the actual hazard scenario, not `AsymptoticDiscreteState`'s R-independent
    state (which never gets close to that leakage)."""
    g = electronic_grid(r_max=16.0, order=8, n_complex=6)
    gb = electronic_grid(r_max=16.0, order=8, n_complex=6, angle_deg=40.0)
    R = np.array([3.0, 2.9, 2.8, 2.75, 2.7, 2.65, 2.6, 2.58, 2.56])
    ds = PhysicalDiscreteState(g, F2, R, elec_grid_b=gb, r_d=10.0)
    return g, ds, R


def test_v_dn_benign_across_the_ecs_leakage_crossing(crossing_setup):
    """V_dn stays finite, and physically small, right through the point
    where P's short-rangedness premise is violated (R=2.6, 1.08% ECS
    leakage; see `discrete_state.PhysicalDiscreteState`'s docstring note).

    A bound `phi_d` is (to numerical precision) an eigenvector of the full
    `H_el`, and `d @ P == 0` exactly (see
    `test_p_projector_is_idempotent_and_annihilates_d`), so Eq. (59)'s
    `V_dn = <phi_d|H_el|phi_n> = E_d <phi_d|phi_n>` vanishes identically for
    every bound-branch R -- including R=2.6, the worst-leakage point. The
    leakage shows up only as a smooth few-times growth in that near-zero
    residual (not a blowup) and the coupling only becomes physically
    meaningful one step inward, where the state crosses onto the
    scattering branch. That is the finding: the projector's stated
    short-rangedness premise is technically violated at R=2.6, but V_dn is
    unaffected in any way that matters here.
    """
    g, ds, R = crossing_setup
    ing = nrm_ingredients(g, F2, ds, R)
    assert np.all(np.isfinite(ing.V_dn))
    assert np.all(np.isfinite(ing.E_n))
    mags = np.abs(ing.V_dn[:, :20]).max(axis=1)
    bound = ~ds.used_scattering
    assert bound.sum() >= 3, "fixture no longer straddles the crossing"
    assert (~bound).sum() >= 2, "fixture no longer straddles the crossing"
    assert np.all(mags[bound] < 1e-8), (
        f"bound-branch V_dn should stay near machine-zero even at the "
        f"leakiest point; got max {mags[bound].max():.3g}"
    )
    assert np.all(mags[~bound] > 1e-3), (
        "scattering-branch V_dn should be physically finite, not a residual"
    )
