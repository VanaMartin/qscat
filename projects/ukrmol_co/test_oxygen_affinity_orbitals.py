"""Analytic tensor/rotation invariants for the saved-state semicanonical control."""

import unittest
from pathlib import Path

import numpy as np

from projects.ukrmol_co.oxygen_affinity_orbitals import (
    authenticate,
    block_rotation,
    rotate_amplitudes,
)


def test_spin_fock_blocks_diagonalize_without_occupied_virtual_or_irrep_mixing():
    labels = np.array([0, 0, 1, 0, 0, 1])
    fock = np.diag([-2.0, -1.0, -0.5, 0.2, 0.6, 0.9])
    fock[0, 1] = fock[1, 0] = 0.12
    fock[3, 4] = fock[4, 3] = 0.21
    fock[1, 3] = fock[3, 1] = 0.06
    u = block_rotation(fock, 3, labels)
    np.testing.assert_allclose(u.T @ u, np.eye(6), atol=1e-14, rtol=0)
    assert np.max(np.abs(u[:3, 3:])) == 0
    assert np.max(np.abs(u[labels[:, None] != labels[None, :]])) == 0
    transformed = u.T @ fock @ u
    for block in [transformed[:3, :3], transformed[3:, 3:]]:
        np.testing.assert_allclose(block, np.diag(np.diag(block)), atol=1e-14, rtol=0)
    density = np.diag([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    np.testing.assert_allclose(u @ density @ u.T, density, atol=1e-14, rtol=0)


def test_invalid_or_nonhermitian_fock_is_rejected():
    with unittest.TestCase().assertRaises(ValueError):
        block_rotation(np.eye(3), 0, [0, 0, 0])
    with unittest.TestCase().assertRaises(AssertionError):
        block_rotation(np.array([[1.0, 2.0], [0.0, 3.0]]), 1, [0, 0])


def test_amplitude_rotation_roundtrip_and_tensor_scalar_invariance():
    rng = np.random.default_rng(20261103)
    no = (2, 1)
    nv = (3, 4)
    rotations = []
    for occupied, virtual in zip(no, nv, strict=True):
        u = np.zeros((occupied + virtual, occupied + virtual))
        for lo, hi in [(0, occupied), (occupied, occupied + virtual)]:
            q, _ = np.linalg.qr(rng.normal(size=(hi - lo, hi - lo)))
            u[lo:hi, lo:hi] = q
        rotations.append(u)
    t1 = tuple(rng.normal(size=(o, v)) for o, v in zip(no, nv, strict=True))
    t2 = tuple(rng.normal(size=shape) for shape in [(2, 2, 3, 3), (2, 1, 3, 4), (1, 1, 4, 4)])
    new1, new2 = rotate_amplitudes(t1, t2, rotations, no)
    assert all(t.flags.c_contiguous for t in new2)
    back1, back2 = rotate_amplitudes(new1, new2, [u.T for u in rotations], no)
    for actual, expected in zip((*back1, *back2), (*t1, *t2), strict=True):
        np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=0)
    g1 = tuple(rng.normal(size=t.shape) for t in t1)
    g2 = tuple(rng.normal(size=t.shape) for t in t2)
    newg1, newg2 = rotate_amplitudes(g1, g2, rotations, no)
    for a, b, c, d in zip((*t1, *t2), (*g1, *g2), (*new1, *new2), (*newg1, *newg2), strict=True):
        np.testing.assert_allclose(np.sum(a * b), np.sum(c * d), atol=1e-12, rtol=0)


def test_release_rejects_expanded_or_iterative_scope_before_quantum_import():
    base = dict(
        source_sha256={},
        input_sha256={},
        execution_authorized=True,
        states=["O", "O-"],
        new_scf_calls=0,
        new_ccsd_iterations=0,
    )
    for change in [
        dict(states=["O", "O-", "C"]),
        dict(new_scf_calls=1),
        dict(new_ccsd_iterations=1),
        dict(execution_authorized=False),
    ]:
        with unittest.TestCase().assertRaises(ValueError):
            authenticate(Path("."), Path("."), base | change)


if __name__ == "__main__":
    suite = unittest.TestSuite(
        unittest.FunctionTestCase(value)
        for name, value in list(globals().items())
        if name.startswith("test_")
    )
    raise SystemExit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
