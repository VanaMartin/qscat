"""Dependency-free oracle checks; external-driver controls run in the pinned image."""

import numpy as np
import pytest

from projects.ukrmol_co.singlet_a2_reference import (
    spin_blocks,
    synthetic_integrals,
    two_electron_matrix,
)


def test_second_quantized_oracle_matches_both_analytic_spin_blocks():
    h1, eri, labels = synthetic_integrals()
    matrix = two_electron_matrix(h1, eri)
    singlet, triplet = spin_blocks(matrix, labels)
    np.testing.assert_allclose(singlet, [[0.42, 0.13], [0.13, 0.23]], atol=1e-13, rtol=0)
    np.testing.assert_allclose(triplet, [[0.18, 0.05], [0.05, 0.07]], atol=1e-13, rtol=0)


def test_hamiltonian_commutes_with_spin_and_respects_c2v():
    h1, eri, labels = synthetic_integrals()
    matrix = two_electron_matrix(h1, eri)
    swap = np.eye(16).reshape(16, 4, 4).transpose(0, 2, 1).reshape(16, 16)
    spin = np.eye(16) - swap
    np.testing.assert_allclose(matrix, matrix.T, atol=1e-14, rtol=0)
    np.testing.assert_allclose(matrix @ spin, spin @ matrix, atol=1e-14, rtol=0)
    sectors = (labels[:, None] ^ labels).ravel()
    np.testing.assert_allclose(matrix[sectors[:, None] != sectors], 0, atol=0, rtol=0)
    np.testing.assert_allclose(eri, eri.transpose(1, 0, 2, 3), atol=0, rtol=0)
    np.testing.assert_allclose(eri, eri.transpose(2, 3, 0, 1), atol=0, rtol=0)
    assert np.linalg.eigvalsh(eri.reshape(16, 16)).min() > -1e-12


def test_exact_triplet_is_not_a_singlet_and_perturbation_has_residual():
    h1, eri, labels = synthetic_integrals()
    matrix = two_electron_matrix(h1, eri)
    _, triplet = spin_blocks(matrix, labels)
    energy, u = np.linalg.eigh(triplet)
    ci = np.zeros((4, 4))
    for amplitude, (i, j) in zip(u[:, 0], [(0, 3), (1, 2)], strict=True):
        ci[i, j] = amplitude / np.sqrt(2)
        ci[j, i] = -amplitude / np.sqrt(2)
    assert np.linalg.norm(matrix @ ci.ravel() - energy[0] * ci.ravel()) < 1e-12
    np.testing.assert_allclose(np.vdot(ci, ci - ci.T), 2, atol=1e-12, rtol=0)
    wrong = ci.copy()
    wrong[0, 3] += 0.01
    assert np.linalg.norm(matrix @ wrong.ravel() - energy[0] * wrong.ravel()) > 1e-3


def test_bad_integral_shapes_reject():
    with pytest.raises(ValueError):
        two_electron_matrix(np.eye(2), np.zeros((3, 3, 3, 3)))
