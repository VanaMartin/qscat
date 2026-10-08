"""Frozen-core Schur reduction against full occupied Slater determinants."""

from itertools import combinations

import numpy as np
import pytest

from projects.ukrmol_co.ci_overlap import frozen_core_overlaps


def occupations(norb, count):
    return sorted(combinations(range(norb), count), key=lambda s: sum(1 << i for i in s))


def full_determinant_oracle(left, right, metric, ncore, electrons):
    ncas = len(metric) - ncore
    alpha, beta = [occupations(ncas, n) for n in electrons]
    result = np.zeros((len(left), len(right)), dtype=complex)
    core = list(range(ncore))
    for a, occ_a in enumerate(alpha):
        for b, occ_b in enumerate(beta):
            for c, occ_c in enumerate(alpha):
                for d, occ_d in enumerate(beta):
                    ia, ib, ic, id_ = [
                        core + [ncore + i for i in s] for s in (occ_a, occ_b, occ_c, occ_d)
                    ]
                    weight = np.linalg.det(metric[np.ix_(ia, ic)]) * np.linalg.det(
                        metric[np.ix_(ib, id_)]
                    )
                    result += np.outer(left[:, a, b].conj(), right[:, c, d]) * weight
    return result


@pytest.mark.parametrize("ncore,electrons", [(0, (2, 1)), (1, (1, 0)), (2, (2, 1))])
def test_complex_overlap_matches_full_slater_determinants(ncore, electrons):
    rng = np.random.default_rng(39)
    ncas = 4
    size = ncore + ncas
    frames = [
        np.linalg.qr(rng.normal(size=(size + 3, size)) + 1j * rng.normal(size=(size + 3, size)))[0]
        for _ in range(2)
    ]
    metric = frames[0].conj().T @ frames[1]
    shape = tuple(len(occupations(ncas, n)) for n in electrons)
    left = rng.normal(size=(2, *shape)) + 1j * rng.normal(size=(2, *shape))
    right = rng.normal(size=(3, *shape)) + 1j * rng.normal(size=(3, *shape))
    actual = frozen_core_overlaps(left, right, metric, ncore, electrons)
    expected = full_determinant_oracle(left, right, metric, ncore, electrons)
    np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)
    reverse = frozen_core_overlaps(right, left, metric.conj().T, ncore, electrons)
    np.testing.assert_allclose(actual, reverse.conj().T, atol=1e-12, rtol=1e-12)


def test_orthonormal_roots_and_root_phase_covariance():
    rng = np.random.default_rng(18)
    roots = np.linalg.qr(rng.normal(size=(24, 3)))[0].T.reshape(3, 6, 4)
    overlap = frozen_core_overlaps(roots, roots, np.eye(6), 2, (2, 1))
    np.testing.assert_allclose(overlap, np.eye(3), atol=1e-12, rtol=0)
    phases = np.exp(1j * np.array([0.2, -0.4, 0.7]))
    changed = frozen_core_overlaps(roots * phases[:, None, None], roots, np.eye(6), 2, (2, 1))
    np.testing.assert_allclose(changed, np.diag(phases.conj()), atol=1e-12, rtol=0)


def test_core_only_slater_overlap_includes_both_spins():
    metric = np.diag([0.8, 0.9, 1.0])
    ci = np.ones((1, 1, 1))
    actual = frozen_core_overlaps(ci, ci, metric, 2, (0, 0))
    np.testing.assert_allclose(actual, [[(0.8 * 0.9) ** 2]], atol=1e-14, rtol=0)


def test_singular_core_and_wrong_ci_shapes_reject():
    ci = np.ones((1, 2, 2))
    with pytest.raises(ValueError, match="Core cross block"):
        frozen_core_overlaps(ci, ci, np.diag([0.0, 1.0, 1.0]), 1, (1, 1))
    with pytest.raises(ValueError, match="determinant order"):
        frozen_core_overlaps(ci[:, :1], ci, np.eye(3), 1, (1, 1))
