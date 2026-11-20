"""Small full spin-orbital determinants versus core-sector factorization."""

from itertools import combinations, product
from math import comb

import numpy as np
import pytest
from scipy.linalg import block_diag

from projects.ukrmol_co.core_sector_overlap import (
    CoreSectorCI,
    CoreSectorMetric,
    removal_overlaps,
    remove_spin_orbital,
    spin_lower,
)


def strings(n, k):
    return sorted(combinations(range(n), k), key=lambda s: sum(1 << i for i in s))


def full(state, root):
    n = state.ncore + state.nactive
    result = {}
    for i, a in enumerate(strings(state.nactive, state.active_electrons[0])):
        for j, b in enumerate(strings(state.nactive, state.active_electrons[1])):
            occupied = (
                *state.core_alpha,
                *(state.ncore + x for x in a),
                *(n + x for x in state.core_beta),
                *(n + state.ncore + x for x in b),
            )
            if state.coefficients[root, i, j] != 0:
                result[occupied] = state.coefficients[root, i, j]
    return result


def delete(state, orbital):
    result = {}
    for occupied, value in state.items():
        if orbital in occupied:
            index = occupied.index(orbital)
            result[occupied[:index] + occupied[index + 1 :]] = (-1) ** index * value
    return result


def inner(left, right, spatial):
    metric = block_diag(spatial, spatial)
    return sum(
        np.conj(a) * b * np.linalg.det(metric[np.ix_(i, j)])
        for i, a in left.items()
        for j, b in right.items()
    )


def random_state(rng, roots=2, ncore=1, nactive=3, counts=(2, 1), cores=None):
    cores = (tuple(range(ncore)),) * 2 if cores is None else cores
    shape = (roots, comb(nactive, counts[0]), comb(nactive, counts[1]))
    coefficients = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    return CoreSectorCI(ncore, nactive, *cores, counts, coefficients)


def test_target_and_every_core_active_spin_removal_against_full_determinants():
    rng = np.random.default_rng(215)
    left, right = random_state(rng), random_state(rng, roots=3)
    frame = [rng.normal(size=(7, 4)) + 1j * rng.normal(size=(7, 4)) for _ in range(2)]
    spatial = frame[0].conj().T @ frame[1] / 8
    metric = CoreSectorMetric(1, 3, spatial)
    expected = np.array(
        [[inner(full(left, i), full(right, j), spatial) for j in range(3)] for i in range(2)]
    )
    np.testing.assert_allclose(metric.overlaps(left, right), expected, rtol=1e-10, atol=1e-10)
    for spin, orbital in product(range(2), range(4)):
        removed = remove_spin_orbital(left, spin, orbital)
        for root in range(2):
            want = delete(full(left, root), spin * 4 + orbital)
            if removed is None:
                assert not any(want.values())
            else:
                actual = full(removed, root)
                assert set(actual) == set(want)
                np.testing.assert_allclose(
                    list(actual.values()), [want[k] for k in actual], rtol=0, atol=0
                )
        other = remove_spin_orbital(right, spin, orbital)
        if removed is not None and other is not None:
            value = np.array(
                [
                    [inner(full(removed, i), full(other, j), spatial) for j in range(3)]
                    for i in range(2)
                ]
            )
            np.testing.assert_allclose(
                metric.overlaps(removed, other), value, rtol=1e-10, atol=1e-10
            )


@pytest.mark.parametrize("spin", [0, 1])
def test_complete_removal_block_and_reverse_for_complex_cross_metric(spin):
    rng = np.random.default_rng(37)
    left, right = (
        random_state(rng, roots=2, nactive=2, counts=(1, 1)),
        random_state(rng, roots=1, nactive=2, counts=(1, 1)),
    )
    spatial = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    actual = removal_overlaps(left, right, spin, CoreSectorMetric(1, 2, spatial))
    for mu, nu, i, j in product(range(3), range(3), range(2), range(1)):
        value = inner(
            delete(full(left, i), spin * 3 + mu), delete(full(right, j), spin * 3 + nu), spatial
        )
        np.testing.assert_allclose(actual[mu, i, nu, j], value, rtol=1e-10, atol=1e-10)
    reverse = removal_overlaps(right, left, spin, CoreSectorMetric(1, 2, spatial.conj().T))
    np.testing.assert_allclose(actual, reverse.transpose(2, 3, 0, 1).conj(), rtol=1e-10, atol=1e-10)


def test_singular_core_and_nonzero_core_hole_active_hole_block():
    state = CoreSectorCI(1, 1, (0,), (), (1, 0), np.ones((1, 1, 1)))
    spatial = np.array([[0.0, 1.0], [1.0, 0.0]])
    metric = CoreSectorMetric(1, 1, spatial)
    np.testing.assert_allclose(metric.overlaps(state, state), [[-1]], rtol=0, atol=1e-14)
    core = remove_spin_orbital(state, 0, 0)
    active = remove_spin_orbital(state, 0, 1)
    np.testing.assert_allclose(metric.overlaps(core, active), [[-1]], rtol=0, atol=1e-14)
    assert remove_spin_orbital(state, 1, 0) is None


@pytest.mark.parametrize("ncore", [0, 1, 2])
def test_spin_lowering_norm_signs_and_singlet_orthogonality(ncore):
    core = tuple(range(ncore))
    state = CoreSectorCI(ncore, 2, core, core, (2, 0), np.ones((1, 1, 1)))
    lowered = spin_lower(state)
    expected = {}
    n = ncore + 2
    source = full(state, 0)
    for spatial in range(n):
        removed = delete(source, spatial)
        for occupied, value in removed.items():
            created = n + spatial
            if created in occupied:
                continue
            insertion = sum(i < created for i in occupied)
            key = tuple(sorted((*occupied, created)))
            expected[key] = expected.get(key, 0) + (-1) ** insertion * value
    actual = full(lowered, 0)
    assert set(actual) == set(expected)
    np.testing.assert_allclose(list(actual.values()), [expected[k] for k in actual], rtol=0, atol=0)
    metric = CoreSectorMetric(ncore, 2, np.eye(n))
    np.testing.assert_allclose(metric.overlaps(lowered, lowered), [[2]], rtol=0, atol=1e-14)
    singlet = np.zeros((1, 2, 2))
    singlet[0, 0, 1] = singlet[0, 1, 0] = 1 / np.sqrt(2)
    singlet = CoreSectorCI(ncore, 2, core, core, (1, 1), singlet)
    np.testing.assert_allclose(metric.overlaps(lowered, singlet), [[0]], rtol=0, atol=1e-14)


def test_spin_changing_removal_between_target_projections_is_not_dropped():
    triplet = CoreSectorCI(0, 2, (), (), (2, 0), np.ones((1, 1, 1)))
    singlet = np.zeros((1, 2, 2))
    singlet[0, 0, 1] = singlet[0, 1, 0] = 1 / np.sqrt(2)
    singlet = CoreSectorCI(0, 2, (), (), (1, 1), singlet)
    metric = CoreSectorMetric(0, 2, np.eye(2))
    np.testing.assert_allclose(metric.overlaps(triplet, singlet), [[0]], rtol=0, atol=0)
    value = metric.overlaps(remove_spin_orbital(triplet, 0, 0), remove_spin_orbital(singlet, 1, 0))
    np.testing.assert_allclose(value, [[-1 / np.sqrt(2)]], rtol=0, atol=1e-14)


@pytest.mark.parametrize("spin", [0, 1])
def test_metric_weighted_removal_trace_is_spin_electron_count_times_overlap(spin):
    rng = np.random.default_rng(827)
    left = random_state(rng, roots=2, nactive=2, counts=(1, 1))
    right = random_state(rng, roots=3, nactive=2, counts=(1, 1))
    spatial = np.eye(3) + (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) / 10
    metric = CoreSectorMetric(1, 2, spatial)
    removed = removal_overlaps(left, right, spin, metric)
    contracted = np.einsum("mn,minj->ij", spatial, removed)
    np.testing.assert_allclose(contracted, 2 * metric.overlaps(left, right), rtol=1e-10, atol=1e-10)


def test_metric_is_frozen_and_inputs_are_not_modified():
    rng = np.random.default_rng(17)
    state = random_state(rng, nactive=2, counts=(1, 1))
    coefficients = state.coefficients.copy()
    spatial = np.eye(3)
    metric = CoreSectorMetric(1, 2, spatial)
    want = metric.overlaps(state, state)
    spatial[:] = 0
    remove_spin_orbital(state, 1, 1)
    spin_lower(state)
    np.testing.assert_allclose(metric.overlaps(state, state), want, rtol=0, atol=0)
    np.testing.assert_array_equal(state.coefficients, coefficients)


def test_root_batches_and_bounded_cache_eviction_preserve_overlaps():
    rng = np.random.default_rng(111)
    left = random_state(rng, roots=9, nactive=2, counts=(1, 1))
    right = random_state(rng, roots=11, nactive=2, counts=(1, 1))
    spatial = np.eye(3) + rng.normal(size=(3, 3)) / 10
    metric = CoreSectorMetric(1, 2, spatial, cache_bytes=64)
    expected = np.array(
        [[inner(full(left, i), full(right, j), spatial) for j in range(11)] for i in range(9)]
    )
    np.testing.assert_allclose(metric.overlaps(left, right), expected, rtol=1e-10, atol=1e-10)
    actual = removal_overlaps(left, right, 0, metric)
    reference = removal_overlaps(left, right, 0, CoreSectorMetric(1, 2, spatial))
    np.testing.assert_allclose(actual, reference, rtol=1e-12, atol=1e-12)
    assert sum(value.nbytes for value in metric._cache.values()) <= 64


def test_vacuum_and_pauli_blocked_spin_lowering_are_exact():
    state = CoreSectorCI(0, 0, (), (), (0, 0), np.array([[[2j]]]))
    metric = CoreSectorMetric(0, 0, np.empty((0, 0)))
    np.testing.assert_allclose(metric.overlaps(state, state), [[4]], rtol=0, atol=0)
    assert spin_lower(state) is None


@pytest.mark.parametrize(
    "case", ["core", "shape", "roots", "nan", "pairs", "cache", "memory", "spin"]
)
def test_malformed_and_resource_caps(case):
    state = CoreSectorCI(1, 2, (0,), (0,), (1, 1), np.ones((1, 2, 2)))
    metric = CoreSectorMetric(1, 2, np.eye(3))
    with pytest.raises(ValueError):
        if case == "core":
            metric.overlaps(CoreSectorCI(1, 2, (0, 0), (0,), (1, 1), np.ones((1, 2, 2))), state)
        elif case == "shape":
            metric.overlaps(CoreSectorCI(1, 2, (0,), (0,), (1, 1), np.ones((1, 1, 2))), state)
        elif case == "roots":
            metric.overlaps(CoreSectorCI(1, 2, (0,), (0,), (1, 1), np.ones((65, 2, 2))), state)
        elif case == "nan":
            CoreSectorMetric(1, 2, np.full((3, 3), np.nan))
        elif case == "pairs":
            CoreSectorMetric(1, 2, np.eye(3), max_pairs=1).overlaps(state, state)
        elif case == "cache":
            CoreSectorMetric(1, 2, np.eye(3), cache_bytes=1).overlaps(state, state)
        elif case == "memory":
            removal_overlaps(state, state, 0, metric, max_left_bytes=1)
        elif case == "spin":
            remove_spin_orbital(state, 2, 0)
