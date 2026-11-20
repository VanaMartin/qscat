"""Small full spin-orbital determinants versus coherent core-sector contractions."""

from dataclasses import replace
from itertools import combinations
from math import comb

import numpy as np
import pytest
from scipy.linalg import block_diag

from projects.ukrmol_co.core_sector_full_overlap import (
    CoreInteriorExpansion,
    TargetProducts,
    core_sector_four_block_overlap,
)
from projects.ukrmol_co.core_sector_overlap import (
    CoreSectorCI,
    CoreSectorMetric,
    removal_overlaps,
    spin_lower,
)
from projects.ukrmol_co.four_block_overlap import InteriorExpansion, four_block_overlap


def occupations(n, k):
    return sorted(combinations(range(n), k), key=lambda s: sum(1 << i for i in s))


def sparse(state, root):
    n = state.ncore + state.nactive
    result = {}
    for i, a in enumerate(occupations(state.nactive, state.active_electrons[0])):
        for j, b in enumerate(occupations(state.nactive, state.active_electrons[1])):
            key = (
                *state.core_alpha,
                *(state.ncore + k for k in a),
                *(n + k for k in state.core_beta),
                *(n + state.ncore + k for k in b),
            )
            result[key] = state.coefficients[root, i, j]
    return result


def state(rng, counts, roots):
    shape = (roots, *(comb(2, k) for k in counts))
    c = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    return CoreSectorCI(1, 2, (0,), (0,), counts, c)


def expansion(rng, continuum=2):
    families = []
    for spin, counts, roots in [(0, (1, 1), 2), (1, (2, 0), 1)]:
        amplitudes = rng.normal(size=(roots, continuum)) + 1j * rng.normal(size=(roots, continuum))
        families.append(TargetProducts(state(rng, counts, roots), spin, amplitudes))
    return CoreInteriorExpansion(state(rng, (2, 1), 1), tuple(families))


def oracle_expansion(value):
    targets = []
    amplitudes = []
    continuum = value.products[0].amplitudes.shape[1]
    for family in value.products:
        for i in range(len(family.targets.coefficients)):
            targets.append(sparse(family.targets, i))
            a = np.zeros(2 * continuum, complex)
            a[family.continuum_spin * continuum : (family.continuum_spin + 1) * continuum] = (
                family.amplitudes[i]
            )
            amplitudes.append(a)
    return InteriorExpansion(sparse(value.bound, 0), tuple(targets), np.array(amplitudes), 6, 4)


def spin_metric(spatial, left_continuum, right_continuum):
    bb = spatial[:3, :3]
    bc = spatial[:3, 3:]
    cb = spatial[3:, :3]
    cc = spatial[3:, 3:]
    assert cc.shape == (left_continuum, right_continuum)
    return np.block(
        [[block_diag(bb, bb), block_diag(bc, bc)], [block_diag(cb, cb), block_diag(cc, cc)]]
    )


def compute(left, right, spatial):
    return core_sector_four_block_overlap(
        left,
        right,
        CoreSectorMetric(1, 2, spatial[:3, :3]),
        spatial[:3, 3:],
        spatial[3:, :3],
        spatial[3:, 3:],
    )


@pytest.mark.parametrize("singular", [False, True])
def test_all_blocks_vs_small_full_determinant_oracle_with_cross_spin(singular):
    rng = np.random.default_rng(817)
    left, right = expansion(rng, 2), expansion(rng, 1)
    spatial = (rng.normal(size=(5, 4)) + 1j * rng.normal(size=(5, 4))) / 3
    if singular:
        spatial[0, 0] = 0
    actual = compute(left, right, spatial)
    expected = four_block_overlap(
        oracle_expansion(left), oracle_expansion(right), spin_metric(spatial, 2, 1)
    )
    for key in expected:
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-10, atol=1e-10)
    reverse = compute(right, left, spatial.conj().T)
    for a, b in [("qq", "qq"), ("pq", "qp"), ("qp", "pq"), ("pp", "pp"), ("total", "total")]:
        np.testing.assert_allclose(actual[a], np.conj(reverse[b]), rtol=1e-10, atol=1e-10)


def test_cross_spin_weighted_removal_trace_matches_spin_lowering():
    rng = np.random.default_rng(319)
    left, right = state(rng, (1, 1), 2), state(rng, (2, 0), 1)
    spatial = np.eye(3) + (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) / 10
    metric = CoreSectorMetric(1, 2, spatial)
    removed = removal_overlaps(left, right, 1, metric, right_spin=0)
    assert np.max(abs(removed)) > 1e-3
    expected = metric.overlaps(left, spin_lower(right))
    np.testing.assert_allclose(
        np.einsum("mn,minj->ij", spatial, removed), expected, rtol=1e-10, atol=1e-10
    )
    reverse = removal_overlaps(
        right, left, 0, CoreSectorMetric(1, 2, spatial.conj().T), right_spin=1
    )
    np.testing.assert_allclose(
        removed, reverse.transpose(2, 3, 0, 1).conj(), rtol=1e-10, atol=1e-10
    )


def test_different_continuum_spins_have_nonzero_exchange_without_direct_target_term():
    rng = np.random.default_rng(723)
    left, right = expansion(rng), expansion(rng)
    left = replace(left, products=(left.products[0],))
    right = replace(right, products=(right.products[1],))
    spatial = (rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))) / 4
    actual = compute(left, right, spatial)
    expected = four_block_overlap(
        oracle_expansion(left), oracle_expansion(right), spin_metric(spatial, 2, 2)
    )
    np.testing.assert_allclose(actual["pp"], expected["pp"], rtol=1e-10, atol=1e-10)
    assert abs(actual["pp"]) > 1e-5


def test_positive_self_metric_gives_complete_state_norm_with_interference():
    rng = np.random.default_rng(98)
    value = expansion(rng)
    frame = (rng.normal(size=(7, 5)) + 1j * rng.normal(size=(7, 5))) / 4
    spatial = frame.conj().T @ frame
    actual = compute(value, value, spatial)
    old = oracle_expansion(value)
    p = {}
    for i, target in enumerate(old.targets):
        for j, a in enumerate(old.amplitudes[i]):
            for occupied, c in target.items():
                key = (*occupied, old.bound_orbitals + j)
                p[key] = p.get(key, 0) + a * c  # Four target electrons: even creation permutation.
    metric = spin_metric(spatial, 2, 2)
    full = dict(old.bound)
    for key, c in p.items():
        full[key] = full.get(key, 0) + c
    norm = sum(
        np.conj(a) * b * np.linalg.det(metric[np.ix_(i, j)])
        for i, a in full.items()
        for j, b in full.items()
    )
    np.testing.assert_allclose(actual["total"], norm, rtol=1e-10, atol=1e-10)
    assert actual["total"].real >= 0 and abs(actual["total"].imag) < 1e-10
    assert abs(actual["pq"]) > 1e-4
    assert abs(actual["total"] - actual["qq"] - actual["pp"]) > 1e-4


def test_target_basis_and_whole_state_phase_covariance_without_mutation():
    rng = np.random.default_rng(1001)
    left, right = expansion(rng), expansion(rng)
    spatial = (rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))) / 4
    expected = compute(left, right, spatial)
    phases = np.exp(1j * np.array([0.3, -0.7]))
    family = left.products[0]
    transformed = replace(
        family,
        targets=replace(
            family.targets, coefficients=phases[:, None, None] * family.targets.coefficients
        ),
        amplitudes=family.amplitudes / phases[:, None],
    )
    rotated = replace(left, products=(transformed, left.products[1]))
    actual = compute(rotated, right, spatial)
    for key in expected:
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-10, atol=1e-10)
    angle = np.exp(0.41j)
    scaled = replace(
        left,
        bound=replace(left.bound, coefficients=angle * left.bound.coefficients),
        products=tuple(replace(f, amplitudes=angle * f.amplitudes) for f in left.products),
    )
    actual = compute(scaled, right, spatial)
    for key in expected:
        np.testing.assert_allclose(
            actual[key], angle.conjugate() * expected[key], rtol=1e-10, atol=1e-10
        )


@pytest.mark.parametrize(
    "case", ["Qroots", "spin", "counts", "duplicate", "shape", "finite", "right_spin"]
)
def test_malformed_and_admission_guards(case):
    rng = np.random.default_rng(9)
    value = expansion(rng)
    spatial = np.eye(5)
    with pytest.raises(ValueError):
        if case == "Qroots":
            compute(replace(value, bound=state(rng, (2, 1), 2)), value, spatial)
        elif case == "spin":
            compute(
                replace(value, products=(replace(value.products[0], continuum_spin=2),)),
                value,
                spatial,
            )
        elif case == "counts":
            compute(
                replace(value, products=(replace(value.products[0], continuum_spin=1),)),
                value,
                spatial,
            )
        elif case == "duplicate":
            compute(replace(value, products=(value.products[0], value.products[0])), value, spatial)
        elif case == "shape":
            compute(value, value, np.ones((6, 5)))
        elif case == "finite":
            compute(value, value, np.full((5, 5), np.nan))
        else:
            removal_overlaps(
                value.products[0].targets,
                value.products[1].targets,
                1,
                CoreSectorMetric(1, 2, np.eye(3)),
                right_spin=2,
            )
