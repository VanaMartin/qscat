"""Independent antisymmetrized-coordinate and residual-bound controls."""

from itertools import combinations, permutations

import numpy as np
import pytest

from projects.ukrmol_co.ci_overlap import frozen_core_overlaps
from projects.ukrmol_co.interior_state_overlap import (
    bound_sector_projection,
    projected_overlap_interval,
    sparse_state_overlap,
)


def leibniz(matrix):
    total = 0j
    for order in permutations(range(len(matrix))):
        inversions = sum(
            order[i] > order[j] for i in range(len(order)) for j in range(i + 1, len(order))
        )
        total += (-1) ** inversions * np.prod([matrix[i, j] for i, j in enumerate(order)])
    return total


def exterior_coordinates(frame, occupations, coefficients):
    return np.array(
        [
            sum(
                c * leibniz(frame[np.ix_(rows, cols)])
                for cols, c in zip(occupations, coefficients, strict=True)
            )
            for rows in combinations(range(len(frame)), occupations.shape[1])
        ]
    )


@pytest.mark.parametrize("electrons", [0, 1, 2, 3])
@pytest.mark.parametrize("tile_size", [1, 2, 7])
def test_sparse_minors_against_independent_exterior_coordinates(electrons, tile_size):
    rng = np.random.default_rng(52)
    frames = [
        np.linalg.qr(rng.normal(size=(7, n)) + 1j * rng.normal(size=(7, n)))[0] for n in [4, 5]
    ]
    occupations = [
        np.array(list(combinations(range(n), electrons))).reshape(-1, electrons)
        if electrons
        else np.empty((1, 0), int)
        for n in [4, 5]
    ]
    coefficients = [rng.normal(size=len(o)) + 1j * rng.normal(size=len(o)) for o in occupations]
    left, right = [
        exterior_coordinates(f, o, c)
        for f, o, c in zip(frames, occupations, coefficients, strict=True)
    ]
    metric = frames[0].conj().T @ frames[1]
    actual = sparse_state_overlap(
        occupations[0],
        coefficients[0],
        occupations[1],
        coefficients[1],
        metric,
        tile_size=tile_size,
    )
    np.testing.assert_allclose(actual, np.vdot(left, right), rtol=1e-12, atol=1e-12)
    reverse = sparse_state_overlap(
        occupations[1],
        coefficients[1],
        occupations[0],
        coefficients[0],
        metric.conj().T,
        tile_size=tile_size,
    )
    np.testing.assert_allclose(actual, reverse.conjugate(), rtol=1e-12, atol=1e-12)


def test_self_norm_core_cross_terms_and_bound_continuum_exchange():
    occupations = np.array([[0, 2], [0, 3], [1, 2]])
    coefficients = np.array([1, 2j, -3]) / np.sqrt(14)
    own = sparse_state_overlap(occupations, coefficients, occupations, coefficients, np.eye(4))
    np.testing.assert_allclose(own, 1, rtol=0, atol=1e-14)
    metric = np.array([[0.8, 0, 0.2, 0], [0, 0.7, 0, 0.3], [0.4, 0, 0.6, 0], [0, 0.1, 0, 0.5]])
    singleton = np.array([[0, 2]])
    value = sparse_state_overlap(singleton, [1], singleton, [1], metric)
    np.testing.assert_allclose(value, 0.8 * 0.6 - 0.2 * 0.4, rtol=0, atol=1e-14)
    # Exchanging a nominal bound and continuum orbital is visible to the metric.
    exchange = np.eye(4)[:, [0, 3, 2, 1]]
    changed = sparse_state_overlap(occupations, coefficients, occupations, coefficients, exchange)
    assert abs(changed) < 0.9


def test_phase_and_determinant_permutation_covariance():
    occupations = np.array([[0, 1], [0, 2], [1, 2]])
    vector = np.array([1, 2j, 3]) / np.sqrt(14)
    phase = np.exp(0.7j)
    value = sparse_state_overlap(
        occupations[::-1], phase * vector[::-1], occupations, vector, np.eye(3)
    )
    np.testing.assert_allclose(value, phase.conjugate(), rtol=0, atol=1e-14)


def test_orbital_permutation_requires_occupied_parity():
    metric = np.eye(3)[:, [1, 0, 2]]
    occupations = np.array([[0, 1]])
    value = sparse_state_overlap(occupations, [1], occupations, [-1], metric)
    np.testing.assert_allclose(value, 1, rtol=0, atol=1e-14)


def test_unitary_orbital_change_and_spin_adapted_csfs():
    rng = np.random.default_rng(71)
    occupied = np.array(list(combinations(range(4), 2)))
    rotation = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))[0]
    induced = np.array([[leibniz(rotation[np.ix_(a, b)]) for b in occupied] for a in occupied])
    vector = rng.normal(size=6) + 1j * rng.normal(size=6)
    changed = induced.conj().T @ vector
    value = sparse_state_overlap(occupied, vector, occupied, changed, rotation)
    np.testing.assert_allclose(value, np.vdot(vector, vector), rtol=1e-12, atol=1e-12)
    # Closed core plus two open orbitals: singlet and M=0 triplet.
    csfs = np.array([[0, 1, 2, 5], [0, 1, 3, 4]])
    singlet, triplet = np.array([1, -1]) / np.sqrt(2), np.ones(2) / np.sqrt(2)
    spatial = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    spin_metric = np.kron(spatial, np.eye(2))
    forbidden = sparse_state_overlap(csfs, singlet, csfs, triplet, spin_metric)
    np.testing.assert_allclose(forbidden, 0, rtol=0, atol=1e-13)


def test_fifteen_electron_bound_sector_against_full_spin_orbital_minors():
    rng = np.random.default_rng(90)
    reference = np.eye(26)[:, :13]
    moved = np.linalg.qr(reference + 0.02 * rng.normal(size=(26, 13)))[0]
    metric = reference.T @ moved
    active = [
        sorted(combinations(range(11), n), key=lambda c: sum(1 << i for i in c)) for n in (6, 5)
    ]
    indices = rng.choice(462 * 462, size=8, replace=False)
    occupied = []
    for index in indices:
        a, b = divmod(int(index), 462)
        occupied.append(
            [0, 1, *[i + 2 for i in active[0][a]], 13, 14, *[i + 15 for i in active[1][b]]]
        )
    vectors = [rng.normal(size=8) + 1j * rng.normal(size=8) for _ in range(2)]
    tensors = []
    for vector in vectors:
        tensor = np.zeros(462 * 462, complex)
        tensor[indices] = vector
        tensors.append(tensor.reshape(1, 462, 462))
    full_metric = np.kron(np.eye(2), metric)
    expected = sparse_state_overlap(occupied, vectors[0], occupied, vectors[1], full_metric)
    actual = frozen_core_overlaps(tensors[0], tensors[1], metric, 2, (6, 5))[0, 0]
    np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-12)
    projected = bound_sector_projection(tensors[1], metric, 2, (6, 5))
    np.testing.assert_allclose(np.vdot(tensors[0], projected), expected, rtol=1e-11, atol=1e-12)
    assert np.linalg.norm(projected) <= np.linalg.norm(tensors[1]) + 1e-12


def test_three_projection_sizes_enclose_exact_overlap_without_renormalizing():
    rng = np.random.default_rng(83)
    occupations = np.array(list(combinations(range(5), 2)))
    frames = [np.linalg.qr(rng.normal(size=(7, 5)))[0] for _ in range(2)]
    vectors = [rng.normal(size=10) + 1j * rng.normal(size=10) for _ in range(2)]
    vectors = [v / np.linalg.norm(v) for v in vectors]
    metric = frames[0].T @ frames[1]
    full = sparse_state_overlap(occupations, vectors[0], occupations, vectors[1], metric)
    widths = []
    for size in [2, 6, 10]:
        partial = sparse_state_overlap(
            occupations[:size], vectors[0][:size], occupations[:size], vectors[1][:size], metric
        )
        capture = [float(np.vdot(v[:size], v[:size]).real) for v in vectors]
        # Exact unit captures avoid treating tiny summation drift as negative loss.
        if size == 10:
            capture = [1.0, 1.0]
        lower, upper = projected_overlap_interval(partial, *capture, absolute_error=1e-13)
        assert lower <= abs(full) ** 2 <= upper
        widths.append(upper - lower)
    assert widths[0] >= widths[1] >= widths[2]
    assert widths[-1] < 1e-12


def test_projection_bound_self_null_ambiguous_and_omitted_character():
    np.testing.assert_allclose(projected_overlap_interval(1, 1, 1), [1, 1], rtol=0, atol=0)
    np.testing.assert_allclose(projected_overlap_interval(0, 1, 1), [0, 0], rtol=0, atol=0)
    np.testing.assert_allclose(projected_overlap_interval(0, 0, 0), [0, 1], rtol=0, atol=0)
    lower, upper = projected_overlap_interval(0.5, 0.5, 0.5)
    assert lower == 0 and upper == 1  # A common surface projection need not fix the state.
    np.testing.assert_allclose(
        projected_overlap_interval(1 / np.sqrt(2), 1, 1), [0.5, 0.5], rtol=1e-14, atol=1e-14
    )


@pytest.mark.parametrize(
    "overlap,a,b,error",
    [
        (1.1, 1, 1, 0),
        (0.3, 0.1, 0.1, 0),
        (0, -0.1, 1, 0),
        (0, 1.01, 1, 0),
        (np.nan, 1, 1, 0),
        (0, 1, 1, -1),
        (0, 1, 1, np.inf),
    ],
)
def test_invalid_capture_or_partial_metric_rejects(overlap, a, b, error):
    with pytest.raises(ValueError):
        projected_overlap_interval(overlap, a, b, absolute_error=error)


@pytest.mark.parametrize(
    "occupations,coefficients,metric",
    [
        ([[1, 0]], [1], np.eye(2)),
        ([[0, 0]], [1], np.eye(2)),
        ([[0, 1], [0, 1]], [1, 1], np.eye(2)),
        ([[0, 2]], [1], np.eye(2)),
        ([[0.0, 1.0]], [1], np.eye(2)),
        ([[0, 1]], [np.nan], np.eye(2)),
        ([[0, 1]], [1, 2], np.eye(2)),
        ([[0, 1]], [1], np.full((2, 2), np.inf)),
    ],
)
def test_malformed_sparse_inputs_reject(occupations, coefficients, metric):
    with pytest.raises(ValueError):
        sparse_state_overlap(occupations, coefficients, [[0, 1]], [1], metric)


def test_known_common_projection_tightens_the_omitted_norm_bound():
    interval = projected_overlap_interval(
        0.9,
        0.9,
        0.9,
        absolute_error=1e-15,
        left_complement_right_norm=0,
        right_complement_left_norm=0,
    )
    np.testing.assert_allclose(interval, [0.64, 1], rtol=1e-14, atol=1e-14)
    for phase in [-1, 1]:
        a = np.array([np.sqrt(0.9), np.sqrt(0.1)])
        b = np.array([np.sqrt(0.9), phase * np.sqrt(0.1)])
        value = abs(np.vdot(a, b)) ** 2
        assert interval[0] - 1e-14 <= value <= interval[1] + 1e-14


def test_complete_frozen_core_projection_matches_explicit_identity_roots():
    metric = np.array([[0.9, 0.2, 0], [0.1, 0.8, 0.1], [0, 0.1, 0.7]])
    right = np.array([[[1.0], [2j]]]) / np.sqrt(5)
    projected = bound_sector_projection(right, metric, 1, (1, 0))
    explicit = frozen_core_overlaps(np.eye(2).reshape(2, 2, 1), right, metric, 1, (1, 0))
    np.testing.assert_allclose(projected.reshape(1, -1).T, explicit, rtol=0, atol=1e-14)


def test_full_minor_projection_allows_singular_cross_core_and_rejects_oversize():
    right = np.ones((1, 1, 1))
    projected = bound_sector_projection(right, np.diag([0.0, 1.0]), 1, (0, 0))
    np.testing.assert_allclose(projected, 0, rtol=0, atol=0)
    with pytest.raises(ValueError, match="allocation cap"):
        bound_sector_projection(right, np.eye(32), 0, (16, 16))


@pytest.mark.parametrize("size", [1, 3, 5])
def test_distinct_projection_spaces_enclose_an_independent_full_vector_overlap(size):
    rng = np.random.default_rng(45)
    frames = [np.linalg.qr(rng.normal(size=(6, size)))[0] for _ in range(2)]
    states = [rng.normal(size=6) + 1j * rng.normal(size=6) for _ in range(2)]
    states = [state / np.linalg.norm(state) for state in states]
    parts = [frame @ frame.T @ state for frame, state in zip(frames, states, strict=True)]
    captures = [float(np.vdot(p, p).real) for p in parts]
    complements = [
        np.linalg.norm(parts[1] - frames[0] @ frames[0].T @ parts[1]),
        np.linalg.norm(parts[0] - frames[1] @ frames[1].T @ parts[0]),
    ]
    lower, upper = projected_overlap_interval(
        np.vdot(*parts),
        *captures,
        absolute_error=1e-13,
        left_complement_right_norm=complements[0],
        right_complement_left_norm=complements[1],
    )
    exact = abs(np.vdot(*states)) ** 2
    assert lower <= exact <= upper


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(left_complement_right_norm=0),
        dict(left_complement_right_norm=2, right_complement_left_norm=0),
        dict(left_complement_right_norm=0, right_complement_left_norm=np.nan),
    ],
)
def test_unjustified_complementary_norm_arguments_reject(kwargs):
    with pytest.raises(ValueError):
        projected_overlap_interval(0, 0.9, 0.9, **kwargs)


def test_work_caps_and_particle_number_are_checked_before_minors():
    occupations = np.array([[0], [1]])
    for kwargs in [dict(max_pairs=3), dict(tile_size=0), dict(tile_size=257), dict(max_pairs=True)]:
        with pytest.raises(ValueError):
            sparse_state_overlap(occupations, [1, 0], occupations, [1, 0], np.eye(2), **kwargs)
    with pytest.raises(ValueError, match="electron"):
        sparse_state_overlap([[0]], [1], [[0, 1]], [1], np.eye(2))
