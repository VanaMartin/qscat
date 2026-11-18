"""Full exterior coordinates versus the factorized removal-state expression."""

from itertools import combinations, permutations

import numpy as np
import pytest

from projects.ukrmol_co.four_block_overlap import (
    InteriorExpansion,
    coefficient_removal,
    four_block_overlap,
)


def determinant(matrix):
    result = 0j
    for order in permutations(range(len(matrix))):
        swaps = sum(
            order[i] > order[j] for i in range(len(order)) for j in range(i + 1, len(order))
        )
        result += (-1) ** swaps * np.prod([matrix[i, j] for i, j in enumerate(order)])
    return result


def explicit_sectors(state):
    p = {}
    for i, target in enumerate(state.targets):
        for j, amplitude in enumerate(state.amplitudes[i]):
            for occupied, coefficient in target.items():
                key = (*occupied, state.bound_orbitals + j)
                p[key] = p.get(key, 0j) + (-1) ** len(occupied) * amplitude * coefficient
    return dict(q=state.bound, p=p)


def ambient_coordinates(state, frame, electrons):
    return np.array(
        [
            sum(c * determinant(frame[np.ix_(rows, columns)]) for columns, c in state.items())
            for rows in combinations(range(len(frame)), electrons)
        ]
    )


def random_state(rng, bound, continuum, electrons):
    def expansion(count):
        return {c: complex(*rng.normal(size=2)) for c in combinations(range(bound), count)}

    return InteriorExpansion(
        expansion(electrons + 1),
        (expansion(electrons), expansion(electrons)),
        rng.normal(size=(2, continuum)) + 1j * rng.normal(size=(2, continuum)),
        bound,
        electrons,
    )


@pytest.mark.parametrize("electrons", [0, 1, 2, 3])
@pytest.mark.parametrize("orthonormal", [False, True])
def test_all_blocks_against_independent_full_exterior_coordinates(electrons, orthonormal):
    rng = np.random.default_rng(108)
    frames = [rng.normal(size=(7, n)) + 1j * rng.normal(size=(7, n)) for n in [5, 6]]
    if orthonormal:
        frames = [np.linalg.qr(frame)[0] for frame in frames]
    else:
        frames = [frame / 3 for frame in frames]
    left, right = random_state(rng, 3, 2, electrons), random_state(rng, 4, 2, electrons)
    actual = four_block_overlap(left, right, frames[0].conj().T @ frames[1])
    coordinates = [
        {
            kind: ambient_coordinates(state, frame, electrons + 1)
            for kind, state in explicit_sectors(expansion).items()
        }
        for expansion, frame in zip([left, right], frames, strict=True)
    ]
    for bra, ket in [("q", "q"), ("p", "q"), ("q", "p"), ("p", "p")]:
        expected = np.vdot(coordinates[0][bra], coordinates[1][ket])
        np.testing.assert_allclose(actual[bra + ket], expected, rtol=1e-11, atol=1e-11)
    reverse = four_block_overlap(right, left, frames[1].conj().T @ frames[0])
    np.testing.assert_allclose(
        actual["total"], reverse["total"].conjugate(), rtol=1e-11, atol=1e-11
    )
    total = np.vdot(sum(coordinates[0].values()), sum(coordinates[1].values()))
    np.testing.assert_allclose(actual["total"], total, rtol=1e-11, atol=1e-11)


def test_singular_bound_metric_has_nonzero_exchange_and_no_inverse():
    state = InteriorExpansion({}, ({(0, 1): 1},), np.ones((1, 1)), 2, 2)
    metric = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    result = four_block_overlap(state, state, metric)
    np.testing.assert_allclose(result["pp"], -1, rtol=0, atol=1e-14)


def test_core_removal_is_required_not_an_optional_frozen_core_drop():
    state = InteriorExpansion({}, ({(0, 1): 1},), np.ones((1, 1)), 2, 2)
    metric = np.array([[1, 0, 0.2], [0, 1, 0.3], [0.4, 0.5, 0.9]])
    result = four_block_overlap(state, state, metric)
    np.testing.assert_allclose(result["pp"], 0.67, rtol=0, atol=1e-14)
    assert abs(result["pp"] - 0.9) > 0.1
    assert coefficient_removal({(0, 1): 2j}, 1) == {(0,): -2j}


def test_spin_adapted_core_plus_singlet_triplet_remain_orthogonal():
    rng = np.random.default_rng(119)
    spatial = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    metric = np.kron(spatial, np.eye(2))
    targets = ({(0, 1, 2): 1}, {(0, 1, 3): 1})
    singlet = InteriorExpansion({}, targets, np.array([[0, -1], [1, 0]]) / np.sqrt(2), 4, 3)
    triplet = InteriorExpansion({}, targets, np.array([[0, -1], [-1, 0]]) / np.sqrt(2), 4, 3)
    result = four_block_overlap(singlet, triplet, metric)
    np.testing.assert_allclose(result["total"], 0, rtol=0, atol=1e-12)


def test_self_norm_is_whole_state_norm_not_sector_renormalization():
    rng = np.random.default_rng(32)
    state = random_state(rng, 3, 2, 1)
    frame = rng.normal(size=(6, 5)) / 2
    result = four_block_overlap(state, state, frame.T @ frame)
    parts = explicit_sectors(state)
    coordinates = [ambient_coordinates(part, frame, 2) for part in parts.values()]
    expected = np.linalg.norm(sum(coordinates)) ** 2
    np.testing.assert_allclose(result["total"], expected, rtol=1e-12, atol=1e-12)
    assert result["total"].real > 0 and abs(result["pq"]) > 1e-3


def test_work_budget_and_malformed_inputs_stop_before_unbounded_contraction():
    state = InteriorExpansion({(0, 1): 1}, ({(0,): 1},), np.ones((1, 1)), 2, 1)
    with pytest.raises(ValueError, match="pair budget"):
        four_block_overlap(state, state, np.eye(3), max_pairs=1)
    for metric in [np.eye(2), np.full((3, 3), np.nan)]:
        with pytest.raises(ValueError):
            four_block_overlap(state, state, metric)
    bad = InteriorExpansion({}, ({(1, 0): 1},), np.ones((1, 1)), 2, 2)
    with pytest.raises(ValueError, match="canonical"):
        four_block_overlap(bad, bad, np.eye(3))


def test_target_phase_and_order_are_coordinate_choices_not_new_states():
    rng = np.random.default_rng(147)
    left, right = random_state(rng, 3, 2, 2), random_state(rng, 3, 2, 2)
    metric = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
    expected = four_block_overlap(left, right, metric)
    phases = np.exp(1j * np.array([0.3, -0.7]))
    targets = tuple(
        {key: phase * value for key, value in target.items()}
        for target, phase in zip(left.targets, phases, strict=True)
    )
    changed = InteriorExpansion(
        left.bound, targets[::-1], (left.amplitudes / phases[:, None])[::-1], 3, 2
    )
    actual = four_block_overlap(changed, right, metric)
    for key in expected:
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-11, atol=1e-11)


def test_continuum_permutation_and_total_complex_phase_covariance():
    rng = np.random.default_rng(177)
    left, right = random_state(rng, 3, 2, 1), random_state(rng, 3, 2, 1)
    metric = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
    expected = four_block_overlap(left, right, metric)["total"]
    phase = np.exp(0.6j)
    changed = InteriorExpansion(
        {key: phase * value for key, value in left.bound.items()},
        left.targets,
        phase * left.amplitudes[:, ::-1],
        3,
        1,
    )
    actual = four_block_overlap(changed, right, metric[[0, 1, 2, 4, 3], :])["total"]
    np.testing.assert_allclose(actual, phase.conjugate() * expected, rtol=1e-11, atol=1e-11)
