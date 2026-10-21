"""Physical CO ensemble-frame overlaps and retained-subspace diagnostics.

Uses the physical Hermitian inner product, not an ECS c-product or AO transport.
See docs/physics/co-ensemble-frame-overlaps.md for the finite dominance screen.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


def full_minor_overlaps(left, right, metric, ncore, electrons):
    """Independent full occupied minors, with integer-bitstring enumeration."""
    ncas = len(metric) - ncore
    determinants = []
    for count in electrons:
        occupied = np.array(
            [
                list(range(ncore)) + [ncore + i for i in range(ncas) if bits & (1 << i)]
                for bits in range(1 << ncas)
                if bits.bit_count() == count
            ],
            dtype=int,
        )
        minors = metric[occupied[:, None, :, None], occupied[None, :, None, :]]
        determinants.append(np.linalg.det(minors))
    alpha, beta = determinants
    transformed = np.stack([alpha @ ci @ beta.T for ci in right])
    return left.reshape(len(left), -1).conj() @ transformed.reshape(len(right), -1).T


def _assign_retained(scores):
    rows, columns = linear_sum_assignment(-scores[:5])
    details = []
    for row, column in zip(rows, columns, strict=True):
        best, second = np.argsort(scores[row])[::-1][:2]
        score = float(scores[row, column])
        margin = float(scores[row, best] - scores[row, second])
        details.append(
            dict(
                retained_root=int(row + 1),
                assigned_opposite_root=int(column + 1),
                assigned_squared_overlap=score,
                row_best_opposite_root=int(best + 1),
                best_to_second_squared_overlap_margin=margin,
                weak_assignment_trigger=score < 0.5,
                ambiguous_row_trigger=margin < 0.1,
                retained_to_extra_trigger=bool(column >= 5),
            )
        )
    return details


def diagnose_overlap(overlap):
    """Screen the retained five without rescaling any lost physical overlap."""
    overlap = np.asarray(overlap)
    if overlap.ndim != 2 or min(overlap.shape) < 5 or not np.isfinite(overlap).all():
        raise ValueError("Finite overlap with at least five roots on each side required")
    singular = np.linalg.svd(overlap, compute_uv=False)
    if singular[0] > 1 + 1e-9:
        raise ValueError("Physical root overlap is not contractive")
    retained = np.linalg.svd(overlap[:5, :5], compute_uv=False)
    scores = abs(overlap) ** 2
    return dict(
        root_overlap_shape=list(overlap.shape),
        root_buffer_singular_values=singular.tolist(),
        retained_squared_principal_overlaps=(retained**2).tolist(),
        minimum_retained_squared_principal_overlap=float(retained[-1] ** 2),
        retained_subspace_dominance_passed=bool(retained[-1] ** 2 >= 0.5),
        old_retained_projection_into_new_retained=scores[:5, :5].sum(axis=1).tolist(),
        new_retained_projection_into_old_retained=scores[:5, :5].sum(axis=0).tolist(),
        old_retained_projection_into_new_buffer=scores[:5].sum(axis=1).tolist(),
        new_retained_projection_into_old_buffer=scores[:, :5].sum(axis=0).tolist(),
        old_to_new_assignments=_assign_retained(scores),
        new_to_old_assignments=_assign_retained(scores.T),
    )


def decision_controls():
    """Analytic pass/reject cases frozen before any molecular overlaps are read."""
    permutation = np.eye(8)[:, [4, 0, 2, 1, 3, 5, 6, 7]]
    unitary = np.eye(8, dtype=complex)
    unitary[:5, :5] = np.exp(2j * np.pi * np.outer(np.arange(5), np.arange(5)) / 5) / np.sqrt(5)
    cases = dict(
        permutation=(permutation, True),
        unitary=(unitary, True),
        retained_buffer_exchange=(np.eye(8)[:, [0, 1, 2, 3, 5, 4, 6, 7]], False),
        contracted=(0.7 * np.eye(8), False),
    )
    results = {}
    for name, (overlap, expected) in cases.items():
        result = diagnose_overlap(overlap)
        if result["retained_subspace_dominance_passed"] is not expected:
            raise ValueError(f"Analytic subspace control failed: {name}")
        results[name] = result
    return results
