"""Check the structure probe against actual unchanged Numerov solve operands."""

import numpy as np
import pytest

from benchmarks.co_exterior_sparsity import inspect_model, support
from projects.ukrmol_co import numerov_typed_coupled as reference
from projects.ukrmol_co.native_boundary import surface_matrix


def test_probe_matches_actual_coupled_recurrence(monkeypatch):
    amplitudes = np.array([[1.0, 0.5, 0.2, 0.3], [0.2, 0.4, 0.7, 0.6]])
    model = dict(
        energies=np.array([0.2]),
        ground=0.0,
        poles=np.array([1.0, 2.0]),
        amplitudes=amplitudes,
        thresholds=np.array([0.0, 0.0, 0.5, 0.5]),
        angular=np.zeros(4, dtype=int),
        c2=np.array([[1, 0.1, 0, 0], [0.1, 2, 0, 0], [0, 0, 3, 0.2], [0, 0, 0.2, 4]]),
        c3=np.zeros((4, 4)),
    )
    observed = inspect_model(model, 0)
    original = reference.refined_solve
    calls = []

    class EndOfProbe(Exception):
        pass

    def capture(a, b, stats):
        if len(calls) == 9:
            raise EndOfProbe
        x = original(a, b, stats)
        calls.append((a.copy(), b.copy(), x.copy()))
        return x

    monkeypatch.setattr(reference, "refined_solve", capture)
    q = surface_matrix(model["poles"], amplitudes, 0.1)
    with pytest.raises(EndOfProbe):
        reference.typed_coupled_k(
            q,
            model["thresholds"],
            model["angular"],
            model["c2"],
            model["c3"],
            0.2,
            18.0,
            18.0 + 5 * 0.0015625,
            0.0015625,
            8,
        )
    assert observed["initial_delta"] == support(calls[0][2].T)
    for i, row in enumerate(observed["intervals"]):
        potential, ratio = calls[1 + 2 * i : 3 + 2 * i]
        assert row["potential_a"] == support(potential[0])
        assert row["ratio_a"] == support(ratio[0])
        assert row["delta"] == support(ratio[2])
    assert observed["potential_component_sizes"] == [2, 2]
    for role in ["potential_factor", "ratio_factor"]:
        factor = observed["intervals"][0][role]
        assert factor["sparse_backward_error"] < 1e-12
        assert factor["relative_solution_difference"] < 1e-12


def test_support_retains_small_nonzero_coefficients():
    matrix = np.array([[1.0, 1e-100], [0.0, 1.0]])
    assert support(matrix)["nonzeros"] == 3
