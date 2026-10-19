"""Deterministic stress selection, historical coverage and independent toy boundary."""

import json

import numpy as np
import pytest

from projects.ukrmol_co import numerov_interior_selection as selection
from projects.ukrmol_co.numerov_increment_study import exact_free_k
from projects.ukrmol_co.numerov_interior_analytic import boundary


def test_selection_keeps_qualified_peaks_and_ties_without_candidate_error_selection(
    tmp_path, monkeypatch
):
    cases = [
        f"{g}/nine/{s}" for g in ["compressed", "equilibrium", "stretched"] for s in ["B1", "B2"]
    ]
    history = tmp_path / "history.json"
    history.write_text(
        json.dumps(
            {
                "probes": [],
                "remainder": [
                    dict(case=c, energy_index=i, rows=[{}])
                    for c in cases
                    for i in [*range(20), 48, 98]
                ],
            }
        )
    )
    prior = tmp_path / "qualified.json"
    passed = [dict(case=c, energy_index=i, verdict=True) for c in cases for i in [0, 98]]
    passed += [
        dict(case=c, energy_index=18, verdict=True) for c in cases if c.startswith("stretched")
    ]
    prior.write_text(json.dumps(dict(experiment_verdict=True, cases=passed)))

    def read(path):
        k = np.zeros((99, 6, 6))
        k[0, 0, 0] = k[98, 0, 0] = 1000  # excluded corners must not steal an interior peak
        index = 87 if "compressed" in str(path) else (53 if "equilibrium" in str(path) else 18)
        k[index, 0, 0] = 200
        if index == 87:
            k[88, 0, 0] = 200  # first tie remains index87
        return np.arange(99, dtype=float) + 1, k, {"same_header": True}

    monkeypatch.setattr(selection, "read_k", read)
    result = selection.select(tmp_path, history, prior, cases)
    assert [x["energy_index"] for x in result["peaks"]] == [87, 87, 53, 53, 18, 18]
    assert len(result["new_identities"]) == 10
    assert len(result["reused_qualified_peaks"]) == 2
    assert result["new_unique_identities_if_all_computed"] == 4
    assert not any(
        x["case"].startswith("stretched") and x["role"] == "peak" for x in result["new_identities"]
    )


@pytest.mark.parametrize("energy,target", [(0.327075, 320.0), (0.202125, -55.0)])
@pytest.mark.parametrize("coupled", [False, True])
def test_high_momentum_boundary_matches_independent_continuous_oracle(energy, target, coupled):
    q, thresholds = boundary(energy, target, coupled)
    exact = np.array(exact_free_k(q, thresholds, energy, 18.0, 100).tolist(), float)
    eigenvalues = np.linalg.eigvalsh(exact)
    assert np.min(abs(eigenvalues - target)) < 1e-8
    if coupled:
        assert np.min(abs(eigenvalues - 0.3)) < 1e-12
