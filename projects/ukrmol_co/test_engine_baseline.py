"""A same-root control must reject drift without claiming spectral convergence."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from projects.ukrmol_co.engine_baseline import audit_baseline, label_engine_reproduction
from projects.ukrmol_co.sparse_comparison import compare_phases


def settings():
    return json.loads(
        Path(__file__).with_name("cas12-audited-engine-baseline-contract.json").read_text()
    )


def phases(position_ev=2.457, branch=0):
    contract = settings()
    requested = 0.1 + 0.05 * np.arange(99)
    energy = requested * contract["native_rydberg_per_requested_ev"] / 2
    center = position_ev / contract["hartree_to_ev"]
    width = 0.5 / contract["hartree_to_ev"]
    s = (center - energy + 0.5j * width) / (center - energy - 0.5j * width)
    np.testing.assert_allclose(abs(s), 1, atol=3e-15, rtol=0)
    return np.column_stack((requested, np.angle(s) / 2 + branch * np.pi))


@pytest.mark.parametrize("shift,expected", [(0, True), (0.001, False)])
def test_engine_reproduction_is_tighter_than_scientific_pair_budget(shift, expected):
    row = compare_phases(phases(), phases(position_ev=2.457 + shift, branch=2), settings())
    result = label_engine_reproduction(
        {
            "root_pair_stability_verdict": row["root_pair_stability_verdict"],
            "spectral_truncation_verdict": None,
            "sectors": {name: copy.deepcopy(row) for name in ("B1", "B2")},
        }
    )
    assert result["matched_engine_baseline_verdict"] is expected
    assert result["root_pair_stability_verdict"] is None
    assert result["spectral_truncation_verdict"] is None
    assert all(r["root_pair_stability_verdict"] is None for r in result["sectors"].values())


def test_invalid_fit_stays_unqualified_and_changed_root_count_is_rejected(tmp_path):
    flat = phases()
    flat[:, 1] = 0
    row = compare_phases(flat, flat, settings())
    result = label_engine_reproduction(
        {
            "root_pair_stability_verdict": row["root_pair_stability_verdict"],
            "sectors": {"B1": row},
        }
    )
    assert result["matched_engine_baseline_verdict"] is None
    with pytest.raises(ValueError, match="2048 roots on both sides"):
        audit_baseline(tmp_path, tmp_path, settings() | {"candidate_roots": 4096})
