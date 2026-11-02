"""Reader contract and fixed scope controls; no campaign files or integrations."""

import json
from pathlib import Path

from projects.ukrmol_co.verify_stretched_contour_study import CONTRACT, CONTRACT_SHA, checks, sha


def test_frozen_authorized_contract_does_not_relax_original_real_gates():
    raw = Path(CONTRACT).read_bytes()
    assert sha(raw) == CONTRACT_SHA
    c = json.loads(raw)
    assert c["execution_authorized"] is True and c["endpoints_bohr"] == [60, 90, 120]
    old = json.loads(
        Path(CONTRACT).with_name("stretched-real-guard-contract-20261101.json").read_text()
    )
    for key in [
        "real_s_absolute_budget",
        "real_phase_sum_budget_rad",
        "real_invariant_absolute_budget",
        "real_flux_absolute_budget",
    ]:
        assert c[key] == old[key]
    metrics = dict(
        s_absolute_error=0.0,
        phase_sum_error_mod_pi_rad=0.0,
        flux_defect=2e-7,
        wronskian_residual=0.0,
        unitarity_defect=0.0,
        reciprocity_defect=0.0,
    )
    result = checks(metrics, c)
    assert result["flux_defect"] is False and result["s_absolute_error"] is True
