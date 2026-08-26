"""The committed O2 fit report is a gate, not a record: the tiers it claims
must be the ones the note quotes, and the model it describes must rebuild."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from projects.potential_factory.ansatz import TailR, params
from projects.potential_factory.report import FitReport
from validation.factory.targets.o2 import O2_R_INF, o2_model_from_report

REPORT = Path(__file__).parent / "results" / "o2-fit-report.json"


def test_committed_o2_report_meets_t0_t1_and_the_asymptote():
    rep = FitReport.from_json(REPORT)
    status = {t.name: t.status for t in rep.tiers}
    assert status["T0"] == "met" and status["T1"] == "met" and status["asymptote"] == "met"
    assert status["T3"] == "not met"  # the discrete-state inconsistency, documented
    t1 = next(t for t in rep.tiers if t.name == "T1")
    assert t1.rms < 1e-3  # E_res rms under the 20 meV extraction floor
    assert abs(rep.crossing_R - 2.289) < 2e-3 and rep.da_threshold_sign == 1


def test_committed_o2_model_rebuilds_from_its_report():
    rep = FitReport.from_json(REPORT)
    m = o2_model_from_report(rep)
    assert isinstance(m.lam, TailR) and len(m.lam.coeffs) == 9
    back = params(m)
    assert set(back) == set(rep.parameters)
    for k, v in rep.parameters.items():
        assert back[k] == v, k
    # the asymptote the report claims: lam settles to f_inf as R^-4
    R = np.array([O2_R_INF, 5 * O2_R_INF])
    assert abs(float(m.lam(R[1]).real) - m.lam.f_inf) < 1e-4
