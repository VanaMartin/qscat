"""The committed O2 fit report is a gate, not a record: the tiers it claims
must be the ones the note quotes, and the model it describes must rebuild."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from qscat.model import O2
from qscat.model.flexible import TailR, params

from projects.potential_factory.report import FitReport
from validation.factory.targets.o2 import O2_MU, O2_R_INF, o2_model_from_report

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


def test_registry_o2_is_the_committed_report_verbatim():
    """`qscat.model.O2` carries the report's constants key for key -- the
    library cannot read `validation/`, so the lock lives here."""
    rep = FitReport.from_json(REPORT)
    assert params(O2) == rep.parameters
    assert O2.mu == O2_MU and O2.ell == 2 and O2.charge == 0
    m = o2_model_from_report(rep)
    r = np.linspace(0.3, 12.0, 40)
    for R in (1.9, 2.289, 3.0, 6.0, 14.0):
        np.testing.assert_allclose(O2.surface(r, R), m.surface(r, R), rtol=0, atol=1e-15)
