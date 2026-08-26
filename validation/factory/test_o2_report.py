"""The committed O2 fit report is a gate, not a record: the tiers it claims
must be the ones the note quotes, and the model it describes must rebuild."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from qscat.core.dissociation import anion_electronic_states
from qscat.model import O2, O2_SO12, O2_SO32
from qscat.model.flexible import TailR, params

from projects.potential_factory.report import FitReport
from projects.potential_factory.tracker import ElectronicPair
from validation.factory.targets.o2 import O2_MU, O2_R_INF, o2_model_from_report
from validation.factory.targets.o2_data import load_so_split

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


@pytest.mark.parametrize(
    "name,model,so",
    [("o2-so12", O2_SO12, -1), ("o2-so32", O2_SO32, +1)],
)
def test_spin_orbit_components_are_their_reports_and_stay_met(name, model, so):
    rep = FitReport.from_json(REPORT.parent / f"{name}-fit-report.json")
    status = {t.name: t.status for t in rep.tiers}
    assert status == {"T1": "met", "asymptote": "met"}
    assert params(model) == rep.parameters
    # only the well moved: the neutral curve is the parent's, exactly
    assert model.D_e == O2.D_e and model.betas == O2.betas and model.R_e == O2.R_e
    # and it moved the right way, by the right amount: on the bound branch
    # the anion's electronic energy sits so * Delta_SO(R)/2 from the parent's
    # (Fig. 1: 18 meV at 3 bohr -> +-9 meV = +-3.3e-4 Ha)
    g = ElectronicPair().grid_a
    e_par, _ = anion_electronic_states(g, O2, 3.0, 1)
    e_so, _ = anion_electronic_states(g, model, 3.0, 1)
    R_so, d_so = load_so_split()
    expect = so * 0.5 * float(np.interp(3.0, R_so, d_so))
    assert abs((e_so[0] - e_par[0]) - expect) < 0.3 * abs(expect)


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
