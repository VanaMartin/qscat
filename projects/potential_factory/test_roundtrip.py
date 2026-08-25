# projects/potential_factory/test_roundtrip.py
from __future__ import annotations

import numpy as np
import pytest
from qscat.model import F2, N2, NO

from projects.potential_factory.ansatz import from_diatomic, with_params
from projects.potential_factory.extract import extract_target
from projects.potential_factory.fit import fit
from projects.potential_factory.report import FitReport
from projects.potential_factory.tracker import ElectronicPair

# Each `(R_hi, R_lo)` node window below has to leave enough RESONANT nodes in
# the extracted target. `lam` and `alpha` are nearly non-identifiable from a
# single pole: at F2's working point the two Hellmann-Feynman derivatives are
# almost parallel in the complex plane (`dE/dalpha ~ -6.1 dE/dlam`), so the
# per-node 2x2 Newton Jacobian is close to singular and a compensating shift in
# the pair is nearly invisible in the pole. A BOUND node cannot break that -- it
# contributes one real equation and has `alpha` frozen -- so only resonant nodes
# constrain both, and since the ratio above varies with `R`, it takes several of
# them at different widths. F2's (4.0, 2.0) window left just 2 resonant nodes
# and the fit settled into a uniform ~1.3% lam/alpha compensation (E_res rms
# 1.8e-3 Ha, over tolerance).
#
# The windows no longer have to make the T1 fit grid COINCIDE with the target's
# own extraction nodes: `fit_resonance` reads the target ON those nodes
# (`_t1_nodes`). It used to rebuild an independent `linspace(*R_range,
# n_nodes)`, and whenever `extract_target` gated a node out the two grids
# shifted apart, so the fit read the target by interpolation across the
# resonance threshold (F2: `R_c = 2.595`), where `E_res(R)` has a branch point
# and no polynomial interpolant is accurate. Measured on F2 at (3.6, 1.8), the
# window used below, which drops one node: with the independent grid T1 was NOT
# MET (E_res rms 4.6e-3 Ha against a 1e-3 tolerance, Gamma rel max 0.51, `lam`
# 65% off, T3 never attempted, 12.1 s); on the target's own nodes it is met at
# E_res rms 1.0e-11 Ha, Gamma rel max 0.000, `lam`/`alpha` recovered to ~1e-9,
# T3 rms 7.6e-10, in 3.8 s. `fit_neutral` (T0) and `_coupling_eval_grid` (T3)
# read their targets on-node for the same reason.
CASES = {"N2": (N2, (3.0, 1.6)), "NO": (NO, (3.2, 1.7)), "F2": (F2, (3.6, 1.8))}

# Sign of the dissociative-attachment threshold `(-ea) - eps_0` measured from
# the fitted model's own v=0 level: `+1` endothermic, `-1` exothermic.
# Independently recomputed thresholds: N2 +0.5016 Ha, NO +0.1719 Ha (the value
# CLAUDE.md records for NO's DA channel), F2 -0.0691 Ha.
DA_SIGN = {"N2": +1, "NO": +1, "F2": -1}


@pytest.mark.slow
@pytest.mark.parametrize("name", list(CASES))
def test_round_trip_recovers_the_published_model(name, tmp_path):
    model, (R_hi, R_lo) = CASES[name]
    pair = ElectronicPair()
    R_desc = np.linspace(R_hi, R_lo, 10)
    target = extract_target(model, pair=pair, R_desc=R_desc, n_eps=4, name=name)
    seed = with_params(
        from_diatomic(model),
        {
            "D_e": model.D0 * 1.2,
            "R_e": model.R0 * 1.05,
            "beta0": model.alpha0 * 0.9,
            "lam.f_inf": model.lambda_inf * 1.1,
            "alpha.f_inf": model.alpha_c * 1.3,
        },
    )
    fitted, report = fit(target, seed, pair=pair, n_nodes=10)
    assert [t.status for t in report.tiers] == ["met", "met", "met"], [
        t.detail for t in report.tiers
    ]
    assert abs(fitted.D_e - model.D0) < 1e-8 and abs(fitted.R_e - model.R0) < 1e-8
    np.testing.assert_allclose(fitted.lam_R(R_desc).real, model.lam(R_desc).real, rtol=3e-3)
    np.testing.assert_allclose(fitted.alpha_R(R_desc).real, model.alpha_c, rtol=3e-3)
    assert report.crossing_R is not None and abs(report.crossing_R - model.R_c) < 0.1
    assert report.da_threshold_sign == DA_SIGN[name]
    # The ECS block separates the ansatz's ANALYTIC bounds from the angles
    # actually evaluated; only `tail_growth` is a measurement.
    ecs = report.ecs_bounds_deg
    assert set(ecs) == {
        "electronic_max_deg",
        "nuclear_max_deg",
        "probed_electronic_deg",
        "probed_nuclear_deg",
        "tail_growth",
    }
    assert ecs["electronic_max_deg"] == 45.0 and ecs["nuclear_max_deg"] == 90.0
    assert ecs["probed_electronic_deg"] == pair.grid_a.spec.elements[-1].angle_deg
    assert ecs["probed_electronic_deg"] < ecs["electronic_max_deg"]
    assert ecs["probed_nuclear_deg"] < ecs["nuclear_max_deg"]
    assert ecs["tail_growth"] <= 10.0
    report.to_json(tmp_path / "r.json")
    back = FitReport.from_json(tmp_path / "r.json")
    assert back.parameters == report.parameters and back.tiers == report.tiers
