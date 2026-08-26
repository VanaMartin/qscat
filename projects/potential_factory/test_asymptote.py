"""The `asymptote` tier: the anion curve must END at the atom + ion energy
through its declared long-range tail, and the neutral at the free-atom
limit -- verified on the model itself, not read off a table."""

from __future__ import annotations

import numpy as np
from qscat.core.dissociation import anion_electronic_states
from qscat.model import N2

from projects.potential_factory.ansatz import from_diatomic, with_params
from projects.potential_factory.fit import _asymptotic_nodes, check_asymptote
from projects.potential_factory.report import Tolerances
from projects.potential_factory.target import Curve, ResonanceTarget, polarisation_tail
from projects.potential_factory.tracker import ElectronicPair


def _n2_resonance_target(pair: ElectronicPair, R_inf: float = 10.0) -> ResonanceTarget:
    eps_e, _ = anion_electronic_states(pair.grid_a, N2, R_inf, 1)
    ea = -(float(eps_e[0]) - float(N2.v0(R_inf).real))
    R = np.linspace(1.6, 3.0, 8)
    flat = Curve.from_table(R, np.zeros_like(R))
    return ResonanceTarget(flat, flat, ea, (1.6, 3.0), R_inf=R_inf)


def test_polarisation_tail_is_minus_alpha_over_2r4():
    tail = polarisation_tail(5.3)
    assert abs(float(tail(2.0)) - (-5.3 / 32.0)) < 1e-15
    t = ResonanceTarget(
        Curve.from_callable(lambda R: R), Curve.from_callable(lambda R: R), 0.05, (1, 2), tail=tail
    )
    want = [-0.05 - 5.3 / 32, -0.05 - 5.3 / 512]
    np.testing.assert_allclose(t.v_ion_asymptotic([2.0, 4.0]), want)
    t0 = ResonanceTarget(t.v_ion, t.gamma, 0.05, (1, 2))
    np.testing.assert_allclose(t0.v_ion_asymptotic([2.0, 4.0]), [-0.05, -0.05])


def test_asymptotic_nodes_run_from_r_inf_to_far_beyond_it():
    t = ResonanceTarget(
        Curve.from_callable(lambda R: R), Curve.from_callable(lambda R: R), 0.05, (1, 6), R_inf=12.0
    )
    nodes = _asymptotic_nodes(t)
    assert abs(nodes[0] - 12.0) < 1e-12 and abs(nodes[-1] - 60.0) < 1e-12
    assert nodes.size >= 3 and np.all(np.diff(nodes) > 0)


def test_published_n2_meets_its_own_asymptote():
    pair = ElectronicPair()
    t = _n2_resonance_target(pair)
    res = check_asymptote(t, from_diatomic(N2), pair, Tolerances())
    assert res.status == "met", res.detail
    # Not zero: `ea` was read at R_inf = 10, where N2's Morse neutral is still
    # -1.5e-4 Ha below its limit and lam(R) is 3e-4 short of lam_inf -- the
    # tier measures the TOTAL V_ion against -EA, and both offsets show in it.
    assert res.rms < Tolerances().e_res_rms and res.name == "asymptote"


def test_deeper_well_at_large_r_misses_the_asymptote():
    pair = ElectronicPair()
    t = _n2_resonance_target(pair)
    m = from_diatomic(N2)
    m = with_params(m, {"lam.f_inf": m.lam.f_inf * 1.1})  # anion 10 % too deep everywhere
    res = check_asymptote(t, m, pair, Tolerances())
    assert res.status == "not met"
    assert res.rms > Tolerances().e_res_rms and "R_inf=10" in res.detail
