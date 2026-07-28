from __future__ import annotations

import pytest
from qscat.model import N2
from qscat.tuning import IncidentSpec, propose_grid, required_extent, tw_analysis
from qscat.tuning.incident import _interaction_extent


def test_required_extent_grows_past_the_wavepacket_and_its_tail():
    spec = IncidentSpec(position=45.0, impulse=-0.5, sigma=4.0)
    assert spec.required_extent() > 45.0


def test_required_extent_free_function_matches_the_method():
    spec = IncidentSpec(position=45.0, impulse=-0.5, sigma=4.0, observation=50.0)
    assert required_extent(spec) == spec.required_extent()


def test_required_extent_respects_a_far_observation_boundary():
    # observation beyond the wavepacket's own tail must win.
    spec = IncidentSpec(position=10.0, impulse=-0.5, sigma=1.0, observation=100.0)
    assert spec.required_extent() == 100.0


def test_required_extent_defaults_observation_to_zero():
    near = IncidentSpec(position=1.0, impulse=-0.5, sigma=0.5)
    assert near.required_extent() == pytest.approx(1.0 + 5.0 * 0.5)


def test_propose_grid_nuclear_extends_past_a_far_incident_position():
    # Nuclear real cutoffs are ~14-22 bohr by default; an incident placed at
    # 45 bohr must force the real region to extend past it.
    incident = IncidentSpec(position=45.0, impulse=-0.5, sigma=4.0)
    g = propose_grid(N2, "nuclear", (0.04, 0.18), incident=incident)
    assert g.R0 > 45.0


def test_propose_grid_electronic_extends_past_a_far_incident_position():
    incident = IncidentSpec(position=45.0, impulse=-0.5, sigma=4.0)
    g = propose_grid(N2, "electronic", (0.04, 0.18), incident=incident)
    assert g.R0 > 45.0


def test_propose_grid_without_incident_does_not_reach_45_bohr():
    # Sanity check that the extension above is actually incident's doing,
    # not just the default electronic/nuclear cutoff already being > 45.
    g = propose_grid(N2, "electronic", (0.04, 0.18))
    assert g.R0 < 45.0


def test_interaction_extent_is_positive_and_finite_for_n2():
    r_int = _interaction_extent(N2)
    assert 0.0 < r_int < 30.0


def test_tw_analysis_rejects_empty_energy_range():
    with pytest.raises(ValueError):
        tw_analysis(N2, (0.18, 0.04))


def test_tw_analysis_spectrum_brackets_the_energy_range():
    e_min, e_max = 0.04, 0.18
    spec = tw_analysis(N2, (e_min, e_max))

    # Mean energy sits inside (or very close to) the requested range.
    e_centre = 0.5 * (spec.impulse**2)
    assert e_min <= e_centre <= e_max

    # Inward-launched, per the n2_2d_td_cross_section convention.
    assert spec.impulse < 0.0

    # The realized energy spread covers at least half the range on either
    # side of the centre (the bracketing condition `tw_analysis` targets).
    delta_e = abs(spec.impulse) / (2.0 * spec.sigma)
    assert delta_e >= 0.5 * (e_max - e_min)

    # Placement is sane: positive, with an observation boundary at/beyond it.
    assert spec.position > 0.0
    assert spec.observation is not None
    assert spec.observation >= spec.position
