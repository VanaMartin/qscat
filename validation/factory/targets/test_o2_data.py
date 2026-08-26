"""The extracted O2 tables must be the curves of Alt & Houfek (2021) Fig. 2."""

from __future__ import annotations

import numpy as np
import pytest

from validation.factory.targets.o2 import o2_seed, o2_target
from validation.factory.targets.o2_data import EV, PRECISION_HA, load_o2

EA_O_EV = 1.4611  # Table I (expt.), p. 032829-3
D0_CALC_EV = 5.159  # Table I (calc.), p. 032829-3
ZPE_EV = 0.098  # omega_e/2 of O2, Huber & Herzberg


def test_tables_are_on_one_strictly_increasing_grid():
    c = load_o2()
    assert np.all(np.diff(c.R) > 0)
    assert c.R.min() < 1.85 and c.R.max() > 5.9
    assert np.all(np.isfinite(c.v0)) and np.all(np.isfinite(c.v_ion)) and np.all(c.gamma >= 0)


def test_neutral_curve_matches_table_i_well_depth_and_asymptote():
    c = load_o2()
    depth = -c.v0.min() * EV
    assert abs(depth - (D0_CALC_EV + ZPE_EV)) < 3 * PRECISION_HA * EV
    assert abs(c.v0[-1] * EV) < 3 * PRECISION_HA * EV
    assert 2.2 < c.R[np.argmin(c.v0)] < 2.3


def test_crossing_and_width_support():
    c = load_o2()
    # the resonance branch meets the bound anion curve at the crossing
    assert 2.27 < c.R_c < 2.31
    d = c.v_ion - c.v0
    assert d[c.R < 2.2].max() > 0 and d[c.R > 2.4].max() < 0
    # Gamma > 0 only on the resonant side, and not beyond the figure's last width point
    assert c.gamma[c.R < 2.0].min() > 0.1 / EV  # Gamma(2.0) ~ 0.17 eV, rising to 0.67 at 1.8
    assert np.all(c.gamma[c.R > 2.45] == 0.0)


def test_anion_curve_is_bound_below_the_neutral_at_its_minimum():
    c = load_o2()
    j = np.argmin(c.v_ion)
    assert 2.4 < c.R[j] < 2.6
    assert (c.v0.min() - c.v_ion.min()) * EV > 0.4  # adiabatic EA(O2) ~ 0.45 eV (Table I)


def test_target_is_complete_and_the_seed_is_a_model():
    from qscat.model import ResonanceModel

    t = o2_target()
    assert t.ell == 2 and t.coordinates == ("R",)
    assert t.neutral is not None and t.resonance is not None and t.coupling is not None
    assert t.coupling.alpha_exponent == 2.5
    assert abs(t.resonance.ea * EV - EA_O_EV) < 1e-6
    assert isinstance(o2_seed(), ResonanceModel)


def test_extractor_reproduces_the_committed_tables_when_the_pdf_is_present(tmp_path):
    pymupdf = pytest.importorskip("pymupdf")
    from validation.factory.extract_fig2 import PDF, extract

    if not PDF.exists():
        pytest.skip("paper PDF not present (gitignored)")
    del pymupdf
    fresh = extract(out=tmp_path)
    for name, arr in fresh.items():
        committed = np.loadtxt(
            f"validation/factory/data/o2/{name}.csv", delimiter=",", comments="#"
        )
        np.testing.assert_allclose(arr, committed, rtol=0, atol=1e-9)
