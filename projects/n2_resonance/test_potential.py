"""Failing-first tests for the N2 electronic potential (Task 1, sub-project #2).

Cross-checks against `qscat.model.N2` (the layer-legal, independently
implemented and independently gated -- libs/qscat/tests/test_model.py --
`DiatomicResonanceModel` instance locked to the same eMoScat deck constants),
and also checks the concrete numeric assertions from the task brief.
"""

import numpy as np
from qscat.model import N2

from projects.n2_resonance import potential


def test_lambda_at_Rc_matches_config():
    p = potential.PARAMS["potential"]
    assert potential.lam(p["R_c"]) == p["lambda_c"]


def test_v0_at_R0_is_minus_D0():
    p = potential.PARAMS["potential"]
    assert potential.v0(p["R_0"]) == -p["D_0"]


def test_v_int_negative_well():
    R0 = potential.PARAMS["potential"]["R_0"]
    assert potential.v_int(1.0, R0) < 0.0


def test_v_eff_el_matches_explicit_centrifugal_formula():
    R0 = potential.PARAMS["potential"]["R_0"]
    r = 2.0
    expected = potential.v_int(r, R0) + 2 * 3 / (2 * 2**2)
    assert potential.v_eff_el(r, R0) == expected


def test_matches_library_model_to_1e_12():
    """Cross-check against qscat.model.N2 -- the same eMoScat constants,
    independently implemented (DiatomicResonanceModel) and independently
    gated (libs/qscat/tests/test_model.py). potential.v_eff_el excludes
    v0(R) while N2.surface includes it, hence the subtraction."""
    rs = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 10.0])
    Rs = np.array([1.5, 2.01943, 2.405, 3.0, 4.0])
    for R in Rs:
        assert abs(potential.v0(R) - complex(N2.v0(R))) < 1e-12
        assert abs(potential.lam(R) - complex(N2.lam(R))) < 1e-12
        for r in rs:
            assert abs(potential.v_int(r, R) - complex(N2.v_int(r, R))) < 1e-12
            assert (
                abs(potential.v_eff_el(r, R) - (complex(N2.surface(r, R)) - complex(N2.v0(R))))
                < 1e-12
            )


def test_vectorized_over_r_array():
    R0 = potential.PARAMS["potential"]["R_0"]
    r = np.array([1.0, 2.0, 3.0])
    out = potential.v_eff_el(r, R0)
    assert out.shape == r.shape
