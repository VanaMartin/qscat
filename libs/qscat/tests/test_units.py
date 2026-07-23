import numpy as np
import pytest
import qscat
from qscat import units


def test_version_is_string() -> None:
    assert isinstance(qscat.__version__, str)


def test_hartree_to_ev_known_value() -> None:
    # CODATA 2018: 1 Hartree = 27.211386245988 eV
    assert units.hartree_to_ev(1.0) == pytest.approx(27.211386245988, rel=1e-12)


def test_roundtrip_is_identity() -> None:
    x = np.array([0.0, 0.5, 2.5, -1.3])
    np.testing.assert_allclose(units.ev_to_hartree(units.hartree_to_ev(x)), x, rtol=1e-14)
